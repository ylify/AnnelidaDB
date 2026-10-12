import os
import re
import glob
import subprocess
from django.core.management.base import BaseCommand
from django.db import transaction
from api.models import Taxonomy, SequenceFileIndex


EXCLUDE_PATTERN = re.compile(r'longest', re.IGNORECASE)


def _format_species_dir(raw_species):
    if not raw_species:
        return ''
    formatted = raw_species.replace('sp.', 'sp_')
    formatted = formatted.replace(' ', '_')
    formatted = re.sub(r'__+', '_', formatted)
    return formatted


class Command(BaseCommand):
    help = 'Build indexed mappings for protein and CDS FASTA files'

    def add_arguments(self, parser):
        parser.add_argument(
            '--data-dir',
            type=str,
            default='static/data',
            help='Base data directory containing species subdirectories (default: static/data)',
        )
        parser.add_argument(
            '--species',
            type=str,
            default=None,
            help='Only process a single species directory name',
        )
        parser.add_argument(
            '--force',
            action='store_true',
            help='Rebuild the sequence index for every file even if entries already exist',
        )

    def handle(self, *args, **options):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', options['data_dir']))
        if not os.path.isdir(base_dir):
            self.stderr.write(f"Data directory does not exist: {base_dir}")
            return

        taxonomies = Taxonomy.objects.all()
        if options['species']:
            taxonomies = Taxonomy.objects.filter(species__iexact=options['species'])
            if not taxonomies.exists():
                self.stdout.write(self.style.WARNING(f'Species not found in taxonomy table: "{options["species"]}"'))
                return

        for tax in taxonomies.order_by('species'):
            formatted_species = _format_species_dir(tax.species)
            species_path = os.path.join(base_dir, formatted_species)
            genome_dir = os.path.join(species_path, 'genome')
            if not os.path.isdir(genome_dir):
                self.stdout.write(self.style.WARNING(f'Skipping {tax.species}: no genome directory found at {species_path}'))
                continue

            self.stdout.write(self.style.NOTICE(f'Processing species: {tax.species} ({formatted_species})'))
            fasta_files = self._find_sequence_files(genome_dir)
            if not fasta_files:
                self.stdout.write(self.style.WARNING(f'No protein/CDS FASTA files found for {formatted_species}'))
                continue

            for fasta_path in fasta_files:
                self._process_fasta_file(tax, fasta_path, options['force'])

    def _find_sequence_files(self, genome_dir):
        patterns = [
            '*.pep.fasta.gz', '*.pep.fa.gz', '*.protein.fasta.gz',
            '*.pep.fasta', '*.pep.fa', '*.protein.fasta',
            '*.CDS.fasta.gz', '*.cds.fasta.gz',
            '*.CDS.fasta', '*.cds.fasta',
        ]
        files = []
        for pattern in patterns:
            for path in sorted(glob.glob(os.path.join(genome_dir, pattern))):
                filename = os.path.basename(path)
                if EXCLUDE_PATTERN.search(filename):
                    continue
                if path in files:
                    continue
                files.append(path)
        return files

    def _infer_seq_type(self, path):
        file_name = os.path.basename(path).lower()
        if '.pep.' in file_name or '.protein.' in file_name:
            return 'protein'
        if '.cds.' in file_name:
            return 'cds'
        return None

    def _is_bgzf(self, path):
        result = subprocess.run(['bgzip', '-t', path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return result.returncode == 0

    def _bgzip_path(self, path):
        if path.lower().endswith('.bgz'):
            return path
        return path + '.bgz'

    def _ensure_bgzf(self, path):
        if path.lower().endswith('.bgz'):
            return path
        if self._is_bgzf(path):
            return path

        bgz_path = self._bgzip_path(path)
        if os.path.exists(bgz_path) and os.path.getmtime(bgz_path) >= os.path.getmtime(path):
            return bgz_path

        self.stdout.write(self.style.NOTICE(f'Converting to BGZF: {path} -> {bgz_path}'))
        with open(bgz_path, 'wb') as out_file:
            gunzip = subprocess.Popen(['gunzip', '-c', path], stdout=subprocess.PIPE)
            bgzip = subprocess.Popen(['bgzip', '-c'], stdin=gunzip.stdout, stdout=out_file)
            gunzip.stdout.close()
            retcode = bgzip.wait()
            gunzip.wait()
            if retcode != 0:
                raise RuntimeError(f'bgzip failed for {path}')
        return bgz_path

    def _build_faidx(self, fasta_path):
        fai_path = fasta_path + '.fai'
        if os.path.exists(fai_path) and os.path.getmtime(fai_path) >= os.path.getmtime(fasta_path):
            self.stdout.write(self.style.NOTICE(f'Index already exists: {fai_path}'))
            return

        self.stdout.write(self.style.NOTICE(f'Running samtools faidx for {fasta_path}'))
        result = subprocess.run(['samtools', 'faidx', fasta_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f'samtools faidx failed for {fasta_path}: {result.stderr.strip()}')

    def _process_fasta_file(self, tax, fasta_path, force):
        seq_type = self._infer_seq_type(fasta_path)
        if not seq_type:
            self.stdout.write(self.style.WARNING(f'Skipping unsupported file type: {fasta_path}'))
            return

        actual_path = self._ensure_bgzf(fasta_path) if fasta_path.lower().endswith('.gz') else fasta_path
        self._build_faidx(actual_path)

        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
        relative_path = os.path.relpath(actual_path, project_root)
        file_name = os.path.basename(fasta_path)

        with transaction.atomic():
            if force:
                SequenceFileIndex.objects.filter(taxonomy=tax, seq_type=seq_type, file_name=file_name).delete()

            if SequenceFileIndex.objects.filter(taxonomy=tax, seq_type=seq_type, file_name=file_name).exists():
                self.stdout.write(self.style.NOTICE(f'Index entries already exist for {file_name}, skipping'))
                return

            existing_seq_ids = set(SequenceFileIndex.objects.filter(taxonomy=tax, seq_type=seq_type).values_list('seq_id', flat=True))
            self.stdout.write(self.style.NOTICE(f'Loading FASTA index entries from {actual_path}'))
            fai_path = actual_path + '.fai'
            entries = []
            seen_seq_ids = set()
            with open(fai_path, 'r', encoding='utf-8') as fai_handle:
                for line in fai_handle:
                    if not line.strip():
                        continue
                    seq_id = line.split('\t', 1)[0]
                    if seq_id in existing_seq_ids or seq_id in seen_seq_ids:
                        continue
                    seen_seq_ids.add(seq_id)
                    entries.append(SequenceFileIndex(
                        taxonomy=tax,
                        seq_type=seq_type,
                        seq_id=seq_id,
                        file_name=file_name,
                        file_path=relative_path,
                    ))
                    if len(entries) >= 5000:
                        SequenceFileIndex.objects.bulk_create(entries)
                        entries = []
            if entries:
                SequenceFileIndex.objects.bulk_create(entries)

            self.stdout.write(self.style.SUCCESS(f'Indexed {file_name} ({seq_type}) for species {tax.species}'))

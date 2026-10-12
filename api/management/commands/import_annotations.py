import os
import re
import csv
import gzip
import itertools
from django.core.management.base import BaseCommand
from api.models import Taxonomy, EggnogAnnotation, InterproScanAnnotation


def _format_species_dir(raw_species):
    if not raw_species:
        return ''
    formatted = raw_species.replace('sp.', 'sp_')
    formatted = formatted.replace(' ', '_')
    formatted = re.sub(r'__+', '_', formatted)
    return formatted


def _choose_annotation_file(files):
    if not files:
        return None
    annevo_files = [f for f in files if '.annevo.' in f.lower()]
    if annevo_files:
        return sorted(annevo_files)[0]
    return sorted(files)[0]


class Command(BaseCommand):
    help = 'Import EggNOG and InterProScan annotations from species data directories'

    def add_arguments(self, parser):
        parser.add_argument(
            '--data-dir',
            type=str,
            default='static/data',
            help='Directory containing species data (default: static/data)',
        )

    def handle(self, *args, **options):
        data_dir = options['data_dir']
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        full_data_dir = os.path.join(base_dir, data_dir)

        if not os.path.exists(full_data_dir):
            self.stdout.write(self.style.ERROR(f'Data directory does not exist: {full_data_dir}'))
            return

        for taxonomy in Taxonomy.objects.all().order_by('species'):
            species_dir = _format_species_dir(taxonomy.species)
            species_path = os.path.join(full_data_dir, species_dir)
            if not os.path.isdir(species_path):
                continue

            # Only process species that have genome records
            if not taxonomy.genomes.exists():
                self.stdout.write(self.style.WARNING(f'Skipping {taxonomy.species} - no genome records found'))
                continue

            # Find annotation files
            genome_dir = os.path.join(species_path, 'genome')
            if not os.path.exists(genome_dir):
                continue

            eggnog_files = []
            interpro_files = []
            for file in os.listdir(genome_dir):
                lower_name = file.lower()
                if lower_name.endswith('.emapper.annotations') or lower_name.endswith('.emapper.annotations.gz'):
                    eggnog_files.append(file)
                if lower_name.endswith('.pep.fasta.tsv') or lower_name.endswith('.pep.fasta.tsv.gz'):
                    interpro_files.append(file)

            if not eggnog_files and not interpro_files:
                self.stdout.write(self.style.WARNING(f'No annotation files found for {taxonomy.species}'))
                continue

            if eggnog_files:
                EggnogAnnotation.objects.filter(taxonomy=taxonomy).delete()
                for eggnog_file in sorted(eggnog_files):
                    eggnog_path = os.path.join(genome_dir, eggnog_file)
                    self.stdout.write(f'Importing EggNOG annotations for {taxonomy.species} from {eggnog_path}')
                    self._import_eggnog_annotations(taxonomy, eggnog_path)

            if interpro_files:
                InterproScanAnnotation.objects.filter(taxonomy=taxonomy).delete()
                for interpro_file in sorted(interpro_files):
                    interpro_path = os.path.join(genome_dir, interpro_file)
                    self.stdout.write(f'Importing InterProScan annotations for {taxonomy.species} from {interpro_path}')
                    self._import_interpro_annotations(taxonomy, interpro_path)

    def _import_eggnog_annotations(self, taxonomy, annotation_file):
        open_func = gzip.open if annotation_file.lower().endswith('.gz') else open
        with open_func(annotation_file, 'rt', encoding='utf-8', errors='ignore') as f:
            data_lines = (line for line in f if not line.startswith('#'))
            reader = csv.DictReader(data_lines, delimiter='	', fieldnames=[
                'query', 'seed_ortholog', 'evalue', 'score', 'eggnog_ogs',
                'max_annot_lvl', 'cog_category', 'description', 'preferred_name',
                'gos', 'ec', 'kegg_ko', 'kegg_pathway', 'kegg_module',
                'kegg_reaction', 'kegg_rclass', 'brite', 'kegg_tc', 'cazys',
                'bigg_reaction', 'pfams'
            ])

            batch = []
            batch_size = 1000
            total = 0
            for row in reader:
                query = row.get('query')
                if not query:
                    continue

                annotation = EggnogAnnotation(
                    taxonomy=taxonomy,
                    query=query,
                    seed_ortholog=row.get('seed_ortholog') or None,
                    evalue=float(row['evalue']) if row.get('evalue') else None,
                    score=float(row['score']) if row.get('score') else None,
                    eggnog_ogs=row.get('eggnog_ogs') or None,
                    max_annot_lvl=row.get('max_annot_lvl') or None,
                    cog_category=row.get('cog_category') or None,
                    description=row.get('description') or None,
                    preferred_name=row.get('preferred_name') or None,
                    gos=row.get('gos') or None,
                    ec=row.get('ec') or None,
                    kegg_ko=row.get('kegg_ko') or None,
                    kegg_pathway=row.get('kegg_pathway') or None,
                    kegg_module=row.get('kegg_module') or None,
                    kegg_reaction=row.get('kegg_reaction') or None,
                    kegg_rclass=row.get('kegg_rclass') or None,
                    brite=row.get('brite') or None,
                    kegg_tc=row.get('kegg_tc') or None,
                    cazys=row.get('cazys') or None,
                    bigg_reaction=row.get('bigg_reaction') or None,
                    pfams=row.get('pfams') or None,
                )
                batch.append(annotation)
                if len(batch) >= batch_size:
                    EggnogAnnotation.objects.bulk_create(batch)
                    total += len(batch)
                    batch = []

            if batch:
                EggnogAnnotation.objects.bulk_create(batch)
                total += len(batch)

        self.stdout.write(self.style.SUCCESS(f'Imported {total} EggNOG annotations for {taxonomy.species}'))

    def _import_interpro_annotations(self, taxonomy, annotation_file):
        open_func = gzip.open if annotation_file.lower().endswith('.gz') else open
        with open_func(annotation_file, 'rt', encoding='utf-8', errors='ignore') as f:
            def non_comment_lines():
                for line in f:
                    if line.startswith('#'):
                        continue
                    yield line

            line_iter = non_comment_lines()
            first_line = next(line_iter, None)
            if first_line is None:
                return

            if first_line.strip().lower().startswith('protein id'):
                data_iter = line_iter
            else:
                data_iter = itertools.chain([first_line], line_iter)

            reader = csv.DictReader(data_iter, delimiter='	', fieldnames=[
                'Protein ID', 'MD5', 'Length', 'Analysis', 'Signature Accession',
                'Signature Description', 'Start', 'End', 'Score', 'Status', 'Date',
                'IPR Accession', 'IPR Description', 'GO', 'Pathway'
            ])

            batch = []
            batch_size = 1000
            total = 0

            for row in reader:
                protein_id = row.get('Protein ID')
                if not protein_id:
                    continue

                def safe_int(value):
                    try:
                        return int(value)
                    except (TypeError, ValueError):
                        return None

                def safe_float(value):
                    try:
                        return float(value)
                    except (TypeError, ValueError):
                        return None

                annotation = InterproScanAnnotation(
                    taxonomy=taxonomy,
                    query=protein_id,
                    md5=row.get('MD5') or None,
                    length=safe_int(row.get('Length')),
                    analysis=row.get('Analysis') or None,
                    signature_accession=row.get('Signature Accession') or None,
                    signature_description=row.get('Signature Description') or None,
                    start=safe_int(row.get('Start')),
                    end=safe_int(row.get('End')),
                    score=safe_float(row.get('Score')),
                    status=row.get('Status') or None,
                    date=row.get('Date') or None,
                    ipr_accession=row.get('IPR Accession') or None,
                    ipr_description=row.get('IPR Description') or None,
                    gos=row.get('GO') or None,
                    pathway=row.get('Pathway') or None,
                )
                batch.append(annotation)
                if len(batch) >= batch_size:
                    InterproScanAnnotation.objects.bulk_create(batch)
                    total += len(batch)
                    batch = []

            if batch:
                InterproScanAnnotation.objects.bulk_create(batch)
                total += len(batch)

        self.stdout.write(self.style.SUCCESS(f'Imported {total} InterProScan annotations for {taxonomy.species}'))

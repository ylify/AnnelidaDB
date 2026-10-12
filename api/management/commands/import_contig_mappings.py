import os
import re
from django.core.management.base import BaseCommand
from api.models import Taxonomy, ContigMapping


def _format_species_dir(raw_species):
    if not raw_species:
        return ''
    formatted = raw_species.replace('sp.', 'sp_')
    formatted = formatted.replace(' ', '_')
    formatted = re.sub(r'__+', '_', formatted)
    return formatted


class Command(BaseCommand):
    help = 'Import contig mappings from FASTA files'

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

            genome_dir = os.path.join(species_path, 'genome')
            if not os.path.exists(genome_dir):
                continue

            # Process genomic.fasta.masked
            genomic_file = None
            for file in os.listdir(genome_dir):
                if file.endswith('.genomic.fasta.masked.gz'):
                    genomic_file = os.path.join(genome_dir, file)
                    break

            if genomic_file:
                self.stdout.write(f'Processing genomic file for {taxonomy.species}: {genomic_file}')
                self._process_fasta_file(genomic_file, taxonomy, taxonomy.species)

            # Process CDS.fasta and compressed CDS files
            cds_file = None
            for file in os.listdir(genome_dir):
                if file.endswith('.CDS.fasta') or file.endswith('.CDS.fasta.gz'):
                    cds_file = os.path.join(genome_dir, file)
                    break

            if cds_file:
                self.stdout.write(f'Processing CDS file for {taxonomy.species}: {cds_file}')
                self._process_fasta_file(cds_file, taxonomy, taxonomy.species)

            # Process protein/peptide files, including compressed versions
            protein_files = []
            for file in os.listdir(genome_dir):
                if file.endswith('.pep.fasta') or file.endswith('.pep.fasta.gz') or file.endswith('.protein.fasta') or file.endswith('.protein.fasta.gz') or 'pep' in file.lower():
                    protein_files.append(os.path.join(genome_dir, file))

            for protein_file in protein_files:
                self.stdout.write(f'Processing protein file for {taxonomy.species}: {protein_file}')
                self._process_fasta_file(protein_file, taxonomy, taxonomy.species)

    def _process_fasta_file(self, file_path, taxonomy, species_name):
        import gzip

        open_func = gzip.open if file_path.endswith('.gz') else open
        try:
            with open_func(file_path, 'rt', encoding='utf-8', errors='ignore') as fasta_file:
                count = 0
                for line in fasta_file:
                    line = line.strip()
                    if not line or not line.startswith('>'):
                        continue

                    contig_id = line[1:].split()[0]
                    _, created = ContigMapping.objects.get_or_create(
                        contig_id=contig_id,
                        species=species_name,
                        defaults={'taxonomy': taxonomy}
                    )
                    if created:
                        count += 1
        except Exception as exc:
            self.stdout.write(self.style.ERROR(f'Failed to read file: {file_path} ({exc})'))
            return

        self.stdout.write(self.style.SUCCESS(f'Processed {count} contig mappings from {file_path}'))
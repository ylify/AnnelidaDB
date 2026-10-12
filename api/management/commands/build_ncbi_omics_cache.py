import os
from django.core.management.base import BaseCommand
from api.views import _build_ncbi_omics_cache, NCBI_OMICS_CACHE_FILE


class Command(BaseCommand):
    help = 'Build and persist the NCBI omics cache for Annelida species using SRA counts.'

    def add_arguments(self, parser):
        parser.add_argument('--api-key', type=str, default=None, help='Optional NCBI API key for higher rate limits.')
        parser.add_argument('--skip-existing', action='store_true', help='Do not rebuild the cache if it already exists.')
        parser.add_argument('--force', action='store_true', help='Force cache rebuild even if it exists.')
        parser.add_argument('--delay', type=float, default=0.34, help='Seconds to wait between NCBI requests.')

    def handle(self, *args, **options):
        api_key = options['api_key'] or os.environ.get('NCBI_API_KEY')
        delay = max(0.0, options['delay'])
        cache_exists = os.path.isfile(NCBI_OMICS_CACHE_FILE)

        if cache_exists and options['skip_existing'] and not options['force']:
            self.stdout.write(self.style.SUCCESS('NCBI omics cache already exists and --skip-existing was specified.'))
            return

        if cache_exists and not options['force']:
            self.stdout.write('Existing NCBI omics cache will be overwritten.')

        # Attach the API key and delay into the environment for shared helper usage.
        if api_key:
            os.environ['NCBI_API_KEY'] = api_key
        os.environ['NCBI_REQUEST_DELAY'] = str(delay)

        self.stdout.write('Building NCBI omics cache for all Annelida taxonomy entries...')
        data = _build_ncbi_omics_cache()
        self.stdout.write(self.style.SUCCESS(f'Built NCBI omics cache with {len(data.get("counts", {}))} entries.'))
        self.stdout.write(self.style.SUCCESS(f'DNA species count: {data.get("dna_species_count")}, RNA species count: {data.get("rna_species_count")}'))

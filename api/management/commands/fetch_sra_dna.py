import os

from django.core.management.base import BaseCommand, CommandError

from api.views import (
    DNA_SRA_CACHE_FILE,
    _ncbi_sra_dna_query,
    _ncbi_sra_runinfo,
    _save_cached_dna_sra_rows,
)


class Command(BaseCommand):
    help = 'Fetch and cache complete NCBI SRA DNA runinfo for Annelida.'

    def add_arguments(self, parser):
        parser.add_argument('--taxid', type=int, default=6340, help='NCBI root taxonomy ID to search.')
        parser.add_argument('--api-key', type=str, default=None, help='Optional NCBI API key.')
        parser.add_argument('--delay', type=float, default=0.34, help='Seconds to wait between NCBI requests.')
        parser.add_argument('--force', action='store_true', help='Overwrite an existing cache.')

    def handle(self, *args, **options):
        if os.path.isfile(DNA_SRA_CACHE_FILE) and not options['force']:
            self.stdout.write(
                self.style.WARNING(
                    f'DNA SRA cache already exists at {DNA_SRA_CACHE_FILE}; use --force to rebuild.'
                )
            )
            return

        api_key = options['api_key'] or os.environ.get('NCBI_API_KEY')
        if api_key:
            os.environ['NCBI_API_KEY'] = api_key
        os.environ['NCBI_REQUEST_DELAY'] = str(max(0.0, options['delay']))

        term = _ncbi_sra_dna_query(options['taxid'])
        self.stdout.write(f'Fetching DNA SRA runinfo for taxid {options["taxid"]}...')
        rows = _ncbi_sra_runinfo(term, api_key=api_key)
        if not rows:
            raise CommandError('NCBI returned no DNA SRA runinfo; existing cache was left unchanged.')
        fetched_at = _save_cached_dna_sra_rows(rows)
        self.stdout.write(self.style.SUCCESS(
            f'Cached {len(rows)} DNA SRA runinfo rows at {DNA_SRA_CACHE_FILE} ({fetched_at}).'
        ))

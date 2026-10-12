import os
from django.core.management.base import BaseCommand
from api.views import _get_statistics_data


class Command(BaseCommand):
    help = 'Refresh the cached AnnelidaDB statistics summary immediately.'

    def add_arguments(self, parser):
        parser.add_argument(
            '--no-cache',
            action='store_true',
            help='Refresh statistics without updating the cache file.',
        )

    def handle(self, *args, **options):
        stats = _get_statistics_data(force_refresh=True)
        self.stdout.write(self.style.SUCCESS('Statistics refreshed successfully.'))
        self.stdout.write(self.style.SUCCESS(f"taxonomy_count={stats.get('taxonomy_count')}, genome_count={stats.get('genome_count')}"))

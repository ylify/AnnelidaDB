from django.core.management.base import BaseCommand
from django.db import connection

class Command(BaseCommand):
    help = 'Update the Full Text Search virtual table for gene annotations'

    def add_arguments(self, parser):
        parser.add_argument(
            '--rebuild',
            action='store_true',
            help='Rebuild the FTS table from scratch',
        )

    def handle(self, *args, **options):
        with connection.cursor() as cursor:
            if options['rebuild']:
                self.stdout.write('Rebuilding FTS table...')
                # Drop and recreate the FTS table
                cursor.execute('DROP TABLE IF EXISTS gene_search_fts')

            # Create FTS virtual table if it doesn't exist
            cursor.execute('''
                CREATE VIRTUAL TABLE IF NOT EXISTS gene_search_fts USING fts5(
                    description, preferred_name, species, taxonomy_id, gene_id,
                    content='api_eggnogannotation',
                    content_rowid='id'
                )
            ''')

            if options['rebuild']:
                # Populate the FTS table
                cursor.execute('''
                    INSERT INTO gene_search_fts(rowid, description, preferred_name, species, taxonomy_id, gene_id)
                    SELECT ea.id, ea.description, ea.preferred_name, t.species, ea.taxonomy_id, ea.query
                    FROM api_eggnogannotation ea
                    LEFT JOIN api_taxonomy t ON ea.taxonomy_id = t.id
                ''')
                self.stdout.write(self.style.SUCCESS('FTS table rebuilt successfully'))
            else:
                # Update existing entries (this is a simple approach - for production, consider triggers)
                cursor.execute('DELETE FROM gene_search_fts')
                cursor.execute('''
                    INSERT INTO gene_search_fts(rowid, description, preferred_name, species, taxonomy_id, gene_id)
                    SELECT ea.id, ea.description, ea.preferred_name, t.species, ea.taxonomy_id, ea.query
                    FROM api_eggnogannotation ea
                    LEFT JOIN api_taxonomy t ON ea.taxonomy_id = t.id
                ''')
                self.stdout.write(self.style.SUCCESS('FTS table updated successfully'))
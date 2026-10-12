import csv
import json
import os
import time
import urllib.parse
import urllib.request

from django.conf import settings
from django.core.management.base import BaseCommand

from api.models import Taxonomy
from api.views import (
    _fetch_ncbi_taxonomy,
    _get_ncbid_omics_for_taxid,
    _ncbi_sra_search_url,
    _ncbi_sra_dna_query,
    _ncbi_sra_rna_query,
)


class Command(BaseCommand):
    help = 'Generate taxonomy.1.csv by fetching NCBI taxonomy classification and joining GBIF/WoRMS IDs.'

    def add_arguments(self, parser):
        parser.add_argument(
            '--output',
            type=str,
            default='taxonomy.1.csv',
            help='Output CSV filename relative to BASE_DIR.',
        )
        parser.add_argument(
            '--delay',
            type=float,
            default=0.34,
            help='Delay in seconds between NCBI fetch requests.',
        )
        parser.add_argument(
            '--limit',
            type=int,
            default=None,
            help='Limit the number of taxonomy rows to process.',
        )
        parser.add_argument(
            '--source',
            type=str,
            choices=['db', 'csv', 'ncbi'],
            default='ncbi',
            help='Source of NCBI taxids: ncbi for full Annelida subtree, db for Taxonomy table, csv for taxonomy.csv file.',
        )
        parser.add_argument(
            '--taxid-file',
            type=str,
            default='taxonomy.csv',
            help='CSV file with NCBI taxonomy ids when --source=csv.',
        )
        parser.add_argument(
            '--skip-existing',
            action='store_true',
            help='Skip fetching if the output file already exists.',
        )

    def _fetch_taxids_under_annelida(self, max_results=None):
        query = 'Annelida[Subtree] AND species[Rank]'
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
        params = {
            'db': 'taxonomy',
            'term': query,
            'retmode': 'json',
            'retmax': str(max_results or 200000),
        }
        query_url = f"{url}?{urllib.parse.urlencode(params)}"
        with urllib.request.urlopen(query_url, timeout=30) as response:
            data = json.loads(response.read().decode('utf-8'))
        idlist = data.get('esearchresult', {}).get('idlist', [])
        return [str(tid).strip() for tid in idlist if tid]

    def _gbif_name_search(self, name):
        if not name:
            return None
        query = urllib.parse.quote(name)
        url = f'https://api.gbif.org/v1/species/match?name={query}'
        try:
            with urllib.request.urlopen(url, timeout=15) as response:
                data = json.load(response)
        except Exception:
            return None
        gbif_id = data.get('usageKey') or data.get('taxonID')
        return str(gbif_id) if gbif_id else None

    def _worms_name_search(self, name):
        if not name:
            return None
        query = urllib.parse.quote(name)
        url = f'https://www.marinespecies.org/rest/AphiaRecordsByName/{query}?marine_only=false'
        try:
            with urllib.request.urlopen(url, timeout=15) as response:
                data = json.load(response)
        except Exception:
            return None
        if isinstance(data, list) and data:
            aphiaid = data[0].get('AphiaID')
            return str(aphiaid) if aphiaid else None
        return None

    def handle(self, *args, **options):
        output_path = os.path.join(settings.BASE_DIR, options['output'])
        delay = max(0.0, options['delay'])
        limit = options['limit']

        if options['skip_existing'] and os.path.exists(output_path):
            self.stdout.write(self.style.SUCCESS(f'Skipping generation because {output_path} already exists.'))
            return

        if options['source'] == 'csv':
            taxid_file = os.path.join(settings.BASE_DIR, options['taxid_file'])
            if not os.path.exists(taxid_file):
                raise FileNotFoundError(f'CSV taxid file not found: {taxid_file}')
            with open(taxid_file, 'r', encoding='utf-8-sig', newline='') as f:
                reader = csv.DictReader(f)
                taxids = [row.get('NCBI taxonomy id', '').strip() for row in reader if row.get('NCBI taxonomy id', '').strip()]
        elif options['source'] == 'db':
            taxids = list(
                Taxonomy.objects.exclude(ncbi_taxonomy_id__isnull=True)
                .order_by('species')
                .values_list('ncbi_taxonomy_id', flat=True)
            )
        else:
            taxids = self._fetch_taxids_under_annelida(limit)

        if limit and options['source'] != 'csv':
            taxids = taxids[:limit]

        taxids = [str(tid).strip() for tid in taxids if tid]
        unique_taxids = []
        seen = set()
        for tid in taxids:
            if tid not in seen:
                seen.add(tid)
                unique_taxids.append(tid)

        header = [
            'Phylum',
            'Class',
            'Order',
            'Family',
            'Geneus',
            'Species',
            'NCBI taxonomy id',
            'gbif_id',
            'worms_id',
            'genome',
            'mitogenome',
            'DNA',
            'DNA_link',
            'RNA',
            'RNA_link',
        ]

        self.stdout.write(f'Generating {output_path} for {len(unique_taxids)} taxonomy records...')
        with open(output_path, 'w', encoding='utf-8', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header, extrasaction='ignore')
            writer.writeheader()

            for idx, ncbi_id in enumerate(unique_taxids, start=1):
                ncbi_data, error = _fetch_ncbi_taxonomy(ncbi_id)
                if error:
                    self.stdout.write(self.style.WARNING(f'[{idx}] NCBI taxonomy fetch failed for {ncbi_id}: {error}'))
                classification = {}
                if ncbi_data and ncbi_data.get('classification'):
                    classification = {entry.get('rank', '').strip().lower(): entry.get('scientific_name', '').strip()
                                      for entry in ncbi_data['classification'] if entry.get('rank')}

                phylum = classification.get('phylum', '')
                class_name = classification.get('class', '')
                order = classification.get('order', '')
                family = classification.get('family', '')
                genus = classification.get('genus', '')
                species = ncbi_data.get('scientific_name', '') if ncbi_data else ''

                tax = Taxonomy.objects.filter(ncbi_taxonomy_id=ncbi_id).first()
                gbif_id = None
                worms_id = None
                if tax:
                    gbif_id = tax.gbif_id
                    worms_id = tax.worms_id

                if not gbif_id and species:
                    gbif_id = self._gbif_name_search(species)
                if not worms_id and species:
                    worms_id = self._worms_name_search(species)

                genome_label = ''
                if tax and tax.genomes.filter(assembly_level__icontains='chromosome').exists():
                    genome_label = 'chromosome'
                elif tax and tax.genomes.exists():
                    genome_label = 'scaffold'

                mitogenome_label = ''
                if tax:
                    mitogenomes = tax.mitogenomes.all()
                    if mitogenomes.filter(topology__icontains='Circular').exists():
                        mitogenome_label = 'Circular'
                    elif mitogenomes.exists():
                        mitogenome_label = 'Yes'

                ncbi_omics = _get_ncbid_omics_for_taxid(ncbi_id)
                dna_count = int(ncbi_omics.get('dna_count', 0) or 0)
                rna_count = int(ncbi_omics.get('rna_count', 0) or 0)
                dna_has = dna_count > 0 or (tax and bool(tax.dna_link))
                rna_has = rna_count > 0 or (tax and bool(tax.rna_link))

                dna_link = (tax.dna_link if tax else None) or _ncbi_sra_search_url(ncbi_id, _ncbi_sra_dna_query(ncbi_id))
                rna_link = (tax.rna_link if tax else None) or _ncbi_sra_search_url(ncbi_id, _ncbi_sra_rna_query(ncbi_id))

                writer.writerow({
                    'Phylum': phylum,
                    'Class': class_name,
                    'Order': order,
                    'Family': family,
                    'Geneus': genus,
                    'Species': species,
                    'NCBI taxonomy id': ncbi_id,
                    'gbif_id': gbif_id or '',
                    'worms_id': worms_id or '',
                    'genome': genome_label,
                    'mitogenome': mitogenome_label,
                    'DNA': '1' if dna_has else '',
                    'DNA_link': dna_link or '',
                    'RNA': '1' if rna_has else '',
                    'RNA_link': rna_link or '',
                })

                if idx < len(unique_taxids):
                    time.sleep(delay)

        self.stdout.write(self.style.SUCCESS(f'Generated taxonomy results at {output_path}'))

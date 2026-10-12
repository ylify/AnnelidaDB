import http.client
import os
import csv
import json
import time
import urllib.request
import urllib.parse
import urllib.error
import xml.etree.ElementTree as ET
from datetime import datetime
from django.core.management.base import BaseCommand


NCBI_EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
NCBI_DEFAULT_REQUEST_DELAY = 0.34
NCBI_DEFAULT_RETRY_COUNT = 5
NCBI_DEFAULT_BACKOFF = 5


class Command(BaseCommand):
    help = 'Fetch NCBI SRA runinfo for a taxonomy ID and cache it for Transcriptome page display.'

    def add_arguments(self, parser):
        parser.add_argument('--taxid', type=str, default='6340', help='NCBI taxonomy ID to query (default: 6340).')
        parser.add_argument('--query-key', type=str, default=None, help='NCBI WebEnv query_key for pre-built search results.')
        parser.add_argument('--webenv', type=str, default=None, help='NCBI WebEnv token for pre-built search results.')
        parser.add_argument('--output-dir', type=str, default='static/data/transcriptome', help='Directory under BASE_DIR to save transcriptome cache.')
        parser.add_argument('--output-name', type=str, default=None, help='Output CSV filename within the output directory.')
        parser.add_argument('--max-results', type=int, default=0, help='Maximum number of SRA records to fetch; 0 means fetch all available results.')
        parser.add_argument('--runinfo-chunk-size', type=int, default=200, help='Batch size for fetching runinfo records.')
        parser.add_argument('--biosample-search-chunk-size', type=int, default=50, help='Batch size for searching BioSample accessions.')
        parser.add_argument('--biosample-fetch-chunk-size', type=int, default=200, help='Batch size for fetching BioSample metadata.')
        parser.add_argument('--skip-biosample', action='store_true', help='Skip BioSample metadata lookup to speed up fetch.')
        parser.add_argument('--incremental', action='store_true', help='Only fetch new SRA records not already cached.')
        parser.add_argument('--rna-only', action='store_true', help='Limit search to RNA/transcriptomic SRA records.')
        parser.add_argument('--user-api-url', type=str, default=None, help='Optional user API endpoint that returns runinfo CSV, used instead of NCBI eutils.')
        parser.add_argument('--user-api-token', type=str, default=None, help='Optional bearer token for user API authentication.')
        parser.add_argument('--request-delay', type=float, default=NCBI_DEFAULT_REQUEST_DELAY, help='Seconds to wait between NCBI requests.')
        parser.add_argument('--retry-count', type=int, default=NCBI_DEFAULT_RETRY_COUNT, help='Number of retries for transient NCBI errors.')
        parser.add_argument('--api-key', type=str, default=None, help='Optional NCBI API key for higher rate limits. Can also be set via NCBI_API_KEY environment variable.')
        parser.add_argument('--skip-omics-cache', action='store_true', help='Do not refresh the NCBI omics cache after fetching transcriptome runinfo.')

    def handle(self, *args, **options):
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        output_dir = os.path.join(base_dir, options['output_dir'])
        os.makedirs(output_dir, exist_ok=True)
        taxid = options['taxid']
        max_results = options['max_results']
        runinfo_chunk_size = options['runinfo_chunk_size']
        biosample_search_chunk_size = options['biosample_search_chunk_size']
        biosample_fetch_chunk_size = options['biosample_fetch_chunk_size']
        skip_biosample = options['skip_biosample']
        incremental = options['incremental']
        rna_only = options['rna_only']
        request_delay = options['request_delay']
        retry_count = options['retry_count']
        user_api_url = options['user_api_url']
        user_api_token = options['user_api_token']
        api_key = options['api_key'] or os.environ.get('NCBI_API_KEY')
        query_key = options['query_key']
        webenv = options['webenv']
        output_name = options['output_name']
        self.request_delay = max(0.0, request_delay)
        self.retry_count = max(1, retry_count)

        self.stdout.write(f'Fetching NCBI SRA runinfo for taxonomy {taxid}...')

        if query_key and webenv:
            csv_path = os.path.join(output_dir, output_name or 'sra_webenv_runinfo.csv')
        else:
            csv_path = os.path.join(output_dir, f'sra_{taxid}_runinfo.csv')

        existing_rows = []
        existing_run_ids = set()
        if incremental and os.path.isfile(csv_path):
            existing_rows = self._load_existing_rows(csv_path)
            existing_run_ids = self._load_existing_run_ids(csv_path)

        if user_api_url:
            self.stdout.write(f'Fetching runinfo from user API: {user_api_url}')
            runinfo_rows = self._fetch_runinfo_from_user_api(user_api_url, user_api_token)
            ids = [row.get('Run', '').strip() for row in runinfo_rows if row.get('Run')]
            if not ids:
                self.stdout.write(self.style.WARNING('User API returned no run IDs; aborting.'))
                return
        elif query_key and webenv:
            self.stdout.write(f'Fetching runinfo from NCBI WebEnv: query_key={query_key}')
            runinfo_rows = self._fetch_runinfo_from_webenv(query_key, webenv, api_key, max_results=max_results)
            ids = [row.get('Run', '').strip() for row in runinfo_rows if row.get('Run')]
            if not ids:
                self.stdout.write(self.style.WARNING('NCBI WebEnv returned no run IDs; aborting.'))
                return
        else:
            ids = self._search_sra_ids(taxid, max_results, api_key, rna_only)
            self.stdout.write(f'Found {len(ids)} SRA IDs')
            if not ids:
                self.stdout.write(self.style.WARNING('No SRA IDs found; aborting.'))
                return
            if incremental and existing_run_ids:
                new_ids = [item for item in ids if item not in existing_run_ids]
                self.stdout.write(f'Incremental mode: skipping {len(ids) - len(new_ids)} existing records, {len(new_ids)} new records')
                ids = new_ids
            if not ids:
                self.stdout.write(self.style.SUCCESS('No new SRA records to fetch.'))
                return
            runinfo_rows = self._fetch_runinfo(ids, api_key, runinfo_chunk_size, skip_biosample, biosample_search_chunk_size, biosample_fetch_chunk_size)
            if not runinfo_rows:
                self.stdout.write(self.style.ERROR('No runinfo rows returned.'))
                return

        if user_api_url and incremental and existing_run_ids:
            runinfo_rows = existing_rows + [row for row in runinfo_rows if row.get('Run', '').strip() not in existing_run_ids]

        elif query_key and webenv and incremental and existing_rows:
            runinfo_rows = existing_rows + [row for row in runinfo_rows if row.get('Run', '').strip() not in existing_run_ids]

        elif not user_api_url and not (query_key and webenv) and incremental and existing_rows:
            runinfo_rows = existing_rows + runinfo_rows

        runinfo_rows = self._dedupe_rows_by_experiment(runinfo_rows)

        meta_path = os.path.join(output_dir, f'sra_{taxid}_runinfo_meta.json')

        all_fieldnames = self._collect_fieldnames(runinfo_rows)
        with open(csv_path, 'w', encoding='utf-8', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=all_fieldnames, extrasaction='ignore')
            writer.writeheader()
            for row in runinfo_rows:
                writer.writerow(row)

        meta = {
            'taxid': taxid,
            'fetched_at': datetime.utcnow().isoformat() + 'Z',
            'record_count': len(runinfo_rows),
        }
        with open(meta_path, 'w', encoding='utf-8') as metafile:
            json.dump(meta, metafile, indent=2)

        self.stdout.write(self.style.SUCCESS(f'Saved {len(runinfo_rows)} transcriptome records to {csv_path}'))
        self.stdout.write(self.style.SUCCESS(f'Meta data written to {meta_path}'))

        if not options['skip_omics_cache']:
            self.stdout.write('Refreshing NCBI omics cache...')
            try:
                from api.views import _build_ncbi_omics_cache
                _build_ncbi_omics_cache()
                self.stdout.write(self.style.SUCCESS('NCBI omics cache refreshed.'))
            except Exception as exc:
                self.stdout.write(self.style.WARNING(f'Failed to refresh NCBI omics cache: {exc}'))

    def _load_existing_rows(self, csv_path):
        rows = []
        with open(csv_path, 'r', encoding='utf-8', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                rows.append(row)
        return rows

    def _load_existing_run_ids(self, csv_path):
        ids = set()
        with open(csv_path, 'r', encoding='utf-8', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                experiment = (row.get('Experiment') or '').strip()
                run_id = experiment or (row.get('Run') or '').strip()
                if run_id:
                    ids.add(run_id)
        return ids

    def _collect_fieldnames(self, rows):
        fieldnames = []
        seen = set()
        for row in rows:
            for key in row.keys():
                if key not in seen:
                    seen.add(key)
                    fieldnames.append(key)
        return fieldnames

    def _search_sra_ids(self, taxid, max_results, api_key=None, rna_only=False):
        term = self._build_sra_search_term(taxid, rna_only)
        ids = []
        retmax = 100000
        if max_results and max_results < retmax:
            retmax = max_results

        start = 0
        while True:
            params = {
                'db': 'sra',
                'term': term,
                'retmax': str(retmax),
                'retstart': str(start),
                'retmode': 'json',
            }
            if api_key:
                params['api_key'] = api_key
            url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
            response_data = self._fetch_url(url, timeout=30)
            data = json.loads(response_data.decode('utf-8'))
            page_ids = data.get('esearchresult', {}).get('idlist', [])
            if not page_ids:
                break
            ids.extend(page_ids)
            if max_results and len(ids) >= max_results:
                ids = ids[:max_results]
                break
            count = int(data.get('esearchresult', {}).get('count', '0') or '0')
            start += len(page_ids)
            if start >= count:
                break
        return ids

    def _build_sra_search_term(self, taxid, rna_only=False):
        base_term = f'txid{taxid}[Organism:exp]'
        if not rna_only:
            return base_term
        return f"{base_term} AND (RNA-Seq[LibraryStrategy] OR RNA[LibraryStrategy] OR transcriptomic[LibrarySource] OR RNA-Seq[Strategy] OR RNA[Strategy])"

    def _fetch_runinfo_from_user_api(self, api_url, api_token=None):
        request = urllib.request.Request(api_url)
        if api_token:
            request.add_header('Authorization', f'Bearer {api_token}')
        response_data = self._fetch_url(request, timeout=60)
        text = response_data.decode('utf-8', errors='replace')
        rows = list(csv.DictReader(text.splitlines()))
        return rows

    def _fetch_runinfo(self, ids, api_key=None, runinfo_chunk_size=200, skip_biosample=False, biosample_search_chunk_size=50, biosample_fetch_chunk_size=200):
        if not ids:
            return []
        rows = []
        chunks = [ids[i:i + runinfo_chunk_size] for i in range(0, len(ids), runinfo_chunk_size)]
        for chunk in chunks:
            params = {
                'db': 'sra',
                'id': ','.join(chunk),
                'rettype': 'runinfo',
                'retmode': 'text',
            }
            if api_key:
                params['api_key'] = api_key
            url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
            text = self._fetch_url(url, timeout=60).decode('utf-8', errors='replace')
            chunk_rows = list(csv.DictReader(text.splitlines()))
            rows.extend(chunk_rows)

        rows = self._dedupe_rows_by_experiment(rows)

        if not skip_biosample:
            biosample_accessions = self._extract_biosample_accessions(rows)
            if biosample_accessions:
                accession_to_id = self._resolve_biosample_accessions(biosample_accessions, api_key, biosample_search_chunk_size)
                biosample_meta = self._fetch_biosample_metadata(accession_to_id.values(), api_key, biosample_fetch_chunk_size)
            else:
                biosample_meta = {}
        else:
            biosample_meta = {}

        for row in rows:
            biosample_accession = (row.get('BioSample') or row.get('BioSample_accession') or row.get('BioSampleAccession') or '').strip()
            metadata = biosample_meta.get(biosample_accession, {}) if biosample_accession else {}
            row['BioSample'] = biosample_accession
            row['BioSampleTissue'] = metadata.get('tissue', '') or metadata.get('tissue/organ', '') or metadata.get('source tissue', '') or metadata.get('organism part', '') or metadata.get('body site', '') or metadata.get('cell type', '')
            row['DevelopmentStage'] = metadata.get('development stage', '') or metadata.get('development_stage', '') or metadata.get('developmentstage', '') or metadata.get('dev_stage', '') or metadata.get('age', '')

        return rows

    def _extract_biosample_accessions(self, rows):
        accessions = set()
        for row in rows:
            for key in ('BioSample', 'BioSample_accession', 'BioSampleAccession'):
                value = (row.get(key) or '').strip()
                if value:
                    accessions.add(value)
        return accessions

    def _resolve_biosample_accessions(self, accessions, api_key=None, search_chunk_size=50):
        accession_to_id = {}
        if not accessions:
            return accession_to_id

        accessions = sorted(accessions)
        groups = [accessions[i:i + search_chunk_size] for i in range(0, len(accessions), search_chunk_size)]
        for group in groups:
            term = ' OR '.join([f'{item}[Accession]' for item in group])
            params = {
                'db': 'biosample',
                'term': term,
                'retmode': 'json',
                'retmax': str(len(group)),
            }
            if api_key:
                params['api_key'] = api_key
            url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
            try:
                response_data = self._fetch_url(url, timeout=30)
                data = json.loads(response_data.decode('utf-8'))
            except Exception:
                continue
            ids = data.get('esearchresult', {}).get('idlist', [])
            if not ids:
                continue

            resolved = self._fetch_biosample_accession_map(ids, api_key)
            accession_to_id.update({k: v for k, v in resolved.items() if k not in accession_to_id})

        return accession_to_id

    def _dedupe_rows_by_experiment(self, rows):
        deduped = []
        seen = set()
        for row in rows:
            experiment = (row.get('Experiment') or '').strip()
            key = experiment or (row.get('Run') or '').strip()
            if not key:
                continue
            if key in seen:
                continue
            seen.add(key)
            row['Run'] = key
            deduped.append(row)
        return deduped

    def _fetch_biosample_accession_map(self, ids, api_key=None):
        mapping = {}
        if not ids:
            return mapping

        params = {
            'db': 'biosample',
            'id': ','.join(ids),
            'retmode': 'xml',
        }
        if api_key:
            params['api_key'] = api_key
        url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
        xml_text = self._fetch_url(url, timeout=60).decode('utf-8', errors='replace')
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return mapping

        for biosample in root.findall('.//BioSample'):
            accession = biosample.get('accession')
            if not accession:
                ids_elems = biosample.findall('.//Id')
                for id_elem in ids_elems:
                    if id_elem.get('db') == 'BioSample' and id_elem.text:
                        accession = id_elem.text.strip()
                        break
            numeric_id = biosample.get('id') or None
            if accession and numeric_id:
                mapping[accession] = numeric_id
        return mapping

    def _fetch_biosample_metadata(self, biosample_ids, api_key=None, fetch_chunk_size=200):
        metadata = {}
        if not biosample_ids:
            return metadata

        chunks = [list(biosample_ids)[i:i + fetch_chunk_size] for i in range(0, len(biosample_ids), fetch_chunk_size)]
        for chunk in chunks:
            params = {
                'db': 'biosample',
                'id': ','.join(chunk),
                'retmode': 'xml',
            }
            if api_key:
                params['api_key'] = api_key
            url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
            xml_text = self._fetch_url(url, timeout=60).decode('utf-8', errors='replace')
            try:
                root = ET.fromstring(xml_text)
            except ET.ParseError:
                continue

            for biosample in root.findall('.//BioSample'):
                accession = biosample.get('accession')
                if not accession:
                    ids = biosample.findall('.//Id')
                    for id_elem in ids:
                        if id_elem.get('db') == 'BioSample' and id_elem.text:
                            accession = id_elem.text.strip()
                            break
                if not accession:
                    continue

                attrs = {}
                for attribute in biosample.findall('.//Attribute'):
                    name = (attribute.get('attribute_name') or '').strip().lower()
                    value = (attribute.text or '').strip()
                    if name and value:
                        attrs[name] = value
                        harmonized = (attribute.get('harmonized_name') or '').strip().lower()
                        if harmonized:
                            attrs[harmonized] = value
                metadata[accession] = attrs

        return metadata

    def _fetch_url(self, request_or_url, timeout=60):
        attempt = 0
        while True:
            attempt += 1
            try:
                with urllib.request.urlopen(request_or_url, timeout=timeout) as response:
                    try:
                        data = response.read()
                    except http.client.IncompleteRead as exc:
                        data = exc.partial or b''
                if self.request_delay > 0:
                    time.sleep(self.request_delay)
                request_desc = request_or_url.full_url if hasattr(request_or_url, 'full_url') else request_or_url
                return data
            except urllib.error.HTTPError as exc:
                request_desc = request_or_url.full_url if hasattr(request_or_url, 'full_url') else request_or_url
                if exc.code == 429 or 500 <= exc.code < 600:
                    if attempt >= self.retry_count:
                        raise
                    wait = NCBI_DEFAULT_BACKOFF * attempt
                    self.stdout.write(self.style.WARNING(f'NCBI HTTP {exc.code} on request to {request_desc}. Retry {attempt}/{self.retry_count} after {wait}s.'))
                    time.sleep(wait)
                    continue
                raise
            except urllib.error.URLError as exc:
                request_desc = request_or_url.full_url if hasattr(request_or_url, 'full_url') else request_or_url
                if attempt >= self.retry_count:
                    raise
                self.stdout.write(self.style.WARNING(f'NCBI URL error {exc} on request to {request_desc}. Retry {attempt}/{self.retry_count} after {NCBI_DEFAULT_BACKOFF}s.'))
                time.sleep(NCBI_DEFAULT_BACKOFF)
                continue
            except http.client.IncompleteRead as exc:
                data = exc.partial or b''
                if data:
                    request_desc = request_or_url.full_url if hasattr(request_or_url, 'full_url') else request_or_url
                    self.stdout.write(self.style.WARNING(f'NCBI incomplete chunked read on request to {request_desc}; using partial response.'))
                    if self.request_delay > 0:
                        time.sleep(self.request_delay)
                    return data
                if attempt >= self.retry_count:
                    raise
                wait = NCBI_DEFAULT_BACKOFF * attempt
                request_desc = request_or_url.full_url if hasattr(request_or_url, 'full_url') else request_or_url
                self.stdout.write(self.style.WARNING(f'NCBI incomplete chunked read on request to {request_desc}. Retry {attempt}/{self.retry_count} after {wait}s.'))
                time.sleep(wait)
                continue

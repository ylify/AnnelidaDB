import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime

from django.conf import settings
from django.core.mail import send_mail
from django.core.management.base import BaseCommand

NCBI_EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
NCBI_DEFAULT_REQUEST_DELAY = 0.34
NCBI_DEFAULT_RETRY_COUNT = 5
NCBI_DEFAULT_BACKOFF = 5


class Command(BaseCommand):
    help = 'Check NCBI Assembly for new Annelida genomes and email changes if any are detected.'

    def add_arguments(self, parser):
        parser.add_argument('--taxid', type=str, default='6340', help='NCBI taxonomy ID to query (default: 6340).')
        parser.add_argument('--output-dir', type=str, default='static/data/genome', help='Directory under BASE_DIR in which to store genome cache.')
        parser.add_argument('--cache-file', type=str, default=None, help='Cache filename to save current genome state.')
        parser.add_argument('--request-delay', type=float, default=NCBI_DEFAULT_REQUEST_DELAY, help='Seconds to wait between NCBI requests.')
        parser.add_argument('--retry-count', type=int, default=NCBI_DEFAULT_RETRY_COUNT, help='Number of retries for transient NCBI errors.')
        parser.add_argument('--api-key', type=str, default=None, help='Optional NCBI API key for higher rate limits. Can also be set via NCBI_API_KEY environment variable.')
        parser.add_argument('--recipient', type=str, default='liyunlong@hkbu.edu.hk', help='Email recipient for update notifications.')
        parser.add_argument('--subject-prefix', type=str, default=settings.EMAIL_SUBJECT_PREFIX if hasattr(settings, 'EMAIL_SUBJECT_PREFIX') else '', help='Subject prefix for notification emails.')
        parser.add_argument('--no-email', action='store_true', help='Do not send email even if changes are detected.')

    def handle(self, *args, **options):
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        output_dir = os.path.join(base_dir, options['output_dir'])
        os.makedirs(output_dir, exist_ok=True)

        taxid = options['taxid']
        cache_filename = options['cache_file'] or f'assembly_{taxid}_cache.json'
        cache_path = os.path.join(output_dir, cache_filename)
        self.request_delay = max(0.0, options['request_delay'])
        self.retry_count = max(1, options['retry_count'])
        api_key = options['api_key'] or os.environ.get('NCBI_API_KEY')
        recipient = options['recipient']
        subject_prefix = options['subject_prefix'] or ''
        send_email_flag = not options['no_email']

        self.stdout.write(f'Checking NCBI Assembly for taxid {taxid}...')
        previous_records = self._load_cache(cache_path)
        current_records = self._fetch_current_assemblies(taxid, api_key)

        if current_records is None:
            self.stdout.write(self.style.ERROR('Failed to retrieve current genome records.'))
            return

        if previous_records is None:
            self.stdout.write(self.style.WARNING('No existing cache found. Saving baseline cache and stopping without email.'))
            self._save_cache(cache_path, current_records)
            return

        changes = self._diff_records(previous_records, current_records)
        if not changes['new'] and not changes['removed'] and not changes['updated']:
            self.stdout.write(self.style.SUCCESS('No new or changed genome assemblies detected.'))
            self._save_cache(cache_path, current_records)
            return

        message = self._format_change_message(taxid, changes)
        subject = f"{subject_prefix}NCBI Genome update for taxid {taxid}"

        if send_email_flag:
            self._send_email(subject, message, recipient)
            self.stdout.write(self.style.SUCCESS(f'Notification email sent to {recipient}.'))
        else:
            self.stdout.write(self.style.WARNING('Changes detected but --send-email not provided. No email sent.'))

        self.stdout.write(message)
        self._save_cache(cache_path, current_records)

    def _search_assembly_ids(self, taxid, api_key=None):
        term = f'txid{taxid}[Organism:exp]'
        ids = []
        retmax = 500
        start = 0

        while True:
            params = {
                'db': 'assembly',
                'term': term,
                'retmode': 'json',
                'retmax': str(retmax),
                'retstart': str(start),
            }
            if api_key:
                params['api_key'] = api_key
            url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
            try:
                response_data = self._fetch_url(url, timeout=30)
            except Exception as exc:
                self.stdout.write(self.style.ERROR(f'Error fetching assembly ids: {exc}'))
                return None

            try:
                data = json.loads(response_data.decode('utf-8'))
            except ValueError:
                self.stdout.write(self.style.ERROR('Unable to parse NCBI esearch response.'))
                return None

            page_ids = data.get('esearchresult', {}).get('idlist', [])
            if not page_ids:
                break

            ids.extend(page_ids)
            count = int(data.get('esearchresult', {}).get('count', '0') or '0')
            start += len(page_ids)
            if start >= count:
                break

        return ids

    def _fetch_assembly_summaries(self, ids, api_key=None):
        if not ids:
            return []

        records = []
        chunk_size = 200
        for i in range(0, len(ids), chunk_size):
            chunk = ids[i:i + chunk_size]
            params = {
                'db': 'assembly',
                'id': ','.join(chunk),
                'retmode': 'json',
            }
            if api_key:
                params['api_key'] = api_key
            url = f"{NCBI_EUTILS_BASE}/esummary.fcgi?{urllib.parse.urlencode(params)}"
            try:
                response_data = self._fetch_url(url, timeout=60)
            except Exception as exc:
                self.stdout.write(self.style.ERROR(f'Error fetching assembly summaries: {exc}'))
                return None

            try:
                data = json.loads(response_data.decode('utf-8'))
            except ValueError:
                self.stdout.write(self.style.ERROR('Unable to parse NCBI esummary response.'))
                return None

            result = data.get('result', {})
            uids = result.get('uids', [])
            for uid in uids:
                item = result.get(uid, {})
                record = self._normalize_record(item)
                if record:
                    records.append(record)

        return records

    def _normalize_record(self, data):
        record = {'uid': str(data.get('uid', ''))}
        for key, value in data.items():
            if key == 'uid':
                continue
            record[key] = self._normalize_value(value)
        return record

    def _normalize_value(self, value):
        if value is None:
            return ''
        if isinstance(value, (str, int, float, bool)):
            return str(value)
        if isinstance(value, list):
            return '; '.join(self._normalize_value(item) for item in value if item is not None)
        if isinstance(value, dict):
            return json.dumps(value, sort_keys=True, ensure_ascii=False)
        return str(value)

    def _fetch_current_assemblies(self, taxid, api_key=None):
        ids = self._search_assembly_ids(taxid, api_key)
        if ids is None:
            return None

        self.stdout.write(f'Found {len(ids)} assembly records for taxid {taxid}.')
        return self._fetch_assembly_summaries(ids, api_key)

    def _load_cache(self, cache_path):
        if not os.path.isfile(cache_path):
            return None
        try:
            with open(cache_path, 'r', encoding='utf-8') as handle:
                return json.load(handle)
        except Exception:
            return None

    def _save_cache(self, cache_path, records):
        meta = {
            'checked_at': datetime.utcnow().isoformat() + 'Z',
            'record_count': len(records),
            'records': records,
        }
        with open(cache_path, 'w', encoding='utf-8') as handle:
            json.dump(meta, handle, indent=2, ensure_ascii=False)

    def _diff_records(self, previous_meta, current_records):
        previous_records = previous_meta.get('records', []) if isinstance(previous_meta, dict) else []
        previous_map = {self._record_key(record): record for record in previous_records}
        current_map = {self._record_key(record): record for record in current_records}

        new_keys = sorted(set(current_map) - set(previous_map))
        removed_keys = sorted(set(previous_map) - set(current_map))
        updated = []

        for key in sorted(set(previous_map) & set(current_map)):
            prev = previous_map[key]
            curr = current_map[key]
            diffs = self._compare_record_fields(prev, curr)
            if diffs:
                updated.append({'key': key, 'diffs': diffs, 'previous': prev, 'current': curr})

        return {'new': [current_map[key] for key in new_keys], 'removed': [previous_map[key] for key in removed_keys], 'updated': updated}

    def _record_key(self, record):
        return record.get('AssemblyAccession') or record.get('uid') or record.get('AssemblyName') or record.get('Title') or ''

    def _compare_record_fields(self, old, new):
        diffs = []
        keys = sorted(set(old.keys()) | set(new.keys()))
        for key in keys:
            old_value = old.get(key, '')
            new_value = new.get(key, '')
            if old_value != new_value:
                diffs.append({'field': key, 'old': old_value, 'new': new_value})
        return diffs

    def _format_change_message(self, taxid, changes):
        lines = [f'NCBI Assembly genome update detected for taxid {taxid}.', '']

        lines.append(f'New genome assemblies: {len(changes["new"])}')
        lines.append(f'Removed genome assemblies: {len(changes["removed"])}')
        lines.append(f'Updated genome assemblies: {len(changes["updated"])}')
        lines.append('')

        if changes['new']:
            lines.append('New assemblies:')
            for record in changes['new']:
                lines.append(self._format_record_summary(record))
            lines.append('')

        if changes['removed']:
            lines.append('Removed assemblies:')
            for record in changes['removed']:
                lines.append(self._format_record_summary(record))
            lines.append('')

        if changes['updated']:
            lines.append('Updated assemblies:')
            for item in changes['updated']:
                lines.append(f'- {item["key"]}')
                for diff in item['diffs']:
                    lines.append(f'    * {diff["field"]}: {diff["old"] or "<empty>"} -> {diff["new"] or "<empty>"}')
                lines.append('')

        lines.append('This report is generated automatically by AnnelidaDB.')
        return '\n'.join(lines)

    def _format_record_summary(self, record):
        accession = record.get('AssemblyAccession') or record.get('uid') or '<unknown accession>'
        title = record.get('Title') or record.get('AssemblyName') or record.get('Organism') or ''
        organism = record.get('Organism', '')
        status = record.get('AssemblyStatus', '')
        return f'- {accession} | {title} | {organism} | Status: {status}'

    def _send_email(self, subject, message, recipient):
        if not recipient:
            self.stdout.write(self.style.WARNING('No recipient provided for notification email.'))
            return

        send_mail(
            subject,
            message,
            settings.DEFAULT_FROM_EMAIL,
            [recipient],
            fail_silently=False,
        )

    def _fetch_url(self, request_or_url, timeout=60):
        attempt = 0
        while True:
            attempt += 1
            try:
                with urllib.request.urlopen(request_or_url, timeout=timeout) as response:
                    data = response.read()
                if self.request_delay > 0:
                    time.sleep(self.request_delay)
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

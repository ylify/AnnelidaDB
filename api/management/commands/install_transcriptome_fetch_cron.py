import os
import shlex
import subprocess
import urllib.parse
from django.conf import settings
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'Install or update a daily cron job to refresh transcriptome runinfo from an NCBI WebEnv/query_key search.'

    def add_arguments(self, parser):
        parser.add_argument('--hour', type=int, default=0, help='Hour of day for the daily refresh (0-23).')
        parser.add_argument('--minute', type=int, default=0, help='Minute of hour for the daily refresh (0-59).')
        parser.add_argument('--query-url', type=str, default=None, help='NCBI query URL containing query_key and WebEnv parameters.')
        parser.add_argument('--query-key', type=str, default=None, help='NCBI WebEnv query_key for pre-built search results.')
        parser.add_argument('--webenv', type=str, default=None, help='NCBI WebEnv token for pre-built search results.')
        parser.add_argument('--taxid', type=str, default='6340', help='NCBI taxonomy ID for the fetch command.')
        parser.add_argument('--output-dir', type=str, default='static/data/transcriptome', help='Directory under BASE_DIR to save transcriptome cache.')
        parser.add_argument('--output-name', type=str, default=None, help='Output CSV filename within the output directory.')
        parser.add_argument('--api-key', type=str, default=None, help='Optional NCBI API key to use in the cron job.')
        parser.add_argument('--python', type=str, default='python3', help='Python interpreter to use in the cron job.')
        parser.add_argument('--skip-omics-cache', action='store_true', help='Pass --skip-omics-cache to the fetch command.')
        parser.add_argument('--force', action='store_true', help='Replace existing cron entry if one already exists.')

    def handle(self, *args, **options):
        hour = options['hour']
        minute = options['minute']
        python_cmd = options['python']
        api_key = options['api_key']
        force = options['force']
        query_url = options['query_url']
        query_key = options['query_key']
        webenv = options['webenv']
        taxid = options['taxid']
        output_dir = options['output_dir']
        output_name = options['output_name']
        skip_omics_cache = options['skip_omics_cache']

        if hour < 0 or hour > 23:
            raise ValueError('Hour must be between 0 and 23.')
        if minute < 0 or minute > 59:
            raise ValueError('Minute must be between 0 and 59.')

        if query_url and (not query_key or not webenv):
            parsed = urllib.parse.urlparse(query_url)
            values = urllib.parse.parse_qs(parsed.query)
            query_key = query_key or (values.get('query_key') or [None])[0]
            webenv = webenv or (values.get('WebEnv') or [None])[0]

        if not query_key or not webenv:
            raise ValueError('query_key and WebEnv are required. Provide them directly or via --query-url.')

        project_dir = str(settings.BASE_DIR)
        scripts_dir = os.path.join(project_dir, 'scripts')
        os.makedirs(scripts_dir, exist_ok=True)
        script_path = os.path.join(scripts_dir, 'fetch_transcriptome_runinfo.sh')
        python_cmd_quoted = shlex.quote(str(python_cmd))
        manage_py = os.path.join(project_dir, 'manage.py')
        manage_py_quoted = shlex.quote(manage_py)
        query_key_quoted = shlex.quote(str(query_key))
        webenv_quoted = shlex.quote(str(webenv))
        output_dir_quoted = shlex.quote(str(output_dir))
        command_parts = [python_cmd_quoted, manage_py_quoted, 'fetch_sra_transcriptomes', f'--query-key={query_key_quoted}', f'--webenv={webenv_quoted}', f'--taxid={shlex.quote(str(taxid))}', f'--output-dir={output_dir_quoted}']
        if output_name:
            command_parts.append(f'--output-name={shlex.quote(str(output_name))}')
        if skip_omics_cache:
            command_parts.append('--skip-omics-cache')

        lines = [
            '#!/usr/bin/env bash',
            'set -e',
            f'cd {shlex.quote(project_dir)}',
        ]
        if api_key:
            lines.append(f'export NCBI_API_KEY={shlex.quote(str(api_key))}')
        lines.append(' '.join(command_parts))
        script_contents = '\n'.join(lines) + '\n'

        with open(script_path, 'w', encoding='utf-8') as handle:
            handle.write(script_contents)
        os.chmod(script_path, 0o755)

        cron_command = shlex.quote(script_path)
        cron_line = f"{minute} {hour} * * * {cron_command} >/dev/null 2>&1"
        marker = '# sra_transcriptome_fetch'
        cron_entry = f"{cron_line} {marker}"

        try:
            existing = subprocess.run(['crontab', '-l'], capture_output=True, text=True)
            current = existing.stdout.splitlines() if existing.returncode == 0 else []
        except FileNotFoundError:
            self.stdout.write(self.style.ERROR('crontab is not available on this system.'))
            return

        updated = []
        for line in current:
            if marker in line:
                continue
            updated.append(line)

        if current and not force and any(marker in line for line in current):
            self.stdout.write(self.style.SUCCESS('Existing transcriptome fetch cron entry already present. Use --force to replace it.'))
            return

        updated.append(cron_entry)
        new_crontab = '\n'.join(updated).strip() + '\n'
        proc = subprocess.run(['crontab', '-'], input=new_crontab, text=True, capture_output=True)
        if proc.returncode != 0:
            self.stdout.write(self.style.ERROR(f'Failed to install crontab: {proc.stderr.strip()}'))
            return

        self.stdout.write(self.style.SUCCESS('Installed daily transcriptome fetch cron job.'))
        self.stdout.write(self.style.SUCCESS(f'Scheduled command: {cron_entry}'))

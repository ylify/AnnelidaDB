import os
import shlex
import subprocess
from django.conf import settings
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = 'Install or update a daily cron job to refresh the NCBI omics cache.'

    def add_arguments(self, parser):
        parser.add_argument('--hour', type=int, default=2, help='Hour of day for the daily refresh (0-23).')
        parser.add_argument('--minute', type=int, default=0, help='Minute of hour for the daily refresh (0-59).')
        parser.add_argument('--api-key', type=str, default=None, help='Optional NCBI API key to use in the cron job.')
        parser.add_argument('--python', type=str, default='python3', help='Python interpreter to use in the cron job.')
        parser.add_argument('--force', action='store_true', help='Replace existing cron entry if one already exists.')

    def handle(self, *args, **options):
        hour = options['hour']
        minute = options['minute']
        python_cmd = options['python']
        api_key = options['api_key']
        force = options['force']

        if hour < 0 or hour > 23:
            raise ValueError('Hour must be between 0 and 23.')
        if minute < 0 or minute > 59:
            raise ValueError('Minute must be between 0 and 59.')

        project_dir = str(settings.BASE_DIR)
        scripts_dir = os.path.join(project_dir, 'scripts')
        os.makedirs(scripts_dir, exist_ok=True)
        script_path = os.path.join(scripts_dir, 'update_ncbi_omics_cache.sh')
        python_cmd_quoted = shlex.quote(str(python_cmd))
        manage_py = os.path.join(project_dir, 'manage.py')
        manage_py_quoted = shlex.quote(manage_py)

        lines = [
            '#!/usr/bin/env bash',
            'set -e',
            f'cd {shlex.quote(project_dir)}',
        ]
        if api_key:
            lines.append(f'export NCBI_API_KEY={shlex.quote(api_key)}')
        lines.append(f'{python_cmd_quoted} {manage_py_quoted} build_ncbi_omics_cache --force')
        script_contents = '\n'.join(lines) + '\n'

        with open(script_path, 'w', encoding='utf-8') as handle:
            handle.write(script_contents)
        os.chmod(script_path, 0o755)

        cron_command = shlex.quote(script_path)
        cron_line = f"{minute} {hour} * * * {cron_command} >/dev/null 2>&1"
        marker = '# ncbi_omics_cache_refresh'
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
            self.stdout.write(self.style.SUCCESS('Existing NCBI omics refresh cron entry already present. Use --force to replace it.'))
            return

        updated.append(cron_entry)
        new_crontab = '\n'.join(updated).strip() + '\n'
        proc = subprocess.run(['crontab', '-'], input=new_crontab, text=True, capture_output=True)
        if proc.returncode != 0:
            self.stdout.write(self.style.ERROR(f'Failed to install crontab: {proc.stderr.strip()}'))
            return

        self.stdout.write(self.style.SUCCESS('Installed daily NCBI omics cache refresh cron job.'))
        self.stdout.write(self.style.SUCCESS(f'Scheduled command: {cron_entry}'))

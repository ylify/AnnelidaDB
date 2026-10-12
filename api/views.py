from django.db.models import Count, Q
from django.http import JsonResponse, HttpResponse, HttpResponseNotFound, HttpResponseForbidden, FileResponse, StreamingHttpResponse
from django.shortcuts import render, get_object_or_404, redirect
from django.urls import reverse
from django.utils.html import format_html
from django.utils.safestring import mark_safe
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_GET, require_POST
from django.conf import settings
from django.core.cache import cache
from django.core.mail import send_mail
from django.core.validators import validate_email
from django.views.decorators.cache import cache_page
from django.core.exceptions import ValidationError
from django.utils.text import get_valid_filename
from .graph_email import is_graph_email_enabled, send_graph_email
import collections
from collections import Counter
from functools import lru_cache
from .models import Taxonomy, Genome, Mitogenome, Symbiont, SingleCell, ContigMapping, SequenceFileIndex, EggnogAnnotation, SubmissionRequest
import ast
import json
import html
import os
import glob
import fnmatch
import subprocess
import sys
import tempfile
import zipfile
import io
import csv
import re
import uuid
import hashlib
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
import shutil
import gzip
import mimetypes
import time
import traceback
import threading
from datetime import datetime, timedelta

from django.db import close_old_connections
from django.views.decorators.http import require_GET


def get_taxonomy_by_identifier(identifier):
    identifier = (identifier or '').strip()
    if not identifier:
        return None
    tax = Taxonomy.objects.filter(ncbi_taxonomy_id=identifier).first()
    if tax:
        return tax
    return Taxonomy.objects.filter(species__iexact=identifier).first()


def get_taxonomy_or_404(identifier):
    identifier = (identifier or '').strip()
    return get_object_or_404(
        Taxonomy,
        Q(ncbi_taxonomy_id=identifier) | Q(species__iexact=identifier)
    )


def get_taxonomy_for_species_name(species_name):
    if not species_name:
        return None
    species_name = species_name.strip()
    tax = Taxonomy.objects.filter(species__iexact=species_name).first()
    if tax:
        return tax
    return Taxonomy.objects.filter(species__icontains=species_name).first()

GBIF_CACHE_DIR = os.path.join(settings.BASE_DIR, 'gbif_cache')
GBIF_CACHE_TTL_SECONDS = 60 * 60 * 24 * 7
BLAST_CMD_TIMEOUT_SECONDS = 1800
MAX_SUBMISSION_UPLOAD_BYTES = int(1.5 * 1024 * 1024 * 1024)

# Global sequential job queue for BLAST and synteny submissions.
_JOB_QUEUE = collections.deque()
_JOB_QUEUE_LOCK = threading.Lock()
_JOB_QUEUE_CONDITION = threading.Condition(_JOB_QUEUE_LOCK)
_JOB_QUEUE_WORKER = None


def _start_job_queue_worker():
    global _JOB_QUEUE_WORKER
    if _JOB_QUEUE_WORKER is None or not _JOB_QUEUE_WORKER.is_alive():
        _JOB_QUEUE_WORKER = threading.Thread(target=_job_queue_worker, daemon=True)
        _JOB_QUEUE_WORKER.start()


def _job_queue_worker():
    while True:
        with _JOB_QUEUE_CONDITION:
            while not _JOB_QUEUE:
                _JOB_QUEUE_CONDITION.wait()
            func, args, kwargs = _JOB_QUEUE.popleft()

        try:
            func(*args, **kwargs)
        except Exception as exc:
            print(f'Job queue task failed: {exc}')


def _enqueue_job(func, *args, **kwargs):
    with _JOB_QUEUE_CONDITION:
        _JOB_QUEUE.append((func, args, kwargs))
        _start_job_queue_worker()
        _JOB_QUEUE_CONDITION.notify()

DOWNLOAD_JOBS = {}
DOWNLOAD_JOBS_LOCK = threading.Lock()
DOWNLOAD_JOB_DIR = os.path.join(settings.BASE_DIR, 'tmp', 'download_jobs')


def _download_job_state_path(job_id):
    return os.path.join(DOWNLOAD_JOB_DIR, f"{job_id}.json")


def _write_download_job_state(job_id, job):
    os.makedirs(DOWNLOAD_JOB_DIR, exist_ok=True)
    path = _download_job_state_path(job_id)
    tmp_path = path + '.tmp'
    with open(tmp_path, 'w', encoding='utf-8') as f:
        json.dump(job, f)
    os.replace(tmp_path, path)


def _load_download_job_state(job_id):
    path = _download_job_state_path(job_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return None


def _cleanup_old_download_jobs(max_age_seconds=3600):
    if not os.path.exists(DOWNLOAD_JOB_DIR):
        return
    now = time.time()
    for filename in os.listdir(DOWNLOAD_JOB_DIR):
        if not filename.endswith(('.json', '.zip')):
            continue
        file_path = os.path.join(DOWNLOAD_JOB_DIR, filename)
        try:
            if os.path.isfile(file_path) and now - os.path.getmtime(file_path) > max_age_seconds:
                os.remove(file_path)
        except Exception:
            pass


def _set_download_job(job_id, **kwargs):
    with DOWNLOAD_JOBS_LOCK:
        if job_id not in DOWNLOAD_JOBS:
            DOWNLOAD_JOBS[job_id] = {}
        DOWNLOAD_JOBS[job_id].update(kwargs)
        _write_download_job_state(job_id, DOWNLOAD_JOBS[job_id])


def _get_download_job(job_id):
    with DOWNLOAD_JOBS_LOCK:
        job = _load_download_job_state(job_id)
        if job is None:
            DOWNLOAD_JOBS.pop(job_id, None)
            return None
        DOWNLOAD_JOBS[job_id] = job
        return job


def _make_zip_archive(data_type, selected_items, output_path, progress_callback=None):
    if data_type == 'genome':
        selected_files = []
        for item in selected_items:
            tax = _resolve_taxonomy_identifier(item)
            selected_files.extend(_collect_genome_files(tax))
    elif data_type == 'mitogenome':
        selected_files = []
        for item in selected_items:
            tax = _resolve_taxonomy_identifier(item)
            selected_files.extend(_collect_mitogenome_files(tax))
    elif data_type == 'symbiont':
        selected_files = []
        for item in selected_items:
            selected_files.extend(_collect_symbiont_files(item))
    elif data_type == 'taxonomy':
        rows = Taxonomy.objects.filter(ncbi_taxonomy_id__in=selected_items)
        if not rows.exists():
            raise ValueError('No selected taxonomy entries found.')
        csv_buffer = io.StringIO()
        writer = csv.writer(csv_buffer)
        writer.writerow(['Species', 'NCBI Taxonomy ID', 'Class', 'Order', 'Family', 'Updated Clade 1', 'Updated Clade 2', 'Updated Clade 3', 'Updated Clade 4', 'Updated Clade 5', 'Updated Clade 6', 'Genome Available', 'Mito Available', 'DNA Available', 'RNA Available'])
        for row in rows:
            writer.writerow([
                row.species,
                row.ncbi_taxonomy_id,
                row.class_name or '',
                row.order or '',
                row.family or '',
                row.updated_clade1 or '',
                row.updated_clade2 or '',
                row.updated_clade3 or '',
                row.updated_clade4 or '',
                row.updated_clade5 or '',
                row.updated_clade6 or '',
                'Yes' if row.genomes.exists() else 'No',
                'Yes' if row.mitogenomes.exists() else 'No',
                'Yes' if getattr(row, 'dna_link', None) else 'No',
                'Yes' if getattr(row, 'rna_link', None) else 'No',
            ])
        selected_files = [('taxonomy_metadata.csv', csv_buffer.getvalue())]
    elif data_type == 'singlecell':
        rows = SingleCell.objects.select_related('taxonomy').filter(pk__in=[int(i) for i in selected_items if str(i).isdigit()])
        if not rows.exists():
            raise ValueError('No selected single-cell datasets found.')
        csv_buffer = io.StringIO()
        writer = csv.writer(csv_buffer)
        writer.writerow(['Species', 'Directory', 'Tissue', 'Stage', 'Clusters', 'Cells', 'Type', 'Year', 'Platform', 'BioProject', 'Journal', 'DOI', 'Title'])
        for row in rows:
            writer.writerow([
                row.species,
            row.dir or '',
                row.tissue or '',
                row.stage or '',
                row.clusters or '',
                row.cells or '',
                row.type or '',
                row.year or '',
                row.platform or '',
                row.bioproject or '',
                row.journal or '',
                row.doi or '',
                row.title or '',
            ])
        selected_files = [('singlecell_metadata.csv', csv_buffer.getvalue())]
        for row in rows:
            selected_files.extend(_collect_singlecell_rds_files(row))
    elif data_type == 'transcriptome':
        rows = _get_sample_transcriptome_rows()
        selected_indexes = [int(i) for i in selected_items if str(i).isdigit()]
        selected_rows = [rows[i] for i in selected_indexes if 0 <= i < len(rows)]
        if not selected_rows:
            raise ValueError('No selected transcriptome datasets found.')
        csv_buffer = io.StringIO()
        writer = csv.writer(csv_buffer)
        writer.writerow(['Species', 'Tissue', 'Library Type', 'Reads', 'Assembly Size', 'Complete Genes', 'Status'])
        for row in selected_rows:
            writer.writerow([
                row['species'],
                row['tissue'],
                row['library_type'],
                row['reads'],
                row['assembly_size'],
                row['complete_genes'],
                row['status'],
            ])
        selected_files = [('transcriptome_metadata.csv', csv_buffer.getvalue())]
    elif data_type == 'selected_files':
        selected_files = []
        allowed_base = os.path.normpath(os.path.join(settings.BASE_DIR, 'static', 'data'))
        for item in selected_items:
            if not item:
                continue
            normalized = os.path.normpath(os.path.join(settings.BASE_DIR, item.lstrip('/\\')))
            if not (normalized == allowed_base or normalized.startswith(allowed_base + os.sep)):
                raise ValueError('Invalid selected file path.')
            if not os.path.isfile(normalized):
                raise ValueError(f'File not found: {item}')
            archive_name = os.path.relpath(normalized, allowed_base)
            selected_files.append((normalized, archive_name))
        if not selected_files:
            raise ValueError('No selected files found.')
    else:
        raise ValueError('Unsupported download type.')

    if data_type in ('genome', 'mitogenome', 'symbiont', 'selected_files'):
        if not selected_files:
            raise ValueError('No downloadable files found for the selected items.')
        with zipfile.ZipFile(output_path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            total = len(selected_files)
            for idx, (abs_path, arcname) in enumerate(selected_files, start=1):
                archive.write(abs_path, arcname)
                if progress_callback:
                    progress_callback(idx, total, f'Compressing file {idx}/{total}...')
    elif data_type == 'singlecell':
        with zipfile.ZipFile(output_path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            total = len(selected_files)
            for idx, entry in enumerate(selected_files, start=1):
                if idx == 1:
                    archive.writestr(entry[0], entry[1])
                else:
                    archive.write(entry[0], entry[1])
                if progress_callback:
                    progress_callback(idx, total, f'Writing file {idx}/{total}...')
    else:
        with zipfile.ZipFile(output_path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            total = len(selected_files)
            for idx, (filename, content) in enumerate(selected_files, start=1):
                archive.writestr(filename, content)
                if progress_callback:
                    progress_callback(idx, total, f'Writing file {idx}/{total}...')


def _process_download_job(job_id, data_type, selected_items):
    try:
        _set_download_job(job_id, status='running', progress=5, message='Preparing download job...', filename='', file_path='')
        os.makedirs(DOWNLOAD_JOB_DIR, exist_ok=True)
        zip_filename = f"AnnelidaDB_{data_type}_download_{job_id}.zip"
        zip_path = os.path.join(DOWNLOAD_JOB_DIR, zip_filename)

        def progress_callback(index, total, message):
            percent = 10 + int(80 * index / total)
            _set_download_job(job_id, progress=percent, message=message)

        _set_download_job(job_id, progress=10, message='Collecting files...')
        _make_zip_archive(data_type, selected_items, zip_path, progress_callback=progress_callback)
        _set_download_job(job_id, status='completed', progress=100, message='Download package is ready.', filename=zip_filename, file_path=zip_path)
    except Exception as exc:
        _set_download_job(job_id, status='failed', progress=0, message=str(exc), error=str(exc))


def download_batch_job(request):
    if request.method != 'POST':
        return JsonResponse({'error': 'POST request required.'}, status=405)
    _cleanup_old_download_jobs()
    data_type = request.POST.get('data_type', 'genome')
    selected_items = request.POST.getlist('selected')
    if not selected_items:
        return JsonResponse({'error': 'No selected items provided.'}, status=400)
    job_id = uuid.uuid4().hex
    _set_download_job(job_id,
                      status='pending',
                      progress=0,
                      message='Queued for download...',
                      data_type=data_type,
                      selected_items=selected_items,
                      created=time.time())
    _enqueue_job(_process_download_job, job_id, data_type, selected_items)
    return JsonResponse({
        'job_id': job_id,
        'status_url': reverse('download_job_status', args=[job_id]),
    })


def download_job_status(request, job_id):
    _cleanup_old_download_jobs()
    job = _get_download_job(job_id)
    if not job:
        return JsonResponse({'error': 'Job not found.'}, status=404)
    response = {
        'status': job.get('status', 'pending'),
        'progress': job.get('progress', 0),
        'message': job.get('message', ''),
    }
    if job.get('status') == 'completed':
        response['download_url'] = reverse('download_job_file', args=[job_id, job.get('filename', '')])
    if job.get('status') == 'failed':
        response['error'] = job.get('error', 'Download failed.');
    return JsonResponse(response)


def download_job_file(request, job_id, filename):
    _cleanup_old_download_jobs()
    job = _get_download_job(job_id)
    if not job or job.get('status') != 'completed' or job.get('filename') != filename:
        return HttpResponseNotFound('Requested download file not found.')
    file_path = job.get('file_path')
    if not file_path or not os.path.isfile(file_path):
        return HttpResponseNotFound('Requested download file not found.')
    response = FileResponse(open(file_path, 'rb'), as_attachment=True, filename=filename)
    return response


def _ensure_gbif_cache_dir():
    if not os.path.exists(GBIF_CACHE_DIR):
        os.makedirs(GBIF_CACHE_DIR, exist_ok=True)


def _gbif_cache_file_path(kind, gbif_id):
    return os.path.join(GBIF_CACHE_DIR, f"{kind}_{gbif_id}.json")


def _load_gbif_cache(kind, gbif_id, max_age_seconds=GBIF_CACHE_TTL_SECONDS):
    path = _gbif_cache_file_path(kind, gbif_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            payload = json.load(fh)
    except Exception:
        return None

    timestamp = payload.get('timestamp')
    if not isinstance(timestamp, (int, float)):
        return None
    if time.time() - timestamp > max_age_seconds:
        return None
    return payload.get('data')


def _save_gbif_cache(kind, gbif_id, data):
    _ensure_gbif_cache_dir()
    path = _gbif_cache_file_path(kind, gbif_id)
    try:
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump({'timestamp': time.time(), 'data': data}, fh)
    except Exception:
        pass


def _parse_int(value):
    if not value:
        return 0
    try:
        return int(str(value).replace(',', '').strip())
    except ValueError:
        return 0


def _extract_busco_s_d(raw_value):
    if not raw_value:
        return None
    s_match = re.search(r'S:([\d\.]+%)', raw_value)
    d_match = re.search(r'D:([\d\.]+%)', raw_value)
    if s_match and d_match:
        return f"S:{s_match.group(1)},D:{d_match.group(1)}"
    return raw_value.strip() or None


def _normalize_species_name(name):
    if not name:
        return ''
    normalized = re.sub(r'[_\-]+', ' ', name)
    normalized = re.sub(r'[\.\,\;\:\'"\(\)\[\]\{\}]', ' ', normalized)
    normalized = re.sub(r'\s+', ' ', normalized)
    return normalized.strip().lower()


def _sanitize_species_filename(name):
    if not name:
        return ''
    sanitized = str(name).strip()
    sanitized = re.sub(r'\s+', '_', sanitized)
    sanitized = re.sub(r'[^A-Za-z0-9_\-.]', '_', sanitized)
    sanitized = re.sub(r'__+', '_', sanitized)
    sanitized = sanitized.strip('_')
    return sanitized or 'species'


def _send_notification_email(subject, message, recipients, html_message=None):
    if isinstance(recipients, str):
        recipients = [recipients]
    recipients = [r.strip() for r in recipients if isinstance(r, str) and r.strip()]
    valid_recipients = []
    for recipient in recipients:
        try:
            validate_email(recipient)
            valid_recipients.append(recipient)
        except ValidationError:
            continue

    if not valid_recipients:
        return

    def _send():
        try:
            if is_graph_email_enabled():
                content_type = 'HTML' if html_message else 'Text'
                body = html_message if html_message else message
                send_graph_email(subject, body, valid_recipients, content_type=content_type)
            else:
                send_mail(
                    subject,
                    message,
                    settings.DEFAULT_FROM_EMAIL,
                    valid_recipients,
                    fail_silently=False,
                    html_message=html_message,
                )
        except Exception as exc:
            print(f"Failed to send notification email to {valid_recipients}: {exc}")

    thread = threading.Thread(target=_send, daemon=True)
    thread.start()


def _safe_cache_key(prefix, key):
    return f"{prefix}:{hashlib.md5(key.encode('utf-8')).hexdigest()}"


def _generate_species_card_image_url(species_name):
    if not species_name:
        species_name = 'Species'
    normalized = str(species_name).strip()
    display_name = normalized
    if len(display_name) > 22:
        display_name = display_name[:22].rsplit(' ', 1)[0] + '...'
    hash_value = hashlib.md5(normalized.encode('utf-8')).hexdigest()
    color1 = f"#{hash_value[:6]}"
    color2 = f"#{hash_value[6:12]}"
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="400" height="260">
      <defs><linearGradient id="bg" x1="0" x2="1" y1="0" y2="1">
        <stop offset="0%" stop-color="{color1}"/>
        <stop offset="100%" stop-color="{color2}"/>
      </linearGradient></defs>
      <rect width="400" height="260" rx="24" ry="24" fill="url(#bg)" />
      <text x="50%" y="45%" dominant-baseline="middle" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="24" fill="#ffffff" opacity="0.95">Species</text>
      <text x="50%" y="65%" dominant-baseline="middle" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="18" fill="#ffffff" opacity="0.95">{html.escape(display_name)}</text>
    </svg>'''
    return 'data:image/svg+xml;charset=UTF-8,' + urllib.parse.quote(svg, safe='')


def _find_species_card_image_url(species_name):
    formatted_species = _format_species(species_name)
    if not formatted_species:
        return None

    species_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species)
    allowed_extensions = ['png', 'webp', 'gif', 'jpeg', 'jpg']

    if os.path.isdir(species_dir):
        for ext in allowed_extensions:
            image_name = f"{formatted_species}.{ext}"
            image_path = os.path.join(species_dir, image_name)
            if os.path.exists(image_path):
                return f"/static/data/{formatted_species}/{image_name}"
        for filename in os.listdir(species_dir):
            if '.' in filename:
                name, ext = filename.rsplit('.', 1)
                if name.lower() == formatted_species.lower() and ext.lower() in allowed_extensions:
                    return f"/static/data/{formatted_species}/{filename}"

    return None


def _find_species_card_image_url_or_placeholder(species_name):
    url = _find_species_card_image_url(species_name)
    return url if url else _generate_species_card_image_url(species_name)


BLAST_JOB_DIR = os.path.join(settings.BASE_DIR, 'tmp', 'blast_jobs')


def _blast_job_dir(job_id=None):
    if job_id:
        return os.path.join(BLAST_JOB_DIR, job_id)
    return BLAST_JOB_DIR


def _find_symbiont_file(base_dir, formatted_species, suffix):
    candidates = [
        f"{formatted_species}_SOB.{suffix}",
        f"{formatted_species}_SOB.complete.{suffix}",
    ]
    for filename in candidates:
        path = os.path.join(base_dir, filename)
        if os.path.exists(path):
            return filename

    glob_pattern = os.path.join(base_dir, f"{formatted_species}_SOB*.{suffix}")
    matches = sorted(glob.glob(glob_pattern))
    if matches:
        return os.path.basename(matches[0])

    suffix_matches = sorted(glob.glob(os.path.join(base_dir, f"*.{suffix}")))
    sob_matches = [path for path in suffix_matches if '_SOB' in os.path.basename(path)]
    matches = sob_matches or suffix_matches
    return os.path.basename(matches[0]) if matches else None


def _ensure_blast_job_dir(job_id=None):
    os.makedirs(_blast_job_dir(job_id), exist_ok=True)


def _blast_job_status_path(job_id):
    return os.path.join(_blast_job_dir(job_id), 'status.json')


def _blast_job_file_path(job_id, filename):
    return os.path.join(_blast_job_dir(job_id), filename)


def _write_blast_job_status(job_id, payload):
    _ensure_blast_job_dir(job_id)
    path = _blast_job_status_path(job_id)
    try:
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(payload, fh)
    except Exception:
        pass


def _read_blast_job_status(job_id):
    path = _blast_job_status_path(job_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path, 'r', encoding='utf-8') as fh:
            return json.load(fh)
    except Exception:
        return None


def _list_blast_job_files(job_id):
    job_dir = _blast_job_dir(job_id)
    if not os.path.isdir(job_dir):
        return []
    files = []
    for name in sorted(os.listdir(job_dir)):
        if name in {'status.json'}:
            continue
        path = os.path.join(job_dir, name)
        if os.path.isfile(path):
            files.append(name)
    return files


def _blast_job_file_url(request, job_id, filename):
    filename = os.path.basename(filename)
    return request.build_absolute_uri(reverse('blast_job_file', args=[job_id, filename]))


def _fetch_sequences_from_blast_db(db_prefix, subject_ids):
    if not subject_ids:
        return {}

    with tempfile.NamedTemporaryFile('w', delete=False, encoding='utf-8') as tmp_file:
        for subject_id in subject_ids:
            if subject_id:
                tmp_file.write(subject_id + '\n')
        entry_path = tmp_file.name

    try:
        blastdbcmd = ['blastdbcmd', '-db', db_prefix, '-entry_batch', entry_path, '-outfmt', '%f']
        result = subprocess.run(blastdbcmd, capture_output=True, text=True, timeout=BLAST_CMD_TIMEOUT_SECONDS)
        if result.returncode != 0 or not result.stdout:
            return {}

        sequences = {}
        current_id = None
        current_seq = []
        for line in result.stdout.splitlines():
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
        return sequences
    except Exception:
        return {}
    finally:
        try:
            os.remove(entry_path)
        except Exception:
            pass


def _backbone_tsv_path(filename):
    return os.path.join(settings.BASE_DIR, 'static', 'backbone', 'backbone', filename)


def _perform_blast(data):
    program = data.get('program', 'blastn')
    database = data.get('database', 'annelidadb_cds')
    evalue = str(data.get('evalue', '1e-5')).strip() or '1e-5'
    max_target_seqs = _parse_int(data.get('max_target_seqs', 100))
    sequence = str(data.get('sequence', '')).strip()

    if not sequence:
        raise ValueError('Please provide a query sequence or upload a FASTA file.')

    expected_type_map = {
        'blastn': 'nucleotide',
        'blastp': 'protein',
        'blastx': 'nucleotide',
        'tblastn': 'protein',
        'tblastx': 'nucleotide',
    }
    program_definitions = {
        'blastn': 'blastn',
        'blastp': 'blastp',
        'blastx': 'blastx',
        'tblastn': 'tblastn',
        'tblastx': 'tblastx',
    }
    database_bases = {
        'annelidadb_cds': {'prefix': 'cds_nucl/Annelida.CDS.blast_db', 'suffix': '.nal'},
        'annelidadb_genome': {'prefix': 'genome_nucl/Annelida.genome.blast_db', 'suffix': '.nal'},
        'annelidadb_protein': {'prefix': 'protein_prot/Annelida.protein.blast_db', 'suffix': '.pal'},
    }

    if program not in program_definitions:
        raise ValueError('Invalid BLAST program')
    if database not in database_bases:
        raise ValueError(f'Invalid BLAST database: {database}')
    if database not in {
        'blastn': ['annelidadb_cds', 'annelidadb_genome'],
        'blastp': ['annelidadb_protein'],
        'blastx': ['annelidadb_protein'],
        'tblastn': ['annelidadb_cds', 'annelidadb_genome'],
        'tblastx': ['annelidadb_cds', 'annelidadb_genome'],
    }[program]:
        raise ValueError('Invalid BLAST database for selected program')

    if not sequence.startswith('>'):
        sequence = '>UserQuery\n' + sequence

    detected_type = _detect_sequence_type(sequence)
    expected_type = expected_type_map.get(program)
    if detected_type and expected_type and detected_type != expected_type:
        raise ValueError(f'Selected BLAST program {program} expects a {expected_type} query, but the input appears to be {detected_type}.')

    db_info = database_bases[database]
    db_prefix = os.path.join(settings.BASE_DIR, 'static', 'blast_db', db_info['prefix'])
    db_file = db_prefix + db_info['suffix']
    if not os.path.exists(db_file):
        raise FileNotFoundError(f'BLAST database file not found: {os.path.basename(db_file)}')

    job_id = data.get('job_id')
    if job_id:
        temp_dir = _blast_job_dir(job_id)
    else:
        temp_dir = os.path.join(settings.BASE_DIR, 'tmp', 'blast')
    os.makedirs(temp_dir, exist_ok=True)

    timestamp = datetime.utcnow().strftime('%Y%m%dT%H%M%S%fZ')
    filename_base = f'blast_{uuid.uuid4().hex}_{timestamp}'
    query_file = os.path.join(temp_dir, filename_base + '.fasta')
    summary_file = os.path.join(temp_dir, filename_base + '.summary')
    raw_file = os.path.join(temp_dir, filename_base + '.raw')

    with open(query_file, 'w', encoding='utf-8') as f:
        f.write(sequence)

    summary_command = [
        program_definitions[program],
        '-query', query_file,
        '-db', db_prefix,
        '-out', summary_file,
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids',
        '-max_target_seqs', str(max_target_seqs),
        '-num_threads', '1',
        '-evalue', evalue,
    ]
    raw_command = [
        program_definitions[program],
        '-query', query_file,
        '-db', db_prefix,
        '-out', raw_file,
        '-outfmt', '0',
        '-max_target_seqs', str(max_target_seqs),
        '-num_threads', '1',
        '-evalue', evalue,
    ]

    try:
        subprocess.run(summary_command, capture_output=True, text=True, check=True, timeout=BLAST_CMD_TIMEOUT_SECONDS)

        if os.path.exists(summary_file):
            taxonomy_map = {}
            with open(summary_file, 'r', encoding='utf-8') as f:
                summary_lines = [line.rstrip('\n') for line in f if line.strip()]
            processed_lines = []
            for line in summary_lines:
                parts = line.split('\t')
                if len(parts) < 13:
                    continue
                staxids = parts[12].strip()
                taxonomy_id = ''
                species_name = ''
                if staxids:
                    taxonomy_id = staxids.split(';')[0].strip()
                    if taxonomy_id:
                        species_name = taxonomy_map.get(taxonomy_id)
                        if species_name is None:
                            tax = Taxonomy.objects.filter(ncbi_taxonomy_id=taxonomy_id).first()
                            species_name = tax.species if tax else ''
                            taxonomy_map[taxonomy_id] = species_name
                processed_lines.append('\t'.join(parts + [species_name]))
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write('\n'.join(processed_lines) + ('\n' if processed_lines else ''))

        subprocess.run(raw_command, capture_output=True, text=True, check=True, timeout=BLAST_CMD_TIMEOUT_SECONDS)

        def normalize_id(raw_id):
            if not raw_id:
                return ''
            raw_id = str(raw_id).strip()
            if raw_id.startswith('>'):
                raw_id = raw_id[1:]
            raw_id = raw_id.split()[0]
            raw_id = raw_id.replace('\\', '/')
            raw_id = os.path.basename(raw_id)
            return raw_id

        alignments = {}
        subject_index = {}
        if os.path.exists(raw_file):
            with open(raw_file, 'r', encoding='utf-8') as f:
                raw_lines = [line.rstrip('\n') for line in f]
            current_query = None
            current_subject = None
            current_lines = []
            for raw_line in raw_lines:
                stripped_line = raw_line.lstrip()
                if stripped_line.startswith('Query='):
                    if current_query and current_subject and current_lines:
                        alignments[(current_query, current_subject)] = '\n'.join(current_lines).strip()
                        subject_index.setdefault((current_query, normalize_id(current_subject)), []).append(current_subject)
                    current_query = normalize_id(stripped_line[len('Query='):])
                    current_subject = None
                    current_lines = [raw_line]
                elif stripped_line.startswith('>') and current_query:
                    if current_query and current_subject and current_lines:
                        alignments[(current_query, current_subject)] = '\n'.join(current_lines).strip()
                        subject_index.setdefault((current_query, normalize_id(current_subject)), []).append(current_subject)
                    current_subject = normalize_id(stripped_line[1:])
                    current_lines = [raw_line]
                elif current_query and current_subject:
                    current_lines.append(raw_line)
            if current_query and current_subject and current_lines:
                alignments[(current_query, current_subject)] = '\n'.join(current_lines).strip()
                subject_index.setdefault((current_query, normalize_id(current_subject)), []).append(current_subject)

        real_results = []
        if os.path.exists(summary_file):
            with open(summary_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.rstrip('\n')
                    if not line.strip():
                        continue
                    parts = line.split('\t')
                    if len(parts) < 13:
                        continue
                    query_id = parts[0].strip()
                    subject_id = normalize_id(parts[1].strip())
                    alignment_text = alignments.get((query_id, subject_id))
                    if alignment_text is None:
                        normalized_subject = normalize_id(subject_id)
                        candidate_subjects = subject_index.get((query_id, normalized_subject), [])
                        if len(candidate_subjects) == 1:
                            alignment_text = alignments.get((query_id, candidate_subjects[0]))
                        else:
                            alignment_text = 'Alignment text not available.'
                    raw_taxonomy_id = parts[12].strip()
                    taxonomy_id = None
                    if raw_taxonomy_id:
                        taxonomy_id = raw_taxonomy_id.split(';')[0].strip()
                    species_name = parts[13].strip() if len(parts) > 13 else ''
                    real_results.append({
                        'query': query_id,
                        'subject': subject_id,
                        'pident': parts[2].strip() + '%',
                        'length': parts[3].strip(),
                        'mismatch': parts[4].strip(),
                        'gapopen': parts[5].strip(),
                        'qstart': parts[6].strip(),
                        'qend': parts[7].strip(),
                        'sstart': parts[8].strip(),
                        'send': parts[9].strip(),
                        'evalue': parts[10].strip(),
                        'bitscore': parts[11].strip(),
                        'taxonomy_id': taxonomy_id,
                        'species': species_name or None,
                        'alignment': alignment_text,
                        'subject_sequence': None,
                    })

        if database in ['annelidadb_cds', 'annelidadb_protein']:
            db_info = database_bases.get(database)
            if db_info:
                db_prefix = os.path.join(settings.BASE_DIR, 'static', 'blast_db', db_info['prefix'])
                subject_ids = [normalize_id(result['subject']) for result in real_results if result.get('subject')]
                subject_seqs = _fetch_sequences_from_blast_db(db_prefix, list(dict.fromkeys(subject_ids)))
                for result in real_results:
                    normalized_subject = normalize_id(result['subject'])
                    if normalized_subject in subject_seqs:
                        result['subject_sequence'] = subject_seqs.get(normalized_subject)

        for result in real_results:
            subject_id = result['subject'].strip()
            contig_id = subject_id.split()[0]
            if '>' in subject_id:
                contig_id = subject_id.split('>')[1].split()[0]
            contig_id = contig_id.strip()

            if database in ['annelidadb_cds', 'annelidadb_genome']:
                if not result.get('taxonomy_id') or not result.get('species'):
                    try:
                        mapping = ContigMapping.objects.get(contig_id=contig_id)
                        result['species'] = result.get('species') or mapping.species
                        result['taxonomy_id'] = result.get('taxonomy_id') or mapping.taxonomy.ncbi_taxonomy_id
                    except ContigMapping.DoesNotExist:
                        found = False
                        for variation in [contig_id, contig_id.split('.')[0], contig_id.replace('_', '')]:
                            try:
                                mapping = ContigMapping.objects.get(contig_id=variation)
                                result['species'] = result.get('species') or mapping.species
                                result['taxonomy_id'] = result.get('taxonomy_id') or mapping.taxonomy.ncbi_taxonomy_id
                                found = True
                                break
                            except ContigMapping.DoesNotExist:
                                continue
                        if not found:
                            result['species'] = result.get('species') or f'Unknown ({contig_id})'
                            result['taxonomy_id'] = result.get('taxonomy_id') or None
            elif database == 'annelidadb_protein':
                if not result.get('taxonomy_id') or not result.get('species'):
                    try:
                        annotation = EggnogAnnotation.objects.get(query=contig_id)
                        result['species'] = result.get('species') or annotation.taxonomy.species
                        result['taxonomy_id'] = result.get('taxonomy_id') or annotation.taxonomy.ncbi_taxonomy_id
                    except EggnogAnnotation.DoesNotExist:
                        result['species'] = result.get('species') or f'Unknown ({contig_id})'
                        result['taxonomy_id'] = result.get('taxonomy_id') or None

            if not result.get('taxonomy_id') and result.get('species'):
                tax = Taxonomy.objects.filter(species=result['species']).first()
                if tax:
                    result['taxonomy_id'] = tax.ncbi_taxonomy_id

            if not result.get('taxonomy_id'):
                tax = _find_taxonomy_from_contig(contig_id)
                if tax:
                    result['taxonomy_id'] = tax.ncbi_taxonomy_id
                    result['species'] = result.get('species') or tax.species
                else:
                    tax = _find_taxonomy_by_gene_in_annotation_files(contig_id)
                    if tax:
                        result['taxonomy_id'] = getattr(tax, 'ncbi_taxonomy_id', None)
                        result['species'] = result.get('species') or getattr(tax, 'species', None)

            if result.get('taxonomy_id') and not result.get('species'):
                tax = Taxonomy.objects.filter(ncbi_taxonomy_id=result['taxonomy_id']).first()
                if tax:
                    result['species'] = tax.species

            if database in ['annelidadb_cds', 'annelidadb_protein']:
                if not result.get('subject_sequence'):
                    subject_seq = None
                    seq_type = 'cds' if database == 'annelidadb_cds' else 'protein'
                    candidate_ids = [contig_id]
                    if '.' in contig_id:
                        candidate_ids.append(contig_id.split('.')[0])
                    if '_' in contig_id:
                        candidate_ids.append(contig_id.replace('_', ''))
                    if '|' in contig_id:
                        parts = contig_id.split('|')
                        candidate_ids.extend([parts[0], parts[-1]])
                    candidate_ids = [cid for cid in dict.fromkeys([cid.strip() for cid in candidate_ids if cid.strip()])]

                    tax = None
                    species_name = result.get('species')
                    if result.get('taxonomy_id'):
                        tax = Taxonomy.objects.filter(ncbi_taxonomy_id=result['taxonomy_id']).first()
                        if tax:
                            species_name = tax.species
                            sequences = _find_gene_sequences(tax, contig_id)
                            seq_info = sequences.get(seq_type)
                            if seq_info:
                                subject_seq = seq_info.get('sequence')

                    if subject_seq is None:
                        genome_dir = None
                        if tax:
                            genome_dir = _species_genome_dir(tax.species)
                        elif species_name:
                            genome_dir = _species_genome_dir(species_name)
                        if genome_dir and os.path.isdir(genome_dir):
                            sequences = _find_gene_sequences(tax or species_name, contig_id, data_dir=genome_dir)
                            seq_info = sequences.get(seq_type)
                            if seq_info:
                                subject_seq = seq_info.get('sequence')

                    if subject_seq is None:
                        for candidate in candidate_ids:
                            entry = SequenceFileIndex.objects.filter(seq_id__iexact=candidate, seq_type=seq_type)
                            if result.get('taxonomy_id'):
                                entry = entry.filter(taxonomy__ncbi_taxonomy_id=result['taxonomy_id'])
                            entry = entry.order_by('id').first()
                            if entry:
                                fasta_path = _resolve_indexed_fasta_path(entry.file_path)
                                if fasta_path:
                                    header, seq = _find_sequence_in_fasta_file(fasta_path, candidate)
                                    if seq:
                                        subject_seq = seq
                                        break

                    if subject_seq is None:
                        for candidate in candidate_ids:
                            entry = SequenceFileIndex.objects.filter(seq_id__iexact=candidate, seq_type=seq_type).order_by('id').first()
                            if entry:
                                fasta_path = _resolve_indexed_fasta_path(entry.file_path)
                                if fasta_path:
                                    header, seq = _find_sequence_in_fasta_file(fasta_path, candidate)
                                    if seq:
                                        subject_seq = seq
                                        break

                    if subject_seq is None:
                        subject_seq = _find_sequence_in_species_data(seq_type, contig_id, tax=tax, species=species_name)

                    result['subject_sequence'] = subject_seq

            gene_families = _find_gene_family_memberships(subject_id)
            result['gene_family_display'] = ', '.join([f"{m['family_id']} ({m['source']})" for m in gene_families]) if gene_families else ''
            result['gene_family_memberships'] = gene_families

            if 'description' not in result:
                try:
                    annotation = EggnogAnnotation.objects.get(query=contig_id)
                    result['description'] = annotation.description or 'No description'
                    result['has_annotation'] = True
                except EggnogAnnotation.DoesNotExist:
                    result['description'] = 'No annotation available'
                    result['has_annotation'] = False

        return real_results
    finally:
        if not data.get('job_id'):
            for temp_path in (query_file, summary_file, raw_file):
                if temp_path and os.path.exists(temp_path):
                    try:
                        os.remove(temp_path)
                    except Exception:
                        pass


def _blast_job_worker(job_id, data):
    close_old_connections()
    email_address = str(data.get('email_address', '') or '').strip()
    job_url = str(data.get('job_url', '') or '').strip()
    _write_blast_job_status(job_id, {
        'status': 'running',
        'created_at': datetime.utcnow().isoformat() + 'Z',
        'program': data.get('program'),
        'database': data.get('database'),
        'email_address': email_address,
        'job_url': job_url,
    })
    job_status = 'running'
    completion_message = None
    try:
        data['job_id'] = job_id
        results = _perform_blast(data)
        job_status = 'success'
        completion_message = 'Your BLAST job has completed successfully.'
        status_payload = {
            'status': 'success',
            'results': results,
            'completed_at': datetime.utcnow().isoformat() + 'Z',
            'email_address': email_address,
            'job_url': job_url,
            'completion_email_sent': False,
        }
        if email_address:
            subject = 'Your AnnelidaDB BLAST job completed'
            body = [
                completion_message,
                f'Status: {job_status}',
            ]
            if job_url:
                body.append(f'View job status: {job_url}')
            _send_notification_email(subject, '\n\n'.join(body), email_address)
            status_payload['completion_email_sent'] = True
        _write_blast_job_status(job_id, status_payload)
    except subprocess.CalledProcessError as exc:
        job_status = 'error'
        completion_message = 'BLAST failed: ' + (exc.stderr.strip() if getattr(exc, 'stderr', None) else str(exc))
        status_payload = {
            'status': 'error',
            'message': completion_message,
            'completed_at': datetime.utcnow().isoformat() + 'Z',
            'email_address': email_address,
            'job_url': job_url,
            'completion_email_sent': False,
        }
        if email_address:
            subject = 'Your AnnelidaDB BLAST job failed'
            body = [
                completion_message,
                f'Status: {job_status}',
            ]
            if job_url:
                body.append(f'View job status: {job_url}')
            _send_notification_email(subject, '\n\n'.join(body), email_address)
            status_payload['completion_email_sent'] = True
        _write_blast_job_status(job_id, status_payload)
    except subprocess.TimeoutExpired:
        job_status = 'error'
        completion_message = f'BLAST timed out after {BLAST_CMD_TIMEOUT_SECONDS} seconds. Try a smaller query or review the database size.'
        status_payload = {
            'status': 'error',
            'message': completion_message,
            'completed_at': datetime.utcnow().isoformat() + 'Z',
            'email_address': email_address,
            'job_url': job_url,
            'completion_email_sent': False,
        }
        if email_address:
            subject = 'Your AnnelidaDB BLAST job failed'
            body = [
                completion_message,
                f'Status: {job_status}',
            ]
            if job_url:
                body.append(f'View job status: {job_url}')
            _send_notification_email(subject, '\n\n'.join(body), email_address)
            status_payload['completion_email_sent'] = True
        _write_blast_job_status(job_id, status_payload)
    except Exception as exc:
        job_status = 'error'
        completion_message = str(exc)
        status_payload = {
            'status': 'error',
            'message': completion_message,
            'traceback': traceback.format_exc(),
            'completed_at': datetime.utcnow().isoformat() + 'Z',
            'email_address': email_address,
            'job_url': job_url,
            'completion_email_sent': False,
        }
        if email_address:
            subject = 'Your AnnelidaDB BLAST job failed'
            body = [
                completion_message,
                f'Status: {job_status}',
            ]
            if job_url:
                body.append(f'View job status: {job_url}')
            _send_notification_email(subject, '\n\n'.join(body), email_address)
            status_payload['completion_email_sent'] = True
        _write_blast_job_status(job_id, status_payload)
    finally:
        close_old_connections()


def _run_blast_sync(data):
    return _perform_blast(data)


def _find_backbone_taxon_id_from_names(names):
    if not names:
        return None, None

    candidates = []
    for name in names:
        if not name:
            continue
        raw = str(name).strip()
        if raw:
            candidates.append(raw)
    if not candidates:
        return None, None

    lowercase_candidates = sorted({c.lower() for c in candidates})
    cache_key = _safe_cache_key('backbone_taxon_id', '|'.join(lowercase_candidates))
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('taxon_id') or None, cached.get('matched_name') or None

    normalized_candidates = {_normalize_species_name(c): c for c in candidates}
    raw_candidates = set(candidates)
    taxon_id = None
    matched_name = None
    taxon_path = _backbone_tsv_path('Taxon.tsv')
    try:
        csv.field_size_limit(10_000_000)
        with open(taxon_path, 'r', encoding='utf-8', errors='replace', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)
            for row in reader:
                if len(row) < 8:
                    continue
                row_sci = row[5].strip()
                row_canon = row[7].strip()
                if row_sci in raw_candidates or row_canon in raw_candidates:
                    taxon_id = row[0].strip()
                    matched_name = row_sci or row_canon
                    break
                if normalized_candidates and (_normalize_species_name(row_sci) in normalized_candidates or _normalize_species_name(row_canon) in normalized_candidates):
                    taxon_id = row[0].strip()
                    matched_name = row_sci or row_canon
                    break
    except FileNotFoundError:
        pass
    except Exception:
        pass

    cache.set(cache_key, {'taxon_id': taxon_id or '', 'matched_name': matched_name or ''}, 60 * 60 * 24)
    return taxon_id, matched_name


def _find_backbone_taxon_id(scientific_name):
    if not scientific_name:
        return None
    taxon_id, _ = _find_backbone_taxon_id_from_names([scientific_name])
    return taxon_id


def _fetch_backbone_distribution(tax, gbif_data=None):
    if not tax:
        return [], {}

    cache_key = _safe_cache_key('backbone_distribution', str(tax.ncbi_taxonomy_id))
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('records', []), cached.get('summary', {})

    names_to_check = [tax.species]
    if gbif_data:
        if gbif_data.get('canonical_name'):
            names_to_check.append(gbif_data.get('canonical_name'))
        if gbif_data.get('scientific_name'):
            names_to_check.append(gbif_data.get('scientific_name'))

    taxon_id, matched_name = _find_backbone_taxon_id_from_names(names_to_check)
    matched_source = None
    if matched_name:
        if gbif_data and matched_name == gbif_data.get('canonical_name'):
            matched_source = 'GBIF canonical name'
        elif gbif_data and matched_name == gbif_data.get('scientific_name'):
            matched_source = 'GBIF scientific name'
        else:
            matched_source = 'Direct species match'

    summary = {
        'taxonID': taxon_id,
        'matched_name': matched_name or '',
        'matched_source': matched_source,
        'record_count': 0,
        'top_countries': [],
        'top_locations': [],
    }
    records = []
    if not taxon_id:
        cache.set(cache_key, {'records': records, 'summary': summary}, 60 * 60 * 24)
        return records, summary

    distribution_path = _backbone_tsv_path('Distribution.tsv')
    country_counts = Counter()
    location_counts = Counter()
    total = 0
    try:
        csv.field_size_limit(10_000_000)
        with open(distribution_path, 'r', encoding='utf-8', errors='replace', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)
            for row in reader:
                if not row or row[0].strip() != taxon_id:
                    continue
                total += 1
                location_id = row[1].strip() if len(row) > 1 else ''
                locality = row[2].strip() if len(row) > 2 else ''
                country = row[3].strip() if len(row) > 3 else ''
                country_code = row[4].strip() if len(row) > 4 else ''
                location_remarks = row[5].strip() if len(row) > 5 else ''
                establishment_means = row[6].strip() if len(row) > 6 else ''
                life_stage = row[7].strip() if len(row) > 7 else ''
                occurrence_status = row[8].strip() if len(row) > 8 else ''
                threat_status = row[9].strip() if len(row) > 9 else ''
                source = row[10].strip() if len(row) > 10 else ''

                country_counts[country or 'Unknown'] += 1
                location_counts[location_id or locality or 'Unknown'] += 1

                if len(records) < 40:
                    records.append({
                        'locationID': location_id,
                        'locality': locality,
                        'country': country,
                        'countryCode': country_code,
                        'locationRemarks': location_remarks,
                        'establishmentMeans': establishment_means,
                        'lifeStage': life_stage,
                        'occurrenceStatus': occurrence_status,
                        'threatStatus': threat_status,
                        'source': source,
                    })
    except FileNotFoundError:
        pass
    except Exception:
        pass

    summary['record_count'] = total
    summary['top_countries'] = sorted(country_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    summary['top_locations'] = sorted(location_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    cache.set(cache_key, {'records': records, 'summary': summary}, 60 * 60 * 24)
    return records, summary


def _match_taxonomy_by_name(name):
    if not name:
        return None
    raw_name = name.strip()
    if not raw_name:
        return None

    tax = Taxonomy.objects.filter(species__iexact=raw_name).first()
    if tax:
        return tax

    alt_name = _normalize_species_name(raw_name)
    if alt_name and alt_name != raw_name.lower():
        tax = Taxonomy.objects.filter(species__iexact=alt_name).first()
        if tax:
            return tax

    tax = Taxonomy.objects.filter(species__icontains=raw_name).first()
    if tax:
        return tax

    if alt_name:
        tax = Taxonomy.objects.filter(species__icontains=alt_name).first()
        if tax:
            return tax

    return None


_PUBLICATIONS_BY_SPECIES = None
_GENOME_REPEAT_BY_SPECIES = None


def _open_csv_with_encoding(path):
    if not os.path.exists(path):
        return None
    encodings = ['utf-8-sig', 'utf-8', 'mac_roman', 'cp1252', 'latin-1']
    for enc in encodings:
        try:
            f = open(path, 'r', encoding=enc, newline='')
            f.readline()
            f.seek(0)
            return f
        except UnicodeDecodeError:
            try:
                f.close()
            except Exception:
                pass
            continue
    return open(path, 'r', encoding='latin-1', newline='')


def _clean_publication_authors(authors_raw):
    if not authors_raw:
        return ''
    authors = html.unescape(str(authors_raw).strip())
    authors = re.sub(r'\s*,\s*', ', ', authors)
    authors = re.sub(r'\s*\.\s*$', '', authors)
    authors = re.sub(r'\s+', ' ', authors)
    return authors


def _clean_publication_title(title_raw):
    if not title_raw:
        return ''
    title = html.unescape(str(title_raw).strip())
    title = re.sub(r'<\s*(?:i|em)\s*>', '__ITALIC_OPEN__', title, flags=re.I)
    title = re.sub(r'<\s*/\s*(?:i|em)\s*>', '__ITALIC_CLOSE__', title, flags=re.I)
    title = re.sub(r'<[^>]*>', '', title)
    title = title.replace('__ITALIC_OPEN__', '<span class="fst-italic">')
    title = title.replace('__ITALIC_CLOSE__', '</span>')
    title = re.sub(r'\s+', ' ', title)
    return title.strip()


def _load_publications_by_species():
    global _PUBLICATIONS_BY_SPECIES
    if _PUBLICATIONS_BY_SPECIES is not None:
        return _PUBLICATIONS_BY_SPECIES

    _PUBLICATIONS_BY_SPECIES = {}
    path = os.path.join(settings.BASE_DIR, 'publications.csv')
    file_handle = _open_csv_with_encoding(path)
    if not file_handle:
        return _PUBLICATIONS_BY_SPECIES

    with file_handle as f:
        first_line = f.readline()
        if first_line.strip().lower() == 'publications':
            # Skip legacy top-line title row in publications.csv
            pass
        else:
            f.seek(0)
        reader = csv.DictReader(f)
        for row in reader:
            species_key = _normalize_species_name(row.get('Species (____)') or row.get('Species') or '')
            if not species_key:
                continue
            publication = {
                'title': _clean_publication_title(row.get('Title') or ''),
                'journal': (row.get('Journal') or '').strip(),
                'doi': (row.get('DOI') or '').strip(),
                'publication_date': (row.get('Publication Date') or '').strip(),
                'authors': _clean_publication_authors(row.get('Authors') or ''),
                'source': (row.get('Source') or '').strip(),
                'pmid': (row.get('PMID') or '').strip(),
            }
            _PUBLICATIONS_BY_SPECIES.setdefault(species_key, []).append(publication)
    return _PUBLICATIONS_BY_SPECIES


def _find_publications_for_taxonomy(tax):
    if not tax or not tax.species:
        return []
    target = _normalize_species_name(tax.species)
    publications = _load_publications_by_species().get(target)
    if publications:
        return publications
    for key, value in _load_publications_by_species().items():
        if key == target or target in key or key in target:
            return value
    return []


def _load_genome_repeat_by_species():
    global _GENOME_REPEAT_BY_SPECIES
    if _GENOME_REPEAT_BY_SPECIES is not None:
        return _GENOME_REPEAT_BY_SPECIES

    _GENOME_REPEAT_BY_SPECIES = {}
    path = os.path.join(settings.BASE_DIR, 'genome.csv')
    file_handle = _open_csv_with_encoding(path)
    if not file_handle:
        return _GENOME_REPEAT_BY_SPECIES

    with file_handle as f:
        reader = csv.DictReader(f)
        for row in reader:
            species_key = _normalize_species_name(row.get('Species') or '')
            if not species_key:
                continue
            repeat_value = (row.get('Repeat (%)') or '').strip()
            if repeat_value:
                _GENOME_REPEAT_BY_SPECIES[species_key] = repeat_value
    return _GENOME_REPEAT_BY_SPECIES


def _find_repeat_for_taxonomy(tax):
    if not tax or not tax.species:
        return None
    target = _normalize_species_name(tax.species)
    repeat = _load_genome_repeat_by_species().get(target)
    if repeat:
        return repeat
    for key, value in _load_genome_repeat_by_species().items():
        if key == target or target in key or key in target:
            return value
    return None


def _is_valid_phylogeny_label(value):
    if isinstance(value, (int, float)):
        return False
    if not isinstance(value, str):
        return False
    text = value.strip()
    if not text:
        return False
    if text.lower() in {'null', 'none', 'nan', 'inf', '-inf'}:
        return False
    if re.fullmatch(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', text):
        return False
    return True


def _extract_phylogeny_label(node):
    if not isinstance(node, dict):
        return None
    for key in ('species_name', 'species', 'taxon', 'taxon_name', 'name', 'label', 'node_name', 'title'):
        value = node.get(key)
        if _is_valid_phylogeny_label(value):
            return value.strip()
    for key, value in node.items():
        if key in {'branch_length', 'length', 'distance', 'weight'}:
            continue
        if _is_valid_phylogeny_label(value):
            return value.strip()
    return None


def _build_phylogeny_node(node):
    if isinstance(node, dict):
        label = _extract_phylogeny_label(node)
        children = []
        for key in ('children', 'child', 'nodes', 'items', 'branches'):
            if key in node and isinstance(node[key], list):
                children = node[key]
                break

        built = {
            'label': label or ('Clade' if children else 'Unknown'),
        }

        if children:
            built['children'] = [_build_phylogeny_node(child) for child in children if child is not None]
            built['taxonomy_id'] = None
        else:
            tax = _match_taxonomy_by_name(label)
            built['taxonomy_id'] = tax.ncbi_taxonomy_id or tax.species if tax else None
        return built
    elif isinstance(node, list):
        return [_build_phylogeny_node(item) for item in node if item is not None]
    else:
        return None


def _detect_sequence_type(sequence):
    if not sequence:
        return None
    seq = re.sub(r'>.*', '', sequence)
    seq = re.sub(r'\s+', '', seq).upper()
    if not seq:
        return None
    seq = re.sub(r'[^A-Z]', '', seq)
    if not seq:
        return None

    atgc_count = sum(seq.count(base) for base in ['A', 'T', 'G', 'C'])
    if atgc_count / len(seq) > 0.5:
        return 'nucleotide'
    return 'protein'


def _read_genome_csv(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            species = row.get('Species', '').strip()
            if not species:
                continue
            taxonomy = None
            tax = get_taxonomy_for_species_name(species)
            if tax:
                taxonomy = {
                    'species': tax.species,
                    'class_name': tax.class_name,
                    'order': tax.order,
                    'family': tax.family,
                    'ncbi_taxonomy_id': tax.ncbi_taxonomy_id,
                }
            else:
                taxonomy = {
                    'species': species,
                    'class_name': None,
                    'order': None,
                    'family': None,
                    'ncbi_taxonomy_id': None,
                }
            completeness_genome = _extract_busco_s_d(row.get('Completeness (genome)', ''))
            completeness_gene = _extract_busco_s_d(row.get('Completeness (gene)', '') or row.get('Completeness (% gene, by ANNEVO)', '')) 
            rows.append({
                'taxonomy': taxonomy,
                'size_mb': row.get('Size (Mb)', '').strip() or None,
                'gc_content': row.get('GC (%)', '').strip() or None,
                'scaffold_n50': row.get('Scaffold N50', '').strip() or None,
                'assembly_level': row.get('Assembly level', '').strip() or None,
                'completeness_genome': completeness_genome,
                'gene_number': row.get('Gene number', '').strip() or None,
                'completeness_gene': _extract_busco_s_d(row.get('Completeness (gene)', '')),
                'gene_number_annevo': row.get('Gene number (by ANNEVO)', '').strip() or None,
                'completeness_gene_annevo': _extract_busco_s_d(row.get('Completeness (% gene, by ANNEVO)', ''))})
    return rows


def _load_eggnog_annotation(path, max_rows=None):
    header_cells = []
    rows = []
    truncated = False
    open_func = gzip.open if path.lower().endswith('.gz') else open
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                if not header_cells:
                    header_cells = [c.strip().lstrip('#') for c in line.split('\t')]
                continue
            if not line.strip():
                continue
            rows.append(line.rstrip('\n').split('\t'))
            if max_rows is not None and len(rows) >= max_rows:
                truncated = True
                break
    return header_cells, rows, truncated


def _stream_eggnog_annotation(path, max_rows=None):
    open_func = gzip.open if path.lower().endswith('.gz') else open

    def generate():
        header_cells = []
        rows_streamed = 0
        truncated = False
        first_row = True
        yield '{'
        yield '"header":'
        with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#'):
                    if not header_cells:
                        header_cells = [c.strip().lstrip('#') for c in line.split('\t')]
                    continue
                if not line.strip():
                    continue
                break
        if not header_cells:
            header_cells = []
        yield json.dumps(header_cells)
        yield ',"rows":['
        with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#'):
                    continue
                if not line.strip():
                    continue
                if not first_row:
                    yield ','
                row = line.rstrip('\n').split('\t')
                yield json.dumps(row)
                first_row = False
                rows_streamed += 1
                if max_rows is not None and rows_streamed >= max_rows:
                    truncated = True
                    break
        yield '],'
        yield '"file":' + json.dumps(os.path.basename(path)) + ','
        yield '"truncated":' + json.dumps(truncated)
        yield '}'

    return StreamingHttpResponse(generate(), content_type='application/json')


def _load_interpro_annotation(path, max_rows=None):
    interpro_col_names = [
        'Protein ID', 'MD5', 'Length', 'Analysis', 'Signature Accession',
        'Signature Description', 'Start', 'End', 'Score', 'Status', 'Date',
        'IPR Accession', 'IPR Description', 'GO', 'Pathway',
    ]
    rows = []
    truncated = False
    open_func = gzip.open if path.lower().endswith('.gz') else open
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip():
                continue
            rows.append(line.rstrip('\n').split('\t'))
            if max_rows is not None and len(rows) >= max_rows:
                truncated = True
                break
    headers = interpro_col_names[:max((len(r) for r in rows), default=0)]
    member_dbs = sorted(set(r[3] for r in rows if len(r) > 3 and r[3].strip()))
    return headers, rows, member_dbs, truncated


def _search_eggnog_annotation(path, query, max_rows=None):
    query_lower = query.lower()
    header_cells = []
    rows = []
    truncated = False
    open_func = gzip.open if path.lower().endswith('.gz') else open
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                if not header_cells:
                    header_cells = [c.strip().lstrip('#') for c in line.split('\t')]
                continue
            if not line.strip():
                continue
            if query_lower not in line.lower():
                continue
            rows.append(line.rstrip('\n').split('\t'))
            if max_rows is not None and len(rows) >= max_rows:
                truncated = True
                break
    return header_cells, rows, truncated


def _search_interpro_annotation(path, query, max_rows=None):
    query_lower = query.lower()
    interpro_col_names = [
        'Protein ID', 'MD5', 'Length', 'Analysis', 'Signature Accession',
        'Signature Description', 'Start', 'End', 'Score', 'Status', 'Date',
        'IPR Accession', 'IPR Description', 'GO', 'Pathway',
    ]
    rows = []
    member_dbs = set()
    truncated = False
    open_func = gzip.open if path.lower().endswith('.gz') else open
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip():
                continue
            if query_lower not in line.lower():
                continue
            row = line.rstrip('\n').split('\t')
            rows.append(row)
            if len(row) > 3 and row[3].strip():
                member_dbs.add(row[3].strip())
            if max_rows is not None and len(rows) >= max_rows:
                truncated = True
                break
    headers = interpro_col_names[:max((len(r) for r in rows), default=0)]
    return headers, rows, sorted(member_dbs), truncated


def _read_eggnog_header(path):
    open_func = gzip.open if path.lower().endswith('.gz') else open
    header_cells = []
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header_cells = [c.strip().lstrip('#') for c in line.split('\t')]
                break
    return header_cells


def _count_eggnog_rows(path, search_query=None):
    search_query = search_query.strip().lower() if search_query else None
    open_func = gzip.open if path.lower().endswith('.gz') else open
    total = 0
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            if search_query and search_query not in line.lower():
                continue
            total += 1
    return total


def _read_eggnog_rows(path, start, length, search_query=None):
    search_query = search_query.strip().lower() if search_query else None
    open_func = gzip.open if path.lower().endswith('.gz') else open
    rows = []
    skipped = 0
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            if search_query and search_query not in line.lower():
                continue
            if skipped < start:
                skipped += 1
                continue
            rows.append(line.rstrip('\n').split('\t'))
            if length is not None and len(rows) >= length:
                break
    return rows


def _read_interpro_header(path):
    return [
        'Protein ID', 'MD5', 'Length', 'Analysis', 'Signature Accession',
        'Signature Description', 'Start', 'End', 'Score', 'Status', 'Date',
        'IPR Accession', 'IPR Description', 'GO', 'Pathway',
    ]


def _count_interpro_rows(path, search_query=None, member_db=None):
    search_query = search_query.strip().lower() if search_query else None
    member_db = member_db.strip().lower() if member_db else None
    open_func = gzip.open if path.lower().endswith('.gz') else open
    total = 0
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if member_db:
                if len(row) <= 3 or row[3].strip().lower() != member_db:
                    continue
            if search_query and search_query not in line.lower():
                continue
            total += 1
    return total


def _read_interpro_rows(path, start, length, search_query=None, member_db=None):
    search_query = search_query.strip().lower() if search_query else None
    member_db = member_db.strip().lower() if member_db else None
    open_func = gzip.open if path.lower().endswith('.gz') else open
    rows = []
    skipped = 0
    with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if member_db:
                if len(row) <= 3 or row[3].strip().lower() != member_db:
                    continue
            if search_query and search_query not in line.lower():
                continue
            if skipped < start:
                skipped += 1
                continue
            rows.append(row)
            if length is not None and len(rows) >= length:
                break
    return rows


def _stream_interpro_annotation(path, max_rows=None):
    interpro_col_names = [
        'Protein ID', 'MD5', 'Length', 'Analysis', 'Signature Accession',
        'Signature Description', 'Start', 'End', 'Score', 'Status', 'Date',
        'IPR Accession', 'IPR Description', 'GO', 'Pathway',
    ]
    open_func = gzip.open if path.lower().endswith('.gz') else open

    def generate():
        rows_streamed = 0
        truncated = False
        member_dbs = set()
        first_row = True
        yield '{'
        headers = interpro_col_names
        yield '"header":' + json.dumps(headers)
        yield ',"rows":['
        with open_func(path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.strip():
                    continue
                row = line.rstrip('\n').split('\t')
                if len(row) > 3 and row[3].strip():
                    member_dbs.add(row[3].strip())
                if not first_row:
                    yield ','
                yield json.dumps(row)
                first_row = False
                rows_streamed += 1
                if max_rows is not None and rows_streamed >= max_rows:
                    truncated = True
                    break
        yield '],'
        yield '"member_dbs":' + json.dumps(sorted(member_dbs)) + ','
        yield '"file":' + json.dumps(os.path.basename(path)) + ','
        yield '"truncated":' + json.dumps(truncated)
        yield '}'

    return StreamingHttpResponse(generate(), content_type='application/json')


def _find_annotation_species_directory(taxonomy_id=None, species=None):
    if taxonomy_id:
        tax = Taxonomy.objects.filter(ncbi_taxonomy_id=taxonomy_id).first()
        if tax:
            return _format_species(tax.species), tax

    if species:
        formatted_species = _format_species(species)
        tax = get_taxonomy_for_species_name(species)
        return formatted_species, tax

    return None, None


def _find_taxonomy_from_contig(gene_id):
    if not gene_id:
        return None

    variations = [gene_id, gene_id.split('.')[0], gene_id.replace('_', ''), gene_id.replace('-', ''), gene_id.replace('.', '')]
    seen = set()
    for candidate in variations:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        try:
            mapping = ContigMapping.objects.get(contig_id=candidate)
            return mapping.taxonomy
        except ContigMapping.DoesNotExist:
            continue
    return None


def _find_taxonomy_by_gene_in_annotation_files(gene_id):
    if not gene_id:
        return None

    base_data_dir = os.path.join(settings.BASE_DIR, 'static', 'data')
    if not os.path.isdir(base_data_dir):
        return None

    normalized_target = gene_id.strip().lower()
    taxonomy_map = { _format_species(tax.species): tax for tax in Taxonomy.objects.all() }

    for species_dir in sorted(os.listdir(base_data_dir)):
        genome_dir = os.path.join(base_data_dir, species_dir, 'genome')
        if not os.path.isdir(genome_dir):
            continue

        for pattern in [
            '*.1.pep.fasta.emapper.annotations',
            '*.1.pep.fasta.emapper.annotations.gz',
            '*.emapper.annotations',
            '*.emapper.annotations.gz',
            '*.1.pep.fasta.tsv',
            '*.1.pep.fasta.tsv.gz',
            '*.pep.fasta.tsv',
            '*.pep.fasta.tsv.gz',
        ]:
            for file_path in sorted(glob.glob(os.path.join(genome_dir, pattern))):
                if not os.path.isfile(file_path):
                    continue
                open_func = gzip.open if file_path.lower().endswith('.gz') else open
                try:
                    with open_func(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
                        for line in f:
                            if not line.strip() or line.startswith('#'):
                                continue
                            row = line.rstrip('\n').split('\t')
                            if row and row[0].strip().lower() == normalized_target:
                                if species_dir in taxonomy_map:
                                    return taxonomy_map[species_dir]
                                matched = next((tax for formatted, tax in taxonomy_map.items() if formatted == species_dir), None)
                                if matched:
                                    return matched
                                species_name = species_dir.replace('_', ' ')
                                return type('TaxonomyFallback', (), {'species': species_name, 'ncbi_taxonomy_id': None})()
                except Exception:
                    continue
    return None


def _auto_select_sequence_files(data_dir, annotation_source=None):
    protein_files = [f for f in _find_sequence_file_paths(data_dir, 'protein') if 'longest' not in f['name'].lower()]
    cds_files = [f for f in _find_sequence_file_paths(data_dir, 'cds') if 'longest' not in f['name'].lower()]

    if annotation_source and 'ANNEVO' in annotation_source.upper():
        protein_with_annevo = [f for f in protein_files if 'ANNEVO' in f['name'].upper() and f['name'].lower().endswith('.gz')]
        cds_with_annevo = [f for f in cds_files if 'ANNEVO' in f['name'].upper() and f['name'].lower().endswith('.gz')]
        if protein_with_annevo:
            protein_files = protein_with_annevo
        else:
            protein_files = [f for f in protein_files if 'ANNEVO' in f['name'].upper()] or protein_files
        if cds_with_annevo:
            cds_files = cds_with_annevo
        else:
            cds_files = [f for f in cds_files if 'ANNEVO' in f['name'].upper()] or cds_files

    selected_protein = protein_files[0] if protein_files else None
    selected_cds = cds_files[0] if cds_files else None
    return selected_protein['path'] if selected_protein else None, selected_cds['path'] if selected_cds else None


def _find_eggnog_annotation_rows_for_gene_dir(formatted_species, gene_id):
    cache_key = f"eggnog_annotation_rows:{formatted_species}:{gene_id.strip().lower()}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    if not formatted_species:
        result = ([], [], None)
        cache.set(cache_key, result, 3600)
        return result
    data_dir = _species_genome_dir(formatted_species)
    if not data_dir or not os.path.isdir(data_dir):
        result = ([], [], None)
        cache.set(cache_key, result, 3600)
        return result

    annotation_patterns = [
        '*.1.pep.fasta.emapper.annotations',
        '*.1.pep.fasta.emapper.annotations.gz',
        '*.emapper.annotations',
        '*.emapper.annotations.gz',
    ]
    for pattern in annotation_patterns:
        for ann_path in sorted(glob.glob(os.path.join(data_dir, pattern))):
            if not os.path.isfile(ann_path):
                continue
            open_func = gzip.open if ann_path.lower().endswith('.gz') else open
            header_cells = []
            matching_rows = []
            with open_func(ann_path, 'rt', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#'):
                        if not header_cells:
                            header_cells = [c.strip().lstrip('#') for c in line.split('\t')]
                        continue
                    if not line.strip():
                        continue
                    row = line.rstrip('\n').split('\t')
                    if row and row[0].strip().lower() == gene_id.strip().lower():
                        matching_rows.append(row)
            if matching_rows:
                result = (header_cells, matching_rows, os.path.basename(ann_path))
                cache.set(cache_key, result, 3600)
                return result
    result = ([], [], None)
    cache.set(cache_key, result, 3600)
    return result


def _find_interproscan_rows_for_gene_dir(formatted_species, gene_id):
    if not formatted_species:
        return [], [], [], None
    data_dir = _species_genome_dir(formatted_species)
    if not data_dir or not os.path.isdir(data_dir):
        return [], [], [], None

    interpro_col_names = [
        'Protein ID', 'MD5', 'Length', 'Analysis', 'Signature Accession',
        'Signature Description', 'Start', 'End', 'Score', 'Status', 'Date',
        'IPR Accession', 'IPR Description', 'GO', 'Pathway',
    ]

    interpro_patterns = [
        '*.1.pep.fasta.tsv',
        '*.1.pep.fasta.tsv.gz',
        '*.pep.fasta.tsv',
        '*.pep.fasta.tsv.gz',
    ]
    for pattern in interpro_patterns:
        for tsv_path in sorted(glob.glob(os.path.join(data_dir, pattern))):
            if not os.path.isfile(tsv_path):
                continue
            matching_rows = []
            open_func = gzip.open if tsv_path.lower().endswith('.gz') else open
            with open_func(tsv_path, 'rt', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if not line.strip():
                        continue
                    row = line.rstrip('\n').split('\t')
                    if not row:
                        continue
                    if row[0].strip().lower() == gene_id.strip().lower():
                        matching_rows.append(row)
            if matching_rows:
                member_dbs = sorted(set(r[3] for r in matching_rows if len(r) > 3 and r[3].strip()))
                headers = interpro_col_names[:max((len(r) for r in matching_rows), default=0)]
                return headers, matching_rows, member_dbs, os.path.basename(tsv_path)
    return [], [], [], None


def _build_annotation_fields_from_file(header_cells, row):
    fields = []
    for index, label in enumerate(header_cells):
        value = row[index] if index < len(row) else ''
        if not value:
            continue
        formatted_value = _format_annotation_value(label, value)
        fields.append((label, formatted_value))
    return fields


def _find_interproscan_rows_for_gene(tax, gene_id):
    cache_key = f"gene_interpro:{tax.ncbi_taxonomy_id}:{gene_id}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    formatted_species = _format_species(tax.species)
    headers, rows, member_dbs, source = _find_interproscan_rows_for_gene_dir(formatted_species, gene_id)
    result = (headers, rows)
    cache.set(cache_key, result, 3600)
    return result


# ──────────────────────────────────────────────
# Home & Index
# ──────────────────────────────────────────────

def _find_data_directory(species, subdir):
    if not species:
        return None
    formatted_species = _format_species(species)
    candidate_dirs = [formatted_species, species]
    for candidate in candidate_dirs:
        directory = os.path.join(settings.BASE_DIR, 'static', 'data', candidate, subdir)
        if os.path.isdir(directory):
            return directory
    return None


GENOME_FASTA_PATTERNS = [
    '*.genomic.fasta.masked', '*.genomic.fasta', '*.fa', '*.fa.gz',
    '*.fasta', '*.fasta.gz', '*.genomic.fasta.gz', '*.genomic.fasta.masked.gz',
]


def _species_has_genome_data(species):
    genome_dir = _find_data_directory(species, 'genome')
    if not genome_dir:
        return False
    for pattern in GENOME_FASTA_PATTERNS:
        if glob.glob(os.path.join(genome_dir, pattern)):
            return True
    return False


def _species_has_chromosome_genome_data(species):
    genome_dir = _find_data_directory(species, 'genome')
    if not genome_dir:
        return False
    for root, _, filenames in os.walk(genome_dir):
        for filename in filenames:
            if '.chr.' in filename.lower():
                abs_path = os.path.join(root, filename)
                if os.path.isfile(abs_path):
                    return True
    return False


def _species_has_gene_model(species):
    genome_dir = _find_data_directory(species, 'genome')
    if not genome_dir:
        return False
    return bool(glob.glob(os.path.join(genome_dir, '*.gff3')) or glob.glob(os.path.join(genome_dir, '*.gff3.gz')))


def _species_has_symbiont_data(species):
    symbiont_files = _collect_symbiont_files(species)
    if not symbiont_files:
        return False
    for _, rel_path in symbiont_files:
        filename = os.path.basename(rel_path).lower()
        if filename.endswith(('.fna', '.fa', '.fasta', '.faa', '.ffn', '.gff3', '.gff3.gz')):
            return True
    return False


def _species_has_mitogenome_data(species):
    if not species:
        return False
    formatted_species = _format_species(species)
    folder_name = f"mitoz.{formatted_species}.mt.fasta.result"
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, folder_name)
    if not os.path.isdir(base_dir):
        return False
    for root, _, filenames in os.walk(base_dir):
        for filename in filenames:
            if filename:
                abs_path = os.path.join(root, filename)
                if os.path.isfile(abs_path):
                    return True
    return False


def _mitogenome_species_count():
    return Taxonomy.objects.filter(
        mitogenomes__isnull=False,
    ).exclude(
        species__isnull=True,
    ).exclude(
        species='',
    ).values('species').distinct().count()


def _species_has_transcriptome_data(species):
    if not species:
        return False
    rows, _ = _load_cached_transcriptome_rows()
    if rows is None:
        rows = _get_sample_transcriptome_rows()
    species_norm = _normalize_species_name(species)
    for row in rows:
        if not _is_rna_transcriptome_row(row):
            continue
        if _normalize_species_name(row.get('species')) == species_norm:
            return True
    return False


def _species_has_epigenome_data(species):
    if not species:
        return False
    tax = get_taxonomy_for_species_name(species)
    return _has_epigenome_csv_data(species) or bool(getattr(tax, 'epigenome_link', None))


def _find_proteome_files(species):
    if not species:
        return []
    formatted_species = _format_species(species)
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, 'proteome')
    files = []
    if os.path.isdir(base_dir):
        files = [path for path in sorted(glob.glob(os.path.join(base_dir, '*'))) if os.path.isfile(path)]
    if files:
        return files
    file_patterns = ['*.pep*', '*.faa*', '*proteome*', '*.protein*']
    for pattern in file_patterns:
        files.extend([path for path in sorted(glob.glob(os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, pattern))) if os.path.isfile(path)])
    return files


def _species_has_proteome_data(species):
    if not species:
        return False
    tax = get_taxonomy_for_species_name(species)
    return _has_proteome_csv_data(species) or bool(getattr(tax, 'proteome_link', None))


STATISTICS_CACHE_TTL_SECONDS = 60 * 60 * 24
STATISTICS_CACHE_DIR = os.path.join(settings.BASE_DIR, 'tmp', 'statistics_cache')
STATISTICS_CACHE_FILE = os.path.join(STATISTICS_CACHE_DIR, 'data_statistics_summary.json')
STATISTICS_CACHE_META = os.path.join(STATISTICS_CACHE_DIR, 'data_statistics_summary_meta.json')
STATISTICS_CACHE_KEY = 'data_statistics_summary'

_STATISTICS_REFRESH_THREAD = None
_STATISTICS_REFRESH_THREAD_LOCK = threading.Lock()

TRANSCRIPTOME_CACHE_DIR = os.path.join(settings.BASE_DIR, 'static', 'data', 'transcriptome')
TRANSCRIPTOME_ROWS_CACHE_KEY = 'transcriptome_rows_parsed'
TRANSCRIPTOME_PARSED_CACHE_FILE = os.path.join(TRANSCRIPTOME_CACHE_DIR, 'transcriptome_rows_cache.json')
TRANSCRIPTOME_PARSED_CACHE_META = os.path.join(TRANSCRIPTOME_CACHE_DIR, 'transcriptome_rows_cache_meta.json')


def _seconds_until_next_midnight():
    now = datetime.utcnow()
    tomorrow = (now + timedelta(days=1)).replace(hour=0, minute=0, second=0, microsecond=0)
    return int((tomorrow - now).total_seconds())


def _statistics_refresh_worker():
    while True:
        try:
            _get_statistics_data(force_refresh=True)
        except Exception:
            pass
        time.sleep(max(_seconds_until_next_midnight(), 60))


def start_statistics_refresh_thread():
    global _STATISTICS_REFRESH_THREAD
    with _STATISTICS_REFRESH_THREAD_LOCK:
        if _STATISTICS_REFRESH_THREAD is None or not _STATISTICS_REFRESH_THREAD.is_alive():
            _STATISTICS_REFRESH_THREAD = threading.Thread(target=_statistics_refresh_worker, daemon=True)
            _STATISTICS_REFRESH_THREAD.start()


_NCBI_OMICS_REFRESH_THREAD = None
_NCBI_OMICS_REFRESH_THREAD_LOCK = threading.Lock()


def _ncbi_omics_refresh_worker():
    while True:
        try:
            _get_ncbi_omics_cache(force_refresh=True)
        except Exception:
            pass
        time.sleep(max(_seconds_until_next_midnight(), 60))


def start_ncbi_omics_refresh_thread():
    global _NCBI_OMICS_REFRESH_THREAD
    with _NCBI_OMICS_REFRESH_THREAD_LOCK:
        if _NCBI_OMICS_REFRESH_THREAD is None or not _NCBI_OMICS_REFRESH_THREAD.is_alive():
            _NCBI_OMICS_REFRESH_THREAD = threading.Thread(target=_ncbi_omics_refresh_worker, daemon=True)
            _NCBI_OMICS_REFRESH_THREAD.start()


def _transcriptome_species_count():
    rows, _ = _load_cached_transcriptome_rows()
    if rows is None:
        rows = _get_sample_transcriptome_rows()
    return len({row.get('species') for row in rows if row.get('species') and _is_rna_transcriptome_row(row)})


def _transcriptome_cache_needs_refresh():
    if not os.path.isfile(TRANSCRIPTOME_PARSED_CACHE_META):
        return True
    try:
        with open(TRANSCRIPTOME_PARSED_CACHE_META, 'r', encoding='utf-8') as handle:
            meta = json.load(handle)
            fetched_at = meta.get('fetched_at')
            if not fetched_at:
                return True
            ts = datetime.fromisoformat(fetched_at.replace('Z', '+00:00'))
            return ts.date() < datetime.utcnow().date()
    except Exception:
        return True


def _load_cached_transcriptome_rows_from_disk():
    if not os.path.isfile(TRANSCRIPTOME_PARSED_CACHE_FILE) or _transcriptome_cache_needs_refresh():
        return None
    try:
        with open(TRANSCRIPTOME_PARSED_CACHE_FILE, 'r', encoding='utf-8') as handle:
            return json.load(handle)
    except Exception:
        return None


def _save_cached_transcriptome_rows(rows):
    os.makedirs(TRANSCRIPTOME_CACHE_DIR, exist_ok=True)
    try:
        with open(TRANSCRIPTOME_PARSED_CACHE_FILE, 'w', encoding='utf-8') as handle:
            json.dump(rows, handle, ensure_ascii=False)
        fetched_at = datetime.utcnow().isoformat() + 'Z'
        with open(TRANSCRIPTOME_PARSED_CACHE_META, 'w', encoding='utf-8') as meta_handle:
            json.dump({'fetched_at': fetched_at}, meta_handle, ensure_ascii=False, indent=2)
    except Exception:
        pass
    return rows


def _get_statistics_data(force_refresh=False):
    cached = None
    if not force_refresh:
        cached = cache.get(STATISTICS_CACHE_KEY)
        if cached is not None:
            cached['mitogenome_count'] = _format_count_if_large(_mitogenome_species_count())
            cached['symbiont_count'] = _format_count_if_large(Symbiont.objects.count())
            cached['chromosome_genome_species_count'] = _format_count_if_large(
                sum(1 for tax in Taxonomy.objects.all() if _species_has_chromosome_genome_data(tax.species))
            )
            cached['dna_species_count'] = _format_count_if_large(_get_dna_species_count())
            return cached

        cached = _load_cached_statistics_summary()
        if cached is not None:
            required_keys = [
                'taxonomy_count', 'genome_count', 'mitogenome_count', 'symbiont_count',
                'chromosome_genome_species_count', 'gene_model_species_count',
                'singlecell_species_count', 'singlecell_dataset_count', 'singlecell_cluster_count',
                'singlecell_number_total', 'genome_cds_sequence_count', 'genome_protein_sequence_count',
                'symbiont_cds_sequence_count', 'symbiont_protein_sequence_count',
                'rna_species_count', 'dna_species_count', 'epigenome_species_count',
                'proteome_species_count', 'oxford_plot_count',
            ]
            for key in required_keys:
                if key not in cached:
                    cached[key] = 0
            cached['mitogenome_count'] = _format_count_if_large(_mitogenome_species_count())
            cached['symbiont_count'] = _format_count_if_large(Symbiont.objects.count())
            cached['chromosome_genome_species_count'] = _format_count_if_large(
                sum(1 for tax in Taxonomy.objects.all() if _species_has_chromosome_genome_data(tax.species))
            )
            cached['dna_species_count'] = _format_count_if_large(_get_dna_species_count())
            cache.set(STATISTICS_CACHE_KEY, cached, _seconds_until_next_midnight())
            return cached

    chromosome_species_count = sum(1 for tax in Taxonomy.objects.all() if _species_has_chromosome_genome_data(tax.species))

    taxonomies = Taxonomy.objects.all()
    genome_count = sum(1 for tax in taxonomies if _species_has_genome_data(tax.species))
    gene_model_species_count = sum(1 for tax in taxonomies if _species_has_gene_model(tax.species))
    symbiont_count = Symbiont.objects.count()
    mitogenome_count = _mitogenome_species_count()
    epigenome_species_count = sum(1 for tax in taxonomies if _species_has_epigenome_data(tax.species))
    proteome_species_count = sum(1 for tax in taxonomies if _species_has_proteome_data(tax.species))

    genome_cds_sequence_count = 0
    genome_protein_sequence_count = 0
    symbiont_cds_sequence_count = 0
    symbiont_protein_sequence_count = 0
    for tax in taxonomies:
        genome_dir = _find_data_directory(tax.species, 'genome')
        if genome_dir:
            genome_cds_sequence_count += _count_sequence_records_in_directory(genome_dir, 'cds')
            genome_protein_sequence_count += _count_sequence_records_in_directory(genome_dir, 'protein')

        symbiont_files = _collect_symbiont_files(tax.species)
        if symbiont_files:
            symbiont_cds_sequence_count += sum(
                _count_fasta_records(abs_path)
                for abs_path, _ in symbiont_files
                if _is_sequence_file_name(os.path.basename(abs_path), 'cds')
            )
            symbiont_protein_sequence_count += sum(
                _count_fasta_records(abs_path)
                for abs_path, _ in symbiont_files
                if _is_sequence_file_name(os.path.basename(abs_path), 'protein')
            )

    transcriptome_species_count = _transcriptome_species_count()
    ncbi_cache = _get_ncbi_omics_cache()
    dna_species_count = _get_dna_species_count(ncbi_cache)
    singlecell_species_count = SingleCell.objects.values('species').distinct().count()
    singlecell_dataset_count = SingleCell.objects.count()
    singlecell_cluster_count = sum(_parse_integer_string(item.clusters) for item in SingleCell.objects.all())
    singlecell_number_total = sum(_parse_integer_string(item.cells) for item in SingleCell.objects.all())

    synteny_dir = os.path.join(settings.BASE_DIR, 'static', 'plots', 'Synteny')
    oxford_plot_count = 0
    if os.path.exists(synteny_dir):
        oxford_plot_count = len([f for f in glob.glob(os.path.join(synteny_dir, '*.AnLGs.oxford_plot.png'))])
        if oxford_plot_count == 0:
            oxford_plot_count = len([f for f in glob.glob(os.path.join(synteny_dir, '*.AnLGs.oxford_plot.pdf'))])

    stats = {
        'taxonomy_count': _format_count_if_large(Taxonomy.objects.count()),
        'genome_count': _format_count_if_large(genome_count),
        'mitogenome_count': _format_count_if_large(mitogenome_count),
        'symbiont_count': _format_count_if_large(symbiont_count),
        'chromosome_genome_species_count': _format_count_if_large(chromosome_species_count),
        'gene_model_species_count': _format_count_if_large(gene_model_species_count),
        'singlecell_species_count': _format_count_if_large(singlecell_species_count),
        'singlecell_dataset_count': _format_count_if_large(singlecell_dataset_count),
        'singlecell_cluster_count': _format_count_if_large(singlecell_cluster_count),
        'singlecell_number_total': _format_count_if_large(singlecell_number_total),
        'genome_cds_sequence_count': _format_count_if_large(genome_cds_sequence_count),
        'genome_protein_sequence_count': _format_count_if_large(genome_protein_sequence_count),
        'symbiont_cds_sequence_count': _format_count_if_large(symbiont_cds_sequence_count),
        'symbiont_protein_sequence_count': _format_count_if_large(symbiont_protein_sequence_count),
        'rna_species_count': _format_count_if_large(transcriptome_species_count),
        'dna_species_count': _format_count_if_large(dna_species_count),
        'epigenome_species_count': _format_count_if_large(epigenome_species_count),
        'proteome_species_count': _format_count_if_large(proteome_species_count),
        'oxford_plot_count': _format_count_if_large(oxford_plot_count),
    }
    cache_timeout = _seconds_until_next_midnight()
    cache.set(STATISTICS_CACHE_KEY, stats, cache_timeout)
    _save_cached_statistics_summary(stats)
    return stats


def _get_client_ip(request):
    forwarded = request.META.get('HTTP_X_FORWARDED_FOR', '')
    if forwarded:
        return forwarded.split(',')[0].strip()
    return request.META.get('REMOTE_ADDR', '').strip() or 'unknown'


def _geolocate_ip(ip):
    if not ip or ip in {'127.0.0.1', '::1', 'unknown'}:
        return None
    try:
        url = f'http://ip-api.com/json/{urllib.parse.quote(ip)}?fields=status,country,countryCode,lat,lon'
        with urllib.request.urlopen(url, timeout=3) as response:
            body = response.read().decode('utf-8', errors='ignore')
            data = json.loads(body)
            if data.get('status') == 'success':
                return {
                    'country': data.get('country', 'Unknown'),
                    'code': (data.get('countryCode') or '').upper(),
                    'lat': float(data.get('lat') or 0),
                    'lon': float(data.get('lon') or 0),
                }
    except (urllib.error.URLError, ValueError, json.JSONDecodeError, Exception):
        pass
    return None


def _save_visit_geo_data(data):
    path = os.path.join(settings.BASE_DIR, 'static', 'data', 'visit_geo.json')
    try:
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except Exception:
        pass


def _update_visit_geo_data(request, visit_geo):
    geo = _geolocate_ip(_get_client_ip(request))
    if not geo:
        return visit_geo
    for item in visit_geo:
        if item.get('code') == geo['code'] and item.get('country', '').lower() == geo['country'].lower():
            item['count'] = item.get('count', 0) + 1
            item['lat'] = geo['lat'] or item.get('lat', 0)
            item['lon'] = geo['lon'] or item.get('lon', 0)
            _save_visit_geo_data(visit_geo)
            return visit_geo
    visit_geo.append({
        'country': geo['country'],
        'code': geo['code'],
        'count': 1,
        'lat': geo['lat'],
        'lon': geo['lon'],
    })
    _save_visit_geo_data(visit_geo)
    return visit_geo


def _load_visit_geo_data():
    path = os.path.join(settings.BASE_DIR, 'static', 'data', 'visit_geo.json')
    if os.path.exists(path):
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                if isinstance(data, list):
                    return [
                        {
                            'country': item.get('country', 'Unknown'),
                            'code': item.get('code', '').upper(),
                            'count': int(item.get('count', 0) or 0),
                            'lat': float(item.get('lat', 0) or 0),
                            'lon': float(item.get('lon', 0) or 0),
                        }
                        for item in data
                        if isinstance(item, dict)
                    ]
        except Exception:
            pass
    return []


def home_page(request):
    visit_geo = _load_visit_geo_data()
    visit_geo = _update_visit_geo_data(request, visit_geo)
    total_visits = sum(item.get('count', 0) for item in visit_geo)
    top_countries = sorted(visit_geo, key=lambda item: item.get('count', 0), reverse=True)[:8]
    featured_species = []

    genome_species = Taxonomy.objects.filter(genomes__isnull=False).distinct().order_by('species')
    for tax in genome_species:
        image_url = _find_species_card_image_url(tax.species)
        if image_url:
            featured_species.append({
                'name': tax.species,
                'taxonomy_id': tax.ncbi_taxonomy_id,
                'image_url': image_url,
            })

    return render(request, 'index.html', {
        'stats': _get_statistics_data(),
        'visit_geo_data': mark_safe(json.dumps(visit_geo)),
        'visit_geo_total': total_visits,
        'visit_geo_top': top_countries,
        'featured_species': featured_species,
    })


def about_page(request):
    return render(request, 'about.html', {
        'page_title': 'About',
        'active_page': 'about',
    })


def user_guide(request):
    guide_path = os.path.join(
        settings.BASE_DIR, 'static', 'documents', 'AnnelidaDB_User_Guide.html'
    )
    return FileResponse(open(guide_path, 'rb'), content_type='text/html; charset=utf-8')


def user_guide_pdf(request):
    pdf_path = os.path.join(
        settings.BASE_DIR, 'static', 'documents', 'AnnelidaDB_User_Guide.pdf'
    )
    return FileResponse(
        open(pdf_path, 'rb'),
        as_attachment=True,
        filename='AnnelidaDB_User_Guide.pdf',
        content_type='application/pdf',
    )


def _get_about_resource_map():
    return {
        'image_source': {
            'label': 'Image Source Records',
            'type': 'csv',
            'path': os.path.join(settings.BASE_DIR, 'static', 'source_file', 'image_source.csv'),
            'download_name': 'image_source.csv',
        },
        'commands': {
            'label': 'Commands Script',
            'type': 'text',
            'path': os.path.join(settings.BASE_DIR, 'static', 'source_file', 'commands.txt'),
            'download_name': 'commands.txt',
        },
        'genomes': {
            'label': 'Genome Source Records',
            'type': 'csv',
            'path': os.path.join(settings.BASE_DIR, 'static', 'source_file', 'Genome_source.csv'),
            'download_name': 'Genome_source.csv',
        },
    }


def _get_management_commands_preview():
    commands_dir = os.path.join(settings.BASE_DIR, 'api', 'management', 'commands')
    rows = []
    if os.path.isdir(commands_dir):
        for filename in sorted(os.listdir(commands_dir)):
            if not filename.endswith('.py') or filename == '__init__.py':
                continue
            command_name = filename[:-3]
            description = '-'
            file_path = os.path.join(commands_dir, filename)
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    for line in f:
                        match = re.search(r"help\s*=\s*['\"](.+?)['\"]", line)
                        if match:
                            description = match.group(1).strip()
                            break
            except Exception:
                pass
            rows.append({'command': command_name, 'description': description})
    return rows


def _load_csv_preview(path, page=1, page_size=20):
    columns = []
    rows = []
    total_rows = 0
    if not os.path.exists(path):
        return columns, rows, total_rows
    if page < 1:
        page = 1
    if page_size < 1:
        page_size = 20
    start_index = (page - 1) * page_size
    end_index = start_index + page_size
    try:
        with open(path, 'r', encoding='utf-8', errors='replace', newline='') as f:
            reader = csv.reader(f)
            for index, row in enumerate(reader):
                if index == 0:
                    columns = row
                    continue
                if total_rows >= start_index and total_rows < end_index:
                    rows.append(row)
                total_rows += 1
    except Exception:
        pass
    return columns, rows, total_rows


def _load_text_preview(path, page=1, page_size=50):
    lines = []
    total_lines = 0
    if not os.path.exists(path):
        return lines, total_lines
    if page < 1:
        page = 1
    if page_size < 1:
        page_size = 50
    start_index = (page - 1) * page_size
    end_index = start_index + page_size
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            for index, line in enumerate(f):
                if index >= start_index and index < end_index:
                    lines.append(line.rstrip('\n'))
                total_lines += 1
    except Exception:
        pass
    return lines, total_lines


def about_file_preview(request, file_key):
    page = 1
    page_size = 20
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1
    try:
        page_size = int(request.GET.get('page_size', '20'))
    except ValueError:
        page_size = 20

    resources = _get_about_resource_map()
    info = resources.get(file_key)
    if info is None:
        return HttpResponseNotFound('Resource not found.')

    if info['type'] == 'csv':
        columns, rows, total_rows = _load_csv_preview(info['path'], page=page, page_size=page_size)
        total_pages = max(1, (total_rows + page_size - 1) // page_size)
        return JsonResponse({
            'key': file_key,
            'label': info['label'],
            'type': 'csv',
            'columns': columns,
            'rows': rows,
            'page': page,
            'page_size': page_size,
            'total_rows': total_rows,
            'total_pages': total_pages,
            'download_url': reverse('about_file_download', args=[file_key]),
        })

    if info['type'] == 'text':
        rows, total_rows = _load_text_preview(info['path'], page=page, page_size=page_size)
        total_pages = max(1, (total_rows + page_size - 1) // page_size)
        return JsonResponse({
            'key': file_key,
            'label': info['label'],
            'type': 'text',
            'rows': rows,
            'page': page,
            'page_size': page_size,
            'total_rows': total_rows,
            'total_pages': total_pages,
            'download_url': reverse('about_file_download', args=[file_key]),
        })

    return HttpResponseNotFound('Unsupported resource type.')


def about_file_download(request, file_key):
    resources = _get_about_resource_map()
    info = resources.get(file_key)
    if info is None or not os.path.exists(info['path']):
        return HttpResponseNotFound('Resource not found.')

    response = FileResponse(open(info['path'], 'rb'), content_type='text/csv')
    response['Content-Disposition'] = 'inline; filename="{}"'.format(info['download_name'])
    return response


def submit_data_page(request):
    success_message = None
    errors = []
    form_values = {
        'data_type': 'other',
        'title': '',
        'submitter_name': '',
        'submitter_email': '',
        'source_url': '',
        'details': '',
    }

    success_link = None
    if request.method == 'POST':
        form_values['data_type'] = request.POST.get('data_type', 'other').strip()
        form_values['title'] = request.POST.get('title', '').strip()
        form_values['submitter_name'] = request.POST.get('submitter_name', '').strip()
        form_values['submitter_email'] = request.POST.get('submitter_email', '').strip()
        form_values['source_url'] = request.POST.get('source_url', '').strip()
        form_values['details'] = request.POST.get('details', '').strip()

        upload_files = request.FILES.getlist('data_files')
        md5_file = request.FILES.get('md5_file')

        if not form_values['title']:
            errors.append('Please provide a submission title.')
        if not form_values['details']:
            errors.append('Please provide submission details.')

        if upload_files:
            total_upload_size = sum(f.size for f in upload_files)
            if total_upload_size > MAX_SUBMISSION_UPLOAD_BYTES:
                errors.append('Total upload size must be less than 1.5 GB.')

        if not errors:
            submission = SubmissionRequest.objects.create(
                data_type=form_values['data_type'],
                title=form_values['title'],
                submitter_name=form_values['submitter_name'],
                submitter_email=form_values['submitter_email'],
                source_url=form_values['source_url'],
                details=form_values['details'],
            )

            saved_files = []
            if upload_files or md5_file:
                try:
                    upload_root = os.path.join(settings.MEDIA_ROOT, 'submission_uploads')
                    os.makedirs(upload_root, exist_ok=True)
                    submission_dir = os.path.join(upload_root, submission.access_token)
                    os.makedirs(submission_dir, exist_ok=True)

                    for uploaded_file in upload_files:
                        safe_name = get_valid_filename(os.path.basename(uploaded_file.name)) or f'file_{uuid.uuid4().hex}'
                        target_path = os.path.join(submission_dir, safe_name)
                        base_name, ext = os.path.splitext(safe_name)
                        counter = 1
                        while os.path.exists(target_path):
                            safe_name = f'{base_name}_{counter}{ext}'
                            target_path = os.path.join(submission_dir, safe_name)
                            counter += 1
                        with open(target_path, 'wb+') as dest_file:
                            for chunk in uploaded_file.chunks():
                                dest_file.write(chunk)
                        saved_files.append(safe_name)

                    if md5_file:
                        safe_name = get_valid_filename(os.path.basename(md5_file.name)) or f'md5_{uuid.uuid4().hex}'
                        target_path = os.path.join(submission_dir, safe_name)
                        base_name, ext = os.path.splitext(safe_name)
                        counter = 1
                        while os.path.exists(target_path):
                            safe_name = f'{base_name}_{counter}{ext}'
                            target_path = os.path.join(submission_dir, safe_name)
                            counter += 1
                        with open(target_path, 'wb+') as dest_file:
                            for chunk in md5_file.chunks():
                                dest_file.write(chunk)
                        saved_files.append(safe_name)

                    if saved_files:
                        submission.uploaded_files = saved_files
                        submission.save(update_fields=['uploaded_files'])
                except Exception as save_error:
                    submission.delete()
                    errors.append('File upload failed. Please try again later.')

            if not errors:
                recipients = ['liyunlong@hkbu.edu.hk']
                if form_values['submitter_email']:
                    recipients.insert(0, form_values['submitter_email'])

                detail_link = request.build_absolute_uri(
                    reverse('submission_request_detail', args=[submission.access_token])
                )

            subject = 'Thank you for your AnnelidaDB data submission'
            plain_message = (
                f'Hello {form_values["submitter_name"] or "Contributor"},\n\n'
                'We have received your data submission request. Thank you for your contribution!\n\n'
                f'Submission title: {form_values["title"]}\n'
                f'Submission type: {form_values["data_type"]}\n'
                f'Source URL: {form_values["source_url"] or "None"}\n\n'
                f'You can check the submission status here: {detail_link}\n\n'
                'We will review and process your submission as soon as possible. After approval, it will be published in the database.\n\n'
                'If there are any questions, we will contact you.\n\n'
                'This email is also copied to liyunlong@hkbu.edu.hk.'
            )
            html_message = f'''
                <html>
                  <body style="font-family: Arial, sans-serif; color: #333; line-height: 1.5;">
                    <h2 style="color: #2c3e50;">Thank you for your AnnelidaDB data submission</h2>
                    <p>Hello <strong>{form_values["submitter_name"] or "Contributor"}</strong>,</p>
                    <p>We have received your data submission request. Thank you for your contribution!</p>
                    <table style="border-collapse: collapse; width: 100%; max-width: 600px;">
                      <tr>
                        <td style="padding: 8px; border: 1px solid #ddd; background: #f7f7f7; width: 140px;">Submission title</td>
                        <td style="padding: 8px; border: 1px solid #ddd;">{form_values["title"]}</td>
                      </tr>
                      <tr>
                        <td style="padding: 8px; border: 1px solid #ddd; background: #f7f7f7;">Submission type</td>
                        <td style="padding: 8px; border: 1px solid #ddd;">{form_values["data_type"]}</td>
                      </tr>
                      <tr>
                        <td style="padding: 8px; border: 1px solid #ddd; background: #f7f7f7;">Source URL</td>
                        <td style="padding: 8px; border: 1px solid #ddd;">{(form_values["source_url"] and f'<a href="{form_values["source_url"]}" target="_blank">{form_values["source_url"]}</a>') or 'None'}</td>
                      </tr>
                    </table>
                    <p>You can check the submission status here:</p>
                    <p><a href="{detail_link}" target="_blank">View submission details</a></p>
                    <p>We will review and process your submission as soon as possible. Once approved, it will be published in the database.</p>
                    <p>If there are any questions, we will contact you.</p>
                    <p>This email is also copied to <a href="mailto:liyunlong@hkbu.edu.hk">liyunlong@hkbu.edu.hk</a>.</p>
                    <p>Best regards,<br>The AnnelidaDB Team</p>
                  </body>
                </html>
            '''
            _send_notification_email(subject, plain_message, recipients, html_message=html_message)

            success_link = request.build_absolute_uri(
                reverse('submission_request_detail', args=[submission.access_token])
            )
            success_message = 'Your submission has been received. Thank you for your contribution. We will review and publish it after approval.'
            form_values = {
                'data_type': 'other',
                'title': '',
                'submitter_name': '',
                'submitter_email': '',
                'source_url': '',
                'details': '',
            }

    return render(request, 'submit.html', {
        'active_page': 'submit',
        'success_message': success_message,
        'success_link': success_link,
        'errors': errors,
        'form_values': form_values,
    })


def submission_request_detail(request, token):
    try:
        submission = SubmissionRequest.objects.get(access_token=token)
    except SubmissionRequest.DoesNotExist:
        return HttpResponseNotFound('Submission not found.')

    update_message = None
    if request.method == 'POST':
        if not request.user.is_staff:
            return HttpResponseForbidden('You are not authorized to update submission status.')
        new_status = request.POST.get('status')
        new_admin_note = request.POST.get('admin_note', '').strip()
        if new_status in dict(SubmissionRequest.STATUS_CHOICES).keys():
            submission.status = new_status
        submission.admin_note = new_admin_note
        submission.save()
        update_message = 'Submission status updated successfully.'

    return render(request, 'submission_detail.html', {
        'submission': submission,
        'update_message': update_message,
        'status_choices': SubmissionRequest.STATUS_CHOICES,
        'show_admin_form': request.user.is_staff,
    })


def _load_json_tolerant(path):
    with open(path, 'r', encoding='utf-8') as f:
        raw_text = f.read()

    try:
        return json.loads(raw_text)
    except json.JSONDecodeError:
        cleaned = re.sub(r'/\*.*?\*/', '', raw_text, flags=re.S)
        cleaned = re.sub(r'//.*', '', cleaned)
        cleaned = re.sub(r',\s*(?=[}\]])', '', cleaned)
        try:
            return json.loads(cleaned)
        except json.JSONDecodeError:
            try:
                return ast.literal_eval(cleaned)
            except Exception as exc:
                raise ValueError(f'Invalid phylogeny JSON file: {exc}')


def _unwrap_phylogeny_root(raw_tree):
    if not isinstance(raw_tree, dict):
        return raw_tree

    for root_key in ('root', 'tree', 'phylogeny'):
        if root_key in raw_tree and isinstance(raw_tree[root_key], (dict, list)):
            return raw_tree[root_key]

    if len(raw_tree) == 1:
        root_key = next(iter(raw_tree))
        if root_key.lower() in {'root', 'tree', 'phylogeny'} and isinstance(raw_tree[root_key], (dict, list)):
            return raw_tree[root_key]

    return raw_tree


# ──────────────────────────────────────────────
# Newick Parser for Phylogenetic Trees
# ──────────────────────────────────────────────

def _parse_newick(content):
    """Parse Newick format tree string into hierarchical structure."""
    content = content.strip()
    if content.endswith(';'):
        content = content[:-1]
    
    pos = [0]
    
    def parse_node():
        if pos[0] >= len(content):
            return None
        
        if content[pos[0]] == '(':
            # Internal node
            pos[0] += 1
            children = []
            while content[pos[0]] != ')':
                child = parse_node()
                if child:
                    children.append(child)
                if content[pos[0]] == ',':
                    pos[0] += 1
            pos[0] += 1
            
            name = ''
            distance = 0.0
            
            # Parse label after closing parenthesis
            while pos[0] < len(content) and content[pos[0]] not in '(),;':
                if content[pos[0]] == ':':
                    pos[0] += 1
                    dist_str = ''
                    while pos[0] < len(content) and content[pos[0]] not in '(),;:':
                        dist_str += content[pos[0]]
                        pos[0] += 1
                    try:
                        distance = float(dist_str) if dist_str else 0.0
                    except ValueError:
                        distance = 0.0
                else:
                    name += content[pos[0]]
                    pos[0] += 1
            
            return {
                'name': name.strip() or 'Clade',
                'children': children,
                'distance': distance
            }
        else:
            # Leaf node
            name = ''
            while pos[0] < len(content) and content[pos[0]] not in '(),;:':
                name += content[pos[0]]
                pos[0] += 1
            
            distance = 0.0
            if pos[0] < len(content) and content[pos[0]] == ':':
                pos[0] += 1
                dist_str = ''
                while pos[0] < len(content) and content[pos[0]] not in '(),;':
                    dist_str += content[pos[0]]
                    pos[0] += 1
                try:
                    distance = float(dist_str) if dist_str else 0.0
                except ValueError:
                    distance = 0.0
            
            return {
                'name': name.strip() or 'Unknown',
                'children': [],
                'distance': distance
            }
    
    return parse_node()


def _convert_newick_to_echarts(node):
    """Convert parsed Newick node to ECharts tree format."""
    if not node:
        return None
    
    name = node.get('name', 'Unknown').strip()
    distance = node.get('distance', 0.0)
    children = node.get('children', [])
    
    # Match with taxonomy if it's a leaf node
    taxonomy_id = None
    if not children:
        tax = _match_taxonomy_by_name(name)
        taxonomy_id = tax.ncbi_taxonomy_id if tax else None
    
    result = {
        'name': name,
        'value': distance if distance > 0 else None,
        'taxonomy_id': taxonomy_id
    }
    
    if children:
        result['children'] = [
            _convert_newick_to_echarts(child) 
            for child in children 
            if child is not None
        ]
    
    return result


def get_phylogeny_tree(request):
    # Try to load from Newick format first
    tree_nwk_path = os.path.join(settings.BASE_DIR, 'static', 'tree', 'tree.nwk')
    tree_json_path = os.path.join(settings.BASE_DIR, 'static', 'tree', 'tree.json')
    
    built_tree = None
    
    if os.path.exists(tree_nwk_path):
        try:
            with open(tree_nwk_path, 'r', encoding='utf-8') as f:
                newick_content = f.read()
            parsed_tree = _parse_newick(newick_content)
            built_tree = _convert_newick_to_echarts(parsed_tree)
        except Exception as exc:
            return JsonResponse({'status': 'error', 'message': f'Failed to parse Newick tree: {exc}'})
    elif os.path.exists(tree_json_path):
        # Fallback to JSON if Newick not available
        try:
            raw_tree = _load_json_tolerant(tree_json_path)
            raw_tree = _unwrap_phylogeny_root(raw_tree)
            built_tree = _build_phylogeny_node(raw_tree)
        except Exception as exc:
            return JsonResponse({'status': 'error', 'message': f'Failed to load phylogeny tree: {exc}'})
    else:
        return JsonResponse({'status': 'error', 'message': 'Phylogeny tree file not found.'})
    
    return JsonResponse({'status': 'success', 'tree': built_tree})


def get_species_list(request):
    taxonomies = Taxonomy.objects.all().prefetch_related('genomes', 'mitogenomes')[:50]
    data = []
    for tax in taxonomies:
        genome_info = None
        first_genome = tax.genomes.first()
        if first_genome:
            genome_info = {
                'size_mb': first_genome.size_mb,
                'gc_content': first_genome.gc_content,
            }
        mito_info = None
        first_mito = tax.mitogenomes.first()
        if first_mito:
            mito_info = {'length_bp': first_mito.length_bp}
        data.append({
            'id': tax.id,
            'ncbi_id': tax.ncbi_taxonomy_id,
            'scientific_name': tax.species,
            'class_name': tax.class_name if tax.class_name else 'Unassigned',
            'family': tax.family if tax.family else 'Unknown Family',
            'genome': genome_info,
            'mitogenome': mito_info,
            'image_url': 'https://via.placeholder.com/300x200/203a43/ffffff?text=Annelida',
        })
    return JsonResponse({'status': 'success', 'data': data})


def query_species_by_name(request):
    name = request.GET.get('name', '').strip()
    if not name:
        return JsonResponse({'status': 'error', 'message': 'Missing species name.'})

    tax = _match_taxonomy_by_name(name)
    if not tax:
        return JsonResponse({'status': 'error', 'message': 'Species not found.'})

    return JsonResponse({
        'status': 'success',
        'species': {
            'name': tax.species,
            'taxonomy_id': tax.ncbi_taxonomy_id,
        }
    })


def search_species(request):
    query = request.GET.get('q', '').strip()
    if not query:
        return JsonResponse({'status': 'error', 'message': 'Missing search query.'})

    normalized_query = query.lower()
    search_filter = (
        Q(species__icontains=query) |
        Q(class_name__icontains=query) |
        Q(order__icontains=query) |
        Q(family__icontains=query)
    )

    taxonomies = Taxonomy.objects.filter(search_filter).distinct()[:50]

    def _match_rank(tax):
        fields = [tax.species, tax.class_name, tax.order, tax.family]
        exact_match = any(str(value).strip().lower() == normalized_query for value in fields if value)
        if exact_match:
            return 0
        starts_match = any(str(value).strip().lower().startswith(normalized_query) for value in fields if value)
        if starts_match:
            return 1
        contains_match = any(normalized_query in str(value).strip().lower() for value in fields if value)
        if contains_match:
            return 2
        return 3

    sorted_taxonomies = sorted(taxonomies, key=lambda tax: (_match_rank(tax), tax.species or ''))[:3]

    results = []
    for tax in sorted_taxonomies:
        results.append({
            'name': tax.species,
            'taxonomy_id': tax.ncbi_taxonomy_id,
            'class_name': tax.class_name,
            'order': tax.order,
            'family': tax.family,
        })

    return JsonResponse({
        'status': 'success',
        'query': query,
        'results': results,
        'count': len(results)
    })


# ──────────────────────────────────────────────
# BLAST
# ──────────────────────────────────────────────

def blast_page(request):
    return render(request, 'blast.html')


def tools_page(request):
    species_list = Taxonomy.objects.filter(genomes__isnull=False).distinct().order_by('species')
    species_options = [
        {
            'species': tax.species,
            'taxonomy_id': tax.ncbi_taxonomy_id,
            'species_dir': _format_species(tax.species),
        }
        for tax in species_list
    ]
    return render(request, 'tools.html', {
        'active_page': 'tools',
        'species_list': species_options,
    })


def synteny_upload_page(request):
    upload_info = None
    upload_error = None
    show_upload_info = False
    species_raw = ''
    token = request.GET.get('token')
    if token and request.method != 'POST':
        upload_info = _get_synteny_upload_results(token, request)
        if upload_info is None:
            upload_error = 'The requested synteny upload token was not found.'
        else:
            show_upload_info = True

    if request.method == 'POST':
        show_upload_info = True
        species_raw = request.POST.get('species', '').strip()
        user_email = request.POST.get('email', '').strip()
        pep_file = request.FILES.get('species_pep')
        bed_file = request.FILES.get('species_bed')
        if not species_raw:
            upload_error = 'Please enter a species name.'
        if user_email and not upload_error:
            try:
                validate_email(user_email)
            except ValidationError:
                upload_error = 'Please enter a valid email address.'
        if not upload_error and (not pep_file or not bed_file):
            upload_error = 'Please upload both .pep and .bed files.'
            upload_error = 'Please upload both .pep and .bed files.'
        else:
            uploaded_pep_name = os.path.basename(pep_file.name)
            uploaded_bed_name = os.path.basename(bed_file.name)
            if not uploaded_pep_name.lower().endswith('.pep'):
                upload_error = 'The first file must be a .pep file.'
            elif not uploaded_bed_name.lower().endswith('.bed'):
                upload_error = 'The second file must be a .bed file.'
            else:
                upload_root = _get_synteny_upload_root()
                upload_token = uuid.uuid4().hex
                request_dir = os.path.join(upload_root, upload_token)
                os.makedirs(request_dir, exist_ok=True)

                species_slug = _sanitize_species_filename(species_raw)
                renamed_pep_name = f"{species_slug}.pep"
                renamed_bed_name = f"{species_slug}.bed"
                uploaded_pep_storage_name = f"uploaded_pep_{uuid.uuid4().hex}.pep"
                uploaded_bed_storage_name = f"uploaded_bed_{uuid.uuid4().hex}.bed"

                uploaded_pep_path = os.path.join(request_dir, uploaded_pep_storage_name)
                uploaded_bed_path = os.path.join(request_dir, uploaded_bed_storage_name)
                renamed_pep_path = os.path.join(request_dir, renamed_pep_name)
                renamed_bed_path = os.path.join(request_dir, renamed_bed_name)
                species_pep_path = os.path.join(request_dir, 'species.pep')
                species_bed_path = os.path.join(request_dir, 'species.bed')

                with open(uploaded_pep_path, 'wb') as f:
                    for chunk in pep_file.chunks():
                        f.write(chunk)
                with open(uploaded_bed_path, 'wb') as f:
                    for chunk in bed_file.chunks():
                        f.write(chunk)

                shutil.copyfile(uploaded_pep_path, renamed_pep_path)
                shutil.copyfile(uploaded_bed_path, renamed_bed_path)
                shutil.copyfile(renamed_pep_path, species_pep_path)
                shutil.copyfile(renamed_bed_path, species_bed_path)

                script_path = '/home/yunlongli/Annelida/Synteny/Annelida-upload.py'
                command = ['python3', script_path, species_slug]
                log_path = os.path.join(request_dir, 'run.log')
                status_path = os.path.join(request_dir, 'status.json')
                try:
                    with open(status_path, 'w', encoding='utf-8') as status_file:
                        json.dump({
                            'status': 'queued',
                            'pid': None,
                            'created': datetime.utcnow().isoformat(),
                            'species': species_raw,
                            'species_slug': species_slug,
                            'page_url': request.build_absolute_uri(f"{request.path}?token={upload_token}"),
                            'email_address': user_email,
                            'completion_email_sent': False,
                            'original_uploaded_files': [uploaded_pep_name, uploaded_bed_name],
                            'uploaded_files': [uploaded_pep_storage_name, uploaded_bed_storage_name],
                            'renamed_files': [renamed_pep_name, renamed_bed_name],
                        }, status_file)

                    _enqueue_job(
                        _run_synteny_upload_task,
                        request_dir,
                        command,
                        log_path,
                        status_path,
                        species_raw,
                        species_slug,
                        user_email,
                        request.build_absolute_uri(f"{request.path}?token={upload_token}"),
                    )

                    upload_info = {
                        'token': upload_token,
                        'base_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/"),
                        'page_url': request.build_absolute_uri(f"{request.path}?token={upload_token}"),
                        'files': [],
                        'status': 'running',
                        'log_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/run.log"),
                        'species': species_raw,
                        'uploaded_files': [uploaded_pep_name, uploaded_bed_name],
                        'renamed_files': [renamed_pep_name, renamed_bed_name],
                    }
                    if user_email:
                        subject = 'Your AnnelidaDB synteny job is queued'
                        message = (
                            f'Your synteny comparison request has been submitted successfully.\n\n'
                            f'You can track it at the following link:\n{upload_info["page_url"]}\n\n'
                            'Please keep this link for future reference.'
                        )
                        _send_notification_email(subject, message, user_email)
                    return redirect(upload_info['page_url'])
                except Exception as exc:
                    upload_error = f'Unable to start background synteny comparison: {str(exc)}'
                    upload_info = {
                        'token': upload_token,
                        'base_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/"),
                        'page_url': request.build_absolute_uri(f"{request.path}?token={upload_token}"),
                        'files': [],
                        'status': 'failed',
                        'log_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/run.log"),
                    }

    if upload_info and upload_info.get('files'):
        png_files = [file for file in upload_info['files'] if file['name'].lower().endswith('.png')]
        upload_info['png_files'] = png_files
        if png_files:
            upload_info['preview_url'] = png_files[0]['url']
            upload_info['preview_name'] = png_files[0]['name']

    return render(request, 'synteny_upload.html', {
        'active_page': 'tools',
        'page_title': 'Synteny Upload',
        'upload_info': upload_info,
        'upload_error': upload_error,
        'show_upload_info': show_upload_info,
        'species_raw': species_raw,
    })


def synteny_upload_edit_seqids_page(request):
    token = request.GET.get('token')
    if not token:
        return HttpResponseNotFound('Missing synteny upload token.')

    upload_root = _get_synteny_upload_root()
    request_dir = os.path.join(upload_root, token)
    if not os.path.isdir(request_dir):
        return HttpResponseNotFound('The requested synteny upload token was not found.')

    seqids_file_name = _find_oxford_plot_seqids_file(request_dir)
    if not seqids_file_name:
        return render(request, 'synteny_upload_edit_seqids.html', {
            'active_page': 'tools',
            'page_title': 'Edit Oxford Plot SeqIDs',
            'token': token,
            'error_message': 'No Annelida.*.oxford_plot.seqids file was found for this upload.',
        })

    seqids_path = os.path.join(request_dir, seqids_file_name)
    seqids_content = ''
    error_message = None

    if request.method == 'POST':
        seqids_content = request.POST.get('seqids_content', '')
        try:
            with open(seqids_path, 'w', encoding='utf-8') as seqids_file:
                seqids_file.write(seqids_content)

            status_path = os.path.join(request_dir, 'status.json')
            log_path = os.path.join(request_dir, 'run.log')
            script_path = '/home/yunlongli/Annelida/Synteny/Annelida-upload.py'
            species_slug = None
            species_raw_value = None
            if os.path.exists(status_path):
                try:
                    with open(status_path, 'r', encoding='utf-8') as status_file:
                        status_data = json.load(status_file)
                        species_slug = status_data.get('species_slug')
                        species_raw_value = status_data.get('species')
                except Exception:
                    species_slug = None
            if not species_slug:
                species_slug = os.path.splitext(seqids_file_name)[0]
            if not species_raw_value:
                species_raw_value = species_slug

            with open(status_path, 'w', encoding='utf-8') as status_file:
                json.dump({
                    'status': 'queued',
                    'pid': None,
                    'created': datetime.utcnow().isoformat(),
                    'species': species_raw_value,
                    'species_slug': species_slug,
                }, status_file)

            _enqueue_job(
                _run_synteny_upload_task,
                request_dir,
                ['python3', script_path, species_slug],
                log_path,
                status_path,
                species_raw_value,
                species_slug,
            )
            return redirect(request.build_absolute_uri(f"/synteny-upload/?token={token}"))
        except Exception as exc:
            error_message = f'Unable to save changes and rerun the task: {exc}'

    else:
        try:
            with open(seqids_path, 'r', encoding='utf-8') as seqids_file:
                seqids_content = seqids_file.read()
        except Exception as exc:
            error_message = f'Unable to read {seqids_file_name}: {exc}'

    return render(request, 'synteny_upload_edit_seqids.html', {
        'active_page': 'tools',
        'page_title': 'Edit Oxford Plot SeqIDs',
        'token': token,
        'seqids_file_name': seqids_file_name,
        'seqids_content': seqids_content,
        'error_message': error_message,
    })


def synteny_upload_download_selected(request):
    token = request.POST.get('token') if request.method == 'POST' else request.GET.get('token')
    if not token:
        return HttpResponseNotFound('Missing synteny upload token.')

    upload_root = _get_synteny_upload_root()
    request_dir = os.path.join(upload_root, token)
    if not os.path.isdir(request_dir):
        return HttpResponseNotFound('The requested synteny upload token was not found.')

    selected_files = request.POST.getlist('file') if request.method == 'POST' else request.GET.getlist('file')
    if not selected_files:
        return HttpResponseNotFound('No files selected for download.')

    valid_files = []
    for name in selected_files:
        if not name or '/' in name or '\\' in name:
            continue
        file_path = os.path.join(request_dir, name)
        if not os.path.isfile(file_path):
            continue
        valid_files.append((name, file_path))

    if not valid_files:
        return HttpResponseNotFound('Selected files were not found in the upload folder.')

    if len(valid_files) == 1:
        file_name, file_path = valid_files[0]
        response = FileResponse(open(file_path, 'rb'), as_attachment=True, filename=file_name)
        return response

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        for file_name, file_path in valid_files:
            zf.write(file_path, arcname=file_name)
    zip_buffer.seek(0)

    response = HttpResponse(zip_buffer.getvalue(), content_type='application/zip')
    file_name = f'synteny_upload_{token}.zip'
    response['Content-Disposition'] = (
        f'attachment; filename="{file_name}"; filename*=UTF-8''{urllib.parse.quote(file_name)}'
    )
    return response


GENE_FAMILY_SOURCE_DIRECTORIES = {
    'ANNEVO': 'ANNEVO',
    'Evidence-based': 'Evidence-based',
    'All': 'All',
}


def _normalize_gene_family_source(source):
    if not source:
        return 'ANNEVO'
    normalized = source.strip().lower()
    if normalized in ('annevo', 'anvevo'):
        return 'ANNEVO'
    if normalized in ('evidence-based', 'evidence_based', 'evidence based', 'evidence'):
        return 'Evidence-based'
    if normalized == 'all':
        return 'All'
    return 'ANNEVO'


def _gene_family_file_path(source):
    if source not in GENE_FAMILY_SOURCE_DIRECTORIES:
        return None
    path = os.path.join(settings.BASE_DIR, 'static', 'gene_family', GENE_FAMILY_SOURCE_DIRECTORIES[source], 'Orthogroups.tsv')
    return path if os.path.isfile(path) else None


def _format_gene_family_species_column(column_name):
    if not column_name:
        return ''
    first_part = column_name.split('.', 1)[0]
    cleaned = first_part.replace('_', ' ').strip()
    return cleaned or column_name


def _map_gene_family_species_to_taxonomy_ids(species_labels):
    taxonomy_map = {}
    queryset = Taxonomy.objects.filter(species__isnull=False)
    for tax in queryset:
        taxonomy_map[_format_species(tax.species)] = tax.ncbi_taxonomy_id

    for label in species_labels:
        if not label:
            continue
        formatted = _format_species(label)
        if formatted in taxonomy_map:
            continue
        tax = get_taxonomy_for_species_name(label)
        if tax:
            taxonomy_map[formatted] = tax.ncbi_taxonomy_id

    return taxonomy_map


def _parse_gene_family_orthogroups(orthogroups_path, source, max_families=None):
    families = []
    if not orthogroups_path or not os.path.isfile(orthogroups_path):
        return families

    try:
        with open(orthogroups_path, 'r', encoding='utf-8', errors='ignore') as f:
            header_line = f.readline()
            if not header_line:
                return families
            header = header_line.rstrip('\r\n').split('\t')
            species_columns = header[1:]
            species_labels = [_format_gene_family_species_column(col) for col in species_columns]

            taxonomy_map = _map_gene_family_species_to_taxonomy_ids(species_labels)
            for row_line in f:
                if not row_line.strip():
                    continue
                cols = row_line.rstrip('\r\n').split('\t')
                if len(cols) < 2:
                    continue
                family_id = cols[0].strip()
                if not family_id:
                    continue

                total_members = 0
                representative_gene = ''
                representative_species = ''
                representative_taxonomy_id = None
                species_with_members = 0

                for idx, cell in enumerate(cols[1:]):
                    cell_text = cell.strip()
                    if not cell_text:
                        continue
                    genes = [gene.strip() for gene in cell_text.split(',') if gene.strip()]
                    if not genes:
                        continue
                    species_with_members += 1
                    total_members += len(genes)
                    if not representative_gene:
                        representative_gene = genes[0]
                        if idx < len(species_labels):
                            representative_species = species_labels[idx]
                            representative_taxonomy_id = taxonomy_map.get(_format_species(species_labels[idx]))

                families.append({
                    'family_id': family_id,
                    'representative_gene': representative_gene,
                    'species': representative_species or 'Multiple species',
                    'taxonomy_id': representative_taxonomy_id,
                    'member_count': total_members,
                    'source': source,
                    'presented_species': f'{species_with_members} species',
                })

                if max_families and len(families) >= max_families:
                    break
    except Exception:
        pass
    return families


def _parse_gene_family_memberships_in_file(file_path, source_name, normalized_query):
    memberships = []
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            header_line = f.readline()
            if not header_line:
                return memberships

            for row_line in f:
                if not row_line.strip():
                    continue
                cols = row_line.rstrip('\r\n').split('\t')
                if len(cols) < 2:
                    continue

                family_id = cols[0].strip()
                if not family_id:
                    continue

                total_members = 0
                species_count = 0
                found = False
                for cell in cols[1:]:
                    cell_text = cell.strip()
                    if not cell_text:
                        continue
                    genes = [gene.strip() for gene in cell_text.split(',') if gene.strip()]
                    if not genes:
                        continue
                    species_count += 1
                    total_members += len(genes)
                    for gene in genes:
                        if gene.lower() == normalized_query:
                            found = True
                if found:
                    membership = {
                        'family_id': family_id,
                        'family_name': f'Orthogroup {family_id}',
                        'source': source_name,
                        'member_count': total_members,
                        'presented_species': f'{species_count} species',
                    }
                    memberships.append(membership)
    except Exception:
        pass
    return memberships


def _get_gene_family_families(source, max_families=None):
    selected_source = _normalize_gene_family_source(source)
    cache_key = f'gene_family_families:{selected_source}:{max_families or "all"}'
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    families = []
    if selected_source == 'All':
        all_file_path = _gene_family_file_path('All')
        if all_file_path:
            families = _parse_gene_family_orthogroups(all_file_path, 'All', max_families=max_families)
        else:
            for source_name in ['ANNEVO', 'Evidence-based']:
                file_path = _gene_family_file_path(source_name)
                families.extend(_parse_gene_family_orthogroups(file_path, source_name, max_families=max_families))
    else:
        file_path = _gene_family_file_path(selected_source)
        families = _parse_gene_family_orthogroups(file_path, selected_source, max_families=max_families)

    cache.set(cache_key, families, 3600)
    return families


def _find_gene_family_memberships(gene_id):
    if not gene_id:
        return []

    memberships = []
    normalized_query = gene_id.strip().lower()
    for source_name in ['ANNEVO', 'Evidence-based']:
        file_path = _gene_family_file_path(source_name)
        if not file_path:
            continue

        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                header_line = f.readline()
                if not header_line:
                    continue

                for row_line in f:
                    if not row_line.strip():
                        continue
                    cols = row_line.rstrip('\r\n').split('\t')
                    if len(cols) < 2:
                        continue

                    family_id = cols[0].strip()
                    if not family_id:
                        continue

                    total_members = 0
                    species_count = 0
                    found = False
                    for cell in cols[1:]:
                        cell_text = cell.strip()
                        if not cell_text:
                            continue
                        genes = [gene.strip() for gene in cell_text.split(',') if gene.strip()]
                        if not genes:
                            continue
                        species_count += 1
                        total_members += len(genes)
                        for gene in genes:
                            if gene.lower() == normalized_query:
                                found = True
                    if found:
                        memberships.append({
                            'family_id': family_id,
                            'family_name': f'Orthogroup {family_id}',
                            'source': source_name,
                            'member_count': total_members,
                            'presented_species': f'{species_count} species',
                        })
        except Exception:
            continue

    return memberships


def _find_gene_family_memberships_by_source(gene_id):
    cache_key = f"gene_family_memberships_by_source:{gene_id.strip().lower()}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    memberships = {
        'ANNEVO': [],
        'Evidence_based': [],
        'All': [],
    }
    if not gene_id:
        cache.set(cache_key, memberships, 3600)
        return memberships

    normalized_query = gene_id.strip().lower()
    all_file_path = _gene_family_file_path('All')
    all_memberships = []
    if all_file_path:
        all_memberships = _parse_gene_family_memberships_in_file(all_file_path, 'All', normalized_query)
        memberships['All'] = all_memberships

    for source_name in ['ANNEVO', 'Evidence-based']:
        file_path = _gene_family_file_path(source_name)
        if not file_path:
            continue
        source_memberships = _parse_gene_family_memberships_in_file(file_path, source_name, normalized_query)
        key_name = 'Evidence_based' if source_name == 'Evidence-based' else source_name
        memberships[key_name].extend(source_memberships)
        if not all_file_path:
            all_memberships.extend(source_memberships)

    if not all_file_path:
        memberships['All'] = all_memberships

    cache.set(cache_key, memberships, 3600)
    return memberships


def gene_family_page(request):
    selected_source = _normalize_gene_family_source(request.GET.get('source', 'ANNEVO'))
    source_options = [
        {'value': 'ANNEVO', 'label': 'ANNEVO'},
        {'value': 'Evidence-based', 'label': 'Evidence-based'},
        {'value': 'All', 'label': 'All'},
    ]
    families = _get_gene_family_families(selected_source, max_families=200)
    return render(request, 'gene_family.html', {
        'active_page': 'tools',
        'families': families,
        'selected_source': selected_source,
        'source_options': source_options,
        'gene_family_api_url': reverse('gene_family_api'),
        'gene_family_preview_limit': 200,
    })


def _get_gene_family_family_detail(source, family_id):
    selected_source = _normalize_gene_family_source(source)
    if selected_source == 'All':
        all_file_path = _gene_family_file_path('All')
        if all_file_path:
            search_sources = ['All']
        else:
            search_sources = ['ANNEVO', 'Evidence-based']
    else:
        search_sources = [selected_source]

    for source_name in search_sources:
        file_path = _gene_family_file_path(source_name)
        if not file_path:
            continue

        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                header_line = f.readline()
                if not header_line:
                    continue
                header = header_line.rstrip('\r\n').split('\t')
                species_columns = header[1:]
                species_labels = [_format_gene_family_species_column(col) for col in species_columns]
                taxonomy_map = _map_gene_family_species_to_taxonomy_ids(species_labels)

                for row_line in f:
                    if not row_line.strip():
                        continue
                    cols = row_line.rstrip('\r\n').split('\t')
                    if not cols or cols[0].strip() != family_id:
                        continue

                    family_members = []
                    total_members = 0
                    representative_gene = ''
                    representative_species = ''
                    representative_taxonomy_id = None

                    for idx, cell in enumerate(cols[1:]):
                        cell_text = cell.strip()
                        if not cell_text:
                            continue
                        genes = [gene.strip() for gene in cell_text.split(',') if gene.strip()]
                        if not genes:
                            continue
                        species_name = species_labels[idx] if idx < len(species_labels) else f'Species {idx + 1}'
                        taxonomy_id = taxonomy_map.get(_format_species(species_name))
                        if not representative_gene:
                            representative_gene = genes[0]
                            representative_species = species_name
                            representative_taxonomy_id = taxonomy_id

                        family_members.append({
                            'species': species_name,
                            'taxonomy_id': taxonomy_id,
                            'genes': genes,
                            'member_count': len(genes),
                        })
                        total_members += len(genes)

                    if total_members == 0:
                        continue

                    return {
                        'family_id': family_id,
                        'family_name': f'Orthogroup {family_id}',
                        'source': source_name,
                        'representative_gene': representative_gene,
                        'species': representative_species or 'Multiple species',
                        'taxonomy_id': representative_taxonomy_id,
                        'member_count': total_members,
                        'presented_species': f'{len(family_members)} species',
                        'description': f'{source_name} orthogroup with {total_members} gene members across {len(family_members)} species.',
                        'members': family_members,
                    }
        except Exception:
            continue
    return None


def gene_family_detail_page(request, family_id):
    selected_source = _normalize_gene_family_source(request.GET.get('source', 'ANNEVO'))
    detail = _get_gene_family_family_detail(selected_source, family_id)
    if not detail:
        return HttpResponseNotFound('Gene family not found.')

    source_options = [
        {'value': 'ANNEVO', 'label': 'ANNEVO'},
        {'value': 'Evidence-based', 'label': 'Evidence-based'},
        {'value': 'All', 'label': 'All'},
    ]
    return render(request, 'gene_family_detail.html', {
        'active_page': 'tools',
        'family': detail,
        'selected_source': selected_source,
        'source_options': source_options,
    })


@require_GET
def gene_family_api(request):
    selected_source = _normalize_gene_family_source(request.GET.get('source', 'ANNEVO'))
    families = _get_gene_family_families(selected_source)
    return JsonResponse({
        'status': 'success',
        'families': families,
    })


def gene_search_page(request):
    # Get species that have genome data
    species_list = Taxonomy.objects.filter(genomes__isnull=False).distinct().order_by('species')
    return render(request, 'gene_search.html', {
        'active_page': 'tools',
        'species_list': species_list
    })


def _ensure_gene_search_fts_table():
    from django.db import connection
    with connection.cursor() as cursor:
        raw_cursor = getattr(cursor, 'cursor', cursor)
        # Create the FTS table if it does not exist
        raw_cursor.execute('''
            CREATE VIRTUAL TABLE IF NOT EXISTS gene_search_fts USING fts5(
                description, preferred_name, species, taxonomy_id, gene_id
            )
        ''')
        # Populate or refresh the index
        raw_cursor.execute('DELETE FROM gene_search_fts')
        raw_cursor.execute('''
            INSERT INTO gene_search_fts(rowid, description, preferred_name, species, taxonomy_id, gene_id)
            SELECT ea.id, ea.description, ea.preferred_name, t.species, ea.taxonomy_id, ea.query
            FROM api_eggnogannotation ea
            LEFT JOIN api_taxonomy t ON ea.taxonomy_id = t.id
        ''')


def _fallback_gene_search(query, species_filter):
    from django.db import connection
    sql_query = '''
        SELECT ea.query, ea.description, ea.preferred_name, t.species, ea.taxonomy_id
        FROM api_eggnogannotation ea
        LEFT JOIN api_taxonomy t ON ea.taxonomy_id = t.id
        WHERE (ea.description LIKE ? OR ea.preferred_name LIKE ?)
    '''
    params = [f'%{query}%', f'%{query}%']
    if species_filter:
        sql_query += ' AND ea.taxonomy_id = ?'
        params.append(species_filter)
    sql_query += ' LIMIT 200'

    with connection.cursor() as cursor:
        raw_cursor = getattr(cursor, 'cursor', cursor)
        raw_cursor.execute(sql_query, params)
        return raw_cursor.fetchall()


def _parse_gff_attributes(attr_str):
    attrs = {}
    for part in attr_str.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            attrs[key.strip()] = urllib.parse.unquote(value.strip())
    return attrs


def _header_matches_gene_id(gene_id, header):
    if not header or not gene_id:
        return False
    header_text = header.strip().lower()
    normalized_header = re.sub(r'[^a-z0-9]+', ' ', header_text)
    normalized_query = re.sub(r'[^a-z0-9]+', ' ', gene_id.strip().lower())
    return normalized_query in normalized_header or gene_id.strip().lower() == header_text


def _gene_attributes_match(gene_id, attributes):
    if not gene_id or not attributes:
        return False
    gene_lower = gene_id.strip().lower()
    for value in attributes.values():
        if not value:
            continue
        if gene_lower == value.strip().lower() or gene_lower in value.strip().lower():
            return True
    return False


def _find_gene_location(tax, gene_id, annotation_source=None):
    cache_key = f"gene_location:{tax.ncbi_taxonomy_id}:{gene_id}:{annotation_source or ''}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    formatted_species = _format_species(tax.species)
    data_dir = _species_genome_dir(formatted_species)
    gff_patterns = ['*.gff3', '*.gff3.gz', '*.gff', '*.gtf', '*.gtf.gz']
    candidates = []
    prefer_annevo = annotation_source and 'ANNEVO' in annotation_source.upper()

    for pattern in gff_patterns:
        for gff_path in glob.glob(os.path.join(data_dir, pattern)):
            try:
                if gff_path.lower().endswith('.gz'):
                    gff_file_handle = gzip.open(gff_path, 'rt', encoding='utf-8', errors='ignore')
                else:
                    gff_file_handle = open(gff_path, 'r', encoding='utf-8', errors='ignore')
                with gff_file_handle as gff_file:
                    for line in gff_file:
                        if not line.strip() or line.startswith('#'):
                            continue
                        cols = line.strip().split('\t')
                        if len(cols) < 9:
                            continue
                        seqid, source, feature_type, start, end, score, strand, phase, attributes = cols[:9]
                        attr_map = _parse_gff_attributes(attributes)
                        if _gene_attributes_match(gene_id, attr_map):
                            candidates.append({
                                'file': os.path.basename(gff_path),
                                'seqid': seqid,
                                'source': source,
                                'feature_type': feature_type,
                                'start': int(start),
                                'end': int(end),
                                'strand': strand,
                                'attributes': attr_map,
                                'is_annevo': prefer_annevo and 'ANNEVO' in os.path.basename(gff_path).upper(),
                            })
            except OSError:
                continue

    if not candidates:
        cache.set(cache_key, None, 3600)
        return None

    priority = {'gene': 0, 'mRNA': 1, 'transcript': 2, 'CDS': 3, 'exon': 4}
    if prefer_annevo:
        candidates.sort(key=lambda item: (0 if item['is_annevo'] else 1, priority.get(item['feature_type'].lower(), 10)))
    else:
        candidates.sort(key=lambda item: priority.get(item['feature_type'].lower(), 10))
    best = candidates[0]
    matching_seqid = best['seqid']

    fasta_lengths = _get_contig_lengths(data_dir)
    seqid_max = fasta_lengths.get(matching_seqid)
    if not seqid_max:
        seqid_max = max(item['end'] for item in candidates if item['seqid'] == matching_seqid)

    best['contig_length'] = seqid_max
    best['relative_start'] = max(0.0, min(100.0, (best['start'] - 1) / seqid_max * 100)) if seqid_max else 0.0
    best['relative_width'] = max(0.8, min(100.0, (best['end'] - best['start'] + 1) / seqid_max * 100)) if seqid_max else 100.0
    cache.set(cache_key, best, 3600)
    return best


def _read_fasta_sequences(fasta_path):
    open_func = gzip.open if fasta_path.lower().endswith('.gz') else open
    with open_func(fasta_path, 'rt', encoding='utf-8', errors='ignore') as fasta_file:
        header = None
        seq_parts = []
        for line in fasta_file:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_parts)


def _count_fasta_records(fasta_path):
    if not fasta_path or not os.path.isfile(fasta_path):
        return 0

    index_count = _count_fasta_records_by_index(fasta_path)
    if index_count is not None:
        return index_count

    open_func = gzip.open if fasta_path.lower().endswith('.gz') else open
    try:
        with open_func(fasta_path, 'rt', encoding='utf-8', errors='ignore') as fasta_file:
            return sum(1 for line in fasta_file if line.startswith('>'))
    except OSError:
        return 0


def _count_fasta_records_by_index(fasta_path):
    if not fasta_path or not os.path.isfile(fasta_path):
        return None

    fasta_indexes = _find_fasta_indexes(fasta_path)
    index_path = fasta_indexes.get('fai')
    if not index_path:
        return None

    open_func = gzip.open if index_path.lower().endswith('.gz') else open
    try:
        with open_func(index_path, 'rt', encoding='utf-8', errors='ignore') as index_file:
            return sum(1 for line in index_file if line.strip())
    except Exception:
        return None


def _count_sequence_records_in_directory(data_dir, seq_type, recursive=True):
    if not os.path.isdir(data_dir):
        return 0
    total = 0
    if recursive:
        for root, _, filenames in os.walk(data_dir):
            for filename in sorted(filenames):
                if _is_sequence_file_name(filename, seq_type):
                    total += _count_fasta_records(os.path.join(root, filename))
    else:
        for filename in sorted(os.listdir(data_dir)):
            if _is_sequence_file_name(filename, seq_type):
                total += _count_fasta_records(os.path.join(data_dir, filename))
    return total


def _normalize_genome_filename(filename):
    if not filename:
        return filename
    return filename.replace('\r', '').replace('\n', '').strip()


def _sanitize_genome_file_names(data_dir):
    if not os.path.isdir(data_dir):
        return
    for filename in os.listdir(data_dir):
        normalized_name = _normalize_genome_filename(filename)
        if normalized_name == filename:
            continue
        if not normalized_name:
            continue
        src = os.path.join(data_dir, filename)
        dst = os.path.join(data_dir, normalized_name)
        if os.path.exists(dst) and os.path.abspath(src) != os.path.abspath(dst):
            base, ext = os.path.splitext(normalized_name)
            index = 1
            while True:
                candidate = f'{base}_{index}{ext}'
                target = os.path.join(data_dir, candidate)
                if not os.path.exists(target):
                    dst = target
                    break
                index += 1
        try:
            os.rename(src, dst)
        except OSError:
            continue


def _is_sequence_file_name(filename, seq_type):
    if not filename:
        return False
    name = filename.lower()
    if 'longest' in name:
        return False
    if seq_type == 'protein':
        return bool(re.search(r'\.pep\.fa(?:sta)?(?:\.gz)?$|\.pep\.faa(?:\.gz)?$|\.protein\.fa(?:sta)?(?:\.gz)?$|\.protein\.faa(?:\.gz)?$|\.faa(?:\.gz)?$', name))
    if seq_type == 'cds':
        return bool(re.search(r'\.cds\.fa(?:sta)?(?:\.gz)?$', name))
    return False


def _find_sequence_file_paths(data_dir, seq_type):
    file_paths = []
    if not os.path.isdir(data_dir):
        return file_paths
    _sanitize_genome_file_names(data_dir)
    for filename in sorted(os.listdir(data_dir)):
        path = os.path.join(data_dir, filename)
        if not os.path.isfile(path):
            continue
        if _is_sequence_file_name(filename, seq_type):
            file_paths.append({'name': filename, 'path': path})
    return file_paths


def _list_sequence_files(data_dir, seq_type):
    return [item['name'] for item in _find_sequence_file_paths(data_dir, seq_type)]


def _find_sequence_in_fasta_files(gene_id, fasta_paths):
    for fasta_path in fasta_paths:
        if not fasta_path or not os.path.isfile(fasta_path):
            continue
        header, seq = _find_sequence_in_fasta_file(fasta_path, gene_id)
        if seq:
            return header, seq, fasta_path
    return None, None, None


def _find_sequence_in_species_data(seq_type, gene_id, tax=None, species=None):
    if not gene_id:
        return None

    search_paths = []
    if isinstance(tax, Taxonomy):
        formatted_species = _format_species(tax.species)
    else:
        formatted_species = _format_species(species or '')

    if formatted_species:
        genome_dir = _species_genome_dir(formatted_species)
        if genome_dir and os.path.isdir(genome_dir):
            for root, _, files in os.walk(genome_dir):
                for filename in files:
                    if _is_sequence_file_name(filename, seq_type):
                        search_paths.append(os.path.join(root, filename))

    if not search_paths:
        data_root = os.path.join(settings.BASE_DIR, 'static', 'data')
        if os.path.isdir(data_root):
            for species_dir in sorted(os.listdir(data_root)):
                genome_dir = os.path.join(data_root, species_dir, 'genome')
                if not os.path.isdir(genome_dir):
                    continue
                for root, _, files in os.walk(genome_dir):
                    for filename in files:
                        if _is_sequence_file_name(filename, seq_type):
                            search_paths.append(os.path.join(root, filename))

    _, seq, _ = _find_sequence_in_fasta_files(gene_id, search_paths)
    return seq


def _read_sequence_from_indexed_fasta(fasta_path, seq_id):
    try:
        result = subprocess.run(['samtools', 'faidx', fasta_path, seq_id], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError:
        return None, None

    lines = result.stdout.strip().splitlines()
    if not lines or not lines[0].startswith('>'):
        return None, None

    header = lines[0][1:].strip()
    sequence = ''.join(line.strip() for line in lines[1:] if line.strip())
    return header, sequence


def _find_sequence_in_fasta_file(fasta_path, gene_id):
    header, seq = _read_sequence_from_indexed_fasta(fasta_path, gene_id)
    if seq:
        return header, seq

    for fasta_header, fasta_seq in _read_fasta_sequences(fasta_path):
        if fasta_header.split()[0].strip().lower() == gene_id.strip().lower():
            return fasta_header, fasta_seq
    return None, None


def _find_gene_sequences(tax, gene_id, selected_files=None, data_dir=None):
    cache_key = f"gene_sequences:{getattr(tax, 'ncbi_taxonomy_id', '')}:{gene_id}"
    if selected_files:
        cache_key += ':' + ':'.join(f"{seq_type}={os.path.basename(path)}" for seq_type, path in sorted(selected_files.items()))

    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    found = {}
    selected_files = selected_files or {}

    def _search_fallback(seq_type, selected_file):
        paths = []
        if selected_file and os.path.isfile(selected_file):
            paths.append(selected_file)
        if data_dir and os.path.isdir(data_dir):
            paths.extend([item['path'] for item in _find_sequence_file_paths(data_dir, seq_type) if item['path'] not in paths])
        header, seq, fasta_path = _find_sequence_in_fasta_files(gene_id, paths)
        if seq:
            found[seq_type] = {
                'sequence': seq,
                'header': header,
                'file_name': os.path.basename(fasta_path),
                'format': 'protein' if seq_type == 'protein' else 'nucleotide',
                'type': seq_type,
            }

    if isinstance(tax, Taxonomy):
        for seq_type in ['protein', 'cds']:
            selected_file = selected_files.get(seq_type)
            query = SequenceFileIndex.objects.filter(
                taxonomy=tax,
                seq_type=seq_type,
                seq_id__iexact=gene_id,
            )
            if selected_file:
                query = query.filter(file_name__iexact=os.path.basename(selected_file))
            entry = query.order_by('id').first()
            if entry:
                fasta_path = _resolve_indexed_fasta_path(entry.file_path)
                if fasta_path:
                    header, seq = _find_sequence_in_fasta_file(fasta_path, gene_id)
                    if seq:
                        found[seq_type] = {
                            'sequence': seq,
                            'header': header,
                            'file_name': entry.file_name,
                            'format': 'protein' if seq_type == 'protein' else 'nucleotide',
                            'type': seq_type,
                        }
            if seq_type not in found:
                _search_fallback(seq_type, selected_file)
    else:
        if data_dir and os.path.isdir(data_dir):
            for seq_type in ['protein', 'cds']:
                _search_fallback(seq_type, selected_files.get(seq_type))
        else:
            for seq_type, selected_file in selected_files.items():
                if not selected_file or not os.path.isfile(selected_file):
                    continue
                header, seq = _find_sequence_in_fasta_file(selected_file, gene_id)
                if seq:
                    found[seq_type] = {
                        'sequence': seq,
                        'header': header,
                        'file_name': os.path.basename(selected_file),
                        'format': 'protein' if seq_type == 'protein' else 'nucleotide',
                        'type': seq_type,
                    }

    cache.set(cache_key, found, 3600)
    return found


def _find_jbrowse_annotation_path(tax):
    formatted_species = _format_species(tax.species)
    data_dir = _species_genome_dir(formatted_species)
    if not data_dir or not os.path.isdir(data_dir):
        return None

    annotation_file = _find_preferred_annotation_file(data_dir)
    if not annotation_file:
        return None

    if annotation_file.lower().endswith('.gtf.gz'):
        return {
            'adapter': 'GtfAdapter',
            'uri': _format_static_url(annotation_file),
            'gtfGzLocation': {'uri': _format_static_url(annotation_file), 'locationType': 'UriLocation'},
        }

    if annotation_file.lower().endswith('.gtf'):
        return {
            'adapter': 'GtfAdapter',
            'uri': _format_static_url(annotation_file),
            'gtfLocation': {'uri': _format_static_url(annotation_file), 'locationType': 'UriLocation'},
        }

    if annotation_file.lower().endswith('.gff3.gz'):
        index_path = _find_annotation_index(annotation_file)
        use_tabix = index_path and _is_bgzip_file(annotation_file)
        if index_path and use_tabix:
            return {
                'adapter': 'Gff3TabixAdapter',
                'uri': _format_static_url(annotation_file),
                'gffGzLocation': {'uri': _format_static_url(annotation_file), 'locationType': 'UriLocation'},
                'index': {'location': {'uri': _format_static_url(index_path), 'locationType': 'UriLocation'}},
            }
        return {
            'adapter': 'Gff3Adapter',
            'uri': _format_static_url(annotation_file),
            'gffGzLocation': {'uri': _format_static_url(annotation_file), 'locationType': 'UriLocation'},
        }

    if annotation_file.lower().endswith('.gff3') or annotation_file.lower().endswith('.gff'):
        return {
            'adapter': 'Gff3Adapter',
            'uri': _format_static_url(annotation_file),
            'gffLocation': _format_static_url(annotation_file),
        }

    return {
        'adapter': 'Gff3Adapter',
        'uri': _format_static_url(annotation_file),
        'gffLocation': _format_static_url(annotation_file),
    }


def _format_sequence_lines(sequence, width=80):
    if not sequence:
        return ''
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))


def _split_terms(value):
    return [term.strip() for term in re.split(r'[;,\n]+', value) if term.strip()]


def _link_go_term(term):
    term = term.strip()
    if not term:
        return ''
    match = re.match(r'^(GO:\d{7})(.*)$', term)
    if match:
        go_id = match.group(1)
        description = match.group(2).strip()
        href = f'https://www.ebi.ac.uk/QuickGO/term/{go_id}'
        if description:
            return format_html('<a href="{}" target="_blank" rel="noopener">{}</a> {}', href, go_id, description)
        return format_html('<a href="{}" target="_blank" rel="noopener">{}</a>', href, go_id)
    return format_html('<span>{}</span>', term)


def _link_kegg_term(term):
    term = term.strip()
    if not term:
        return ''
    match = re.match(r'^(?:ko:)?(K\d{5}|map\d{5}|RC\d{5}|R\d{5}|r\d{5})(.*)$', term, flags=re.IGNORECASE)
    if match:
        code = match.group(1)
        description = match.group(2).strip()
        code_lower = code.lower()
        code_upper = code.upper()
        if code_lower.startswith('k'):
            href = f'https://www.kegg.jp/dbget-bin/www_bget?ko:{code_upper}'
        elif code_lower.startswith('map'):
            href = f'https://www.kegg.jp/dbget-bin/www_bget?pathway:{code_lower}'
        elif code_upper.startswith('RC'):
            href = f'https://www.kegg.jp/dbget-bin/www_bget?rc:{code_upper}'
        elif code_lower.startswith('r'):
            href = f'https://www.kegg.jp/dbget-bin/www_bget?reaction:{code_upper}'
        else:
            href = None
        if href:
            if description:
                return format_html('<a href="{}" target="_blank" rel="noopener">{}</a> {}', href, code_upper, description)
            return format_html('<a href="{}" target="_blank" rel="noopener">{}</a>', href, code_upper)
    return format_html('<span>{}</span>', term)


def _link_ec_term(term):
    term = term.strip()
    if not term:
        return ''
    match = re.match(r'^(\d+\.\d+\.\d+\.\d+)(.*)$', term)
    if match:
        ec_id = match.group(1)
        description = match.group(2).strip()
        href = f'https://enzyme.expasy.org/EC/{ec_id}'
        if description:
            return format_html('<a href="{}" target="_blank" rel="noopener">{}</a> {}', href, ec_id, description)
        return format_html('<a href="{}" target="_blank" rel="noopener">{}</a>', href, ec_id)
    return format_html('<span>{}</span>', term)


def _link_brite_term(term):
    term = term.strip()
    if not term:
        return ''
    match = re.match(r'^(?:BRITE:?\s*)?([^;\s]+)(.*)$', term, flags=re.IGNORECASE)
    if match:
        label = match.group(1).strip()
        description = match.group(2).strip()
        href = f'https://www.kegg.jp/dbget-bin/www_bget?brite:{urllib.parse.quote(label)}'
        if description:
            return format_html('<a href="{}" target="_blank" rel="noopener">{}</a> {}', href, label, description)
        return format_html('<a href="{}" target="_blank" rel="noopener">{}</a>', href, label)
    return format_html('<span>{}</span>', term)


def _link_terms(value, label=''):
    terms = _split_terms(value)
    linked = []
    for term in terms:
        if label == 'GOs':
            linked.append(_link_go_term(term))
        elif label == 'EC number':
            linked.append(_link_ec_term(term))
        elif label.startswith('KEGG'):
            linked.append(_link_kegg_term(term))
        elif label == 'BRITE':
            linked.append(_link_brite_term(term))
        else:
            linked.append(format_html('<span>{}</span>', term))
    return mark_safe(', '.join(str(item) for item in linked))


def _format_go_terms(value):
    terms = _split_terms(value)
    linked_terms = [_link_go_term(term) for term in terms]
    if len(linked_terms) > 12:
        element_id = f'go-full-{uuid.uuid4().hex}'
        summary = mark_safe(', '.join(str(item) for item in linked_terms[:12]))
        full_list = mark_safe(', '.join(str(item) for item in linked_terms))
        return format_html(
            '<div class="annotation-collapsible"><span class="annotation-summary">{}</span>'
            '<button type="button" class="btn btn-link btn-sm p-0 ms-2" onclick="toggleCollapsibleList(\'{}\', this)">Show all</button>'
            '<div id="{}" class="annotation-full d-none mt-2">{}</div></div>',
            summary,
            element_id,
            element_id,
            full_list,
        )
    return mark_safe(', '.join(str(item) for item in linked_terms))


def _format_kegg_terms(value, label):
    terms = _split_terms(value)
    if label == 'KEGG pathway':
        map_terms = [term for term in terms if re.match(r'^(?:map\d{5})', term, flags=re.IGNORECASE)]
        if map_terms:
            terms = map_terms
    # Preserve order while removing exact duplicates.
    seen = set()
    unique_terms = []
    for term in terms:
        normalized = term.strip().lower()
        if normalized and normalized not in seen:
            seen.add(normalized)
            unique_terms.append(term)
    linked_terms = [_link_kegg_term(term) for term in unique_terms]
    if len(linked_terms) > 12:
        element_id = f'kegg-full-{uuid.uuid4().hex}'
        summary = mark_safe(', '.join(str(item) for item in linked_terms[:12]))
        full_list = mark_safe(', '.join(str(item) for item in linked_terms))
        return format_html(
            '<div class="annotation-collapsible"><span class="annotation-summary">{}</span>'
            '<button type="button" class="btn btn-link btn-sm p-0 ms-2" onclick="toggleCollapsibleList(\'{}\', this)">Show all</button>'
            '<div id="{}" class="annotation-full d-none mt-2">{}</div></div>',
            summary,
            element_id,
            element_id,
            full_list,
        )
    return mark_safe(', '.join(str(item) for item in linked_terms))


def _format_annotation_value(label, value):
    if not value:
        return value
    if label == 'E-value':
        try:
            return '{:.1e}'.format(float(value))
        except (ValueError, TypeError):
            return value
    if label == 'GOs':
        return _format_go_terms(value)
    if label == 'EC number':
        return _link_terms(value, label)
    if label.startswith('KEGG'):
        return _format_kegg_terms(value, label)
    if label == 'BRITE':
        return _link_terms(value, label)
    return value


def _collect_annotation_fields(annotation):
    fields = [
        ('Preferred name', annotation.preferred_name),
        ('Description', annotation.description),
        ('Seed ortholog', annotation.seed_ortholog),
        ('E-value', annotation.evalue),
        ('Score', annotation.score),
        ('EggNOG OGs', annotation.eggnog_ogs),
        ('Max annotation level', annotation.max_annot_lvl),
        ('COG category', annotation.cog_category),
        ('GOs', annotation.gos),
        ('EC number', annotation.ec),
        ('KEGG KO', annotation.kegg_ko),
        ('KEGG pathway', annotation.kegg_pathway),
        ('KEGG module', annotation.kegg_module),
        ('KEGG reaction', annotation.kegg_reaction),
        ('KEGG RClass', annotation.kegg_rclass),
        ('BRITE', annotation.brite),
        ('KEGG TC', annotation.kegg_tc),
        ('CAZy', annotation.cazys),
        ('BiGG reaction', annotation.bigg_reaction),
        ('PFAMs', annotation.pfams),
    ]
    formatted = []
    for label, value in fields:
        formatted_value = _format_annotation_value(label, value)
        if formatted_value:
            formatted.append((label, formatted_value))
    return formatted


@csrf_exempt
def gene_search_api(request):
    if request.method != 'POST':
        return JsonResponse({'status': 'error', 'message': 'POST only'})

    try:
        data = json.loads(request.body)
        query = data.get('query', '').strip()
        species_filter = data.get('species', '').strip()

        if not query:
            return JsonResponse({'status': 'error', 'message': 'Missing search query.'})

        # Create cache key
        cache_key = f"gene_search:{query}:{species_filter}"

        # Try to get from cache first
        cached_result = cache.get(cache_key)
        if cached_result:
            return JsonResponse(cached_result)

        # Use SQLite FTS (Full Text Search) for blazing fast search
        from django.db import connection

        # Build FTS query - search in description and preferred_name
        fts_query = f"description:{query} OR preferred_name:{query}"

        sql_query = """
            SELECT gene_id, description, preferred_name, species, taxonomy_id
            FROM gene_search_fts
            WHERE gene_search_fts MATCH ?
        """

        params = [fts_query]

        if species_filter:
            # For species filtering, we need to filter by taxonomy_id
            sql_query += " AND taxonomy_id = ?"
            params.append(species_filter)

        sql_query += " LIMIT 200"

        try:
            with connection.cursor() as cursor:
                raw_cursor = getattr(cursor, 'cursor', cursor)
                raw_cursor.execute(sql_query, params)
                rows = raw_cursor.fetchall()
        except Exception as exc:
            if 'no such table: gene_search_fts' in str(exc).lower():
                _ensure_gene_search_fts_table()
                with connection.cursor() as cursor:
                    raw_cursor = getattr(cursor, 'cursor', cursor)
                    raw_cursor.execute(sql_query, params)
                    rows = raw_cursor.fetchall()
            else:
                # Fall back to basic LIKE search if FTS query fails for any reason
                rows = _fallback_gene_search(query, species_filter)

        # Process results
        results = []
        for row in rows:
            gene_id, description, preferred_name, species, taxonomy_id = row
            results.append({
                'gene_id': gene_id,
                'description': description or preferred_name or 'No description',
                'species': species or 'Unknown',
                'taxonomy_id': taxonomy_id,
            })

        response_data = {
            'status': 'success',
            'query': query,
            'results': results,
            'count': len(results),
            'search_method': 'fts_sqlite'
        }

        # Cache the result for 15 minutes
        cache.set(cache_key, response_data, 900)

        return JsonResponse(response_data)

    except Exception as exc:
        return JsonResponse({'status': 'error', 'message': str(exc)})


def gene_detail_page(request, gene_id, taxonomy_id=None):
    annotation = None
    tax = None
    genome = None
    annotation_fields = []
    annotation_source = None

    if taxonomy_id:
        tax = get_object_or_404(Taxonomy, ncbi_taxonomy_id=taxonomy_id)
        try:
            annotation = EggnogAnnotation.objects.get(query=gene_id, taxonomy=tax)
        except EggnogAnnotation.DoesNotExist:
            annotation = {}
    else:
        try:
            annotation = EggnogAnnotation.objects.get(query=gene_id)
            tax = annotation.taxonomy
        except EggnogAnnotation.DoesNotExist:
            tax = _find_taxonomy_from_contig(gene_id)
            if tax is None:
                tax = _find_taxonomy_by_gene_in_annotation_files(gene_id)
            annotation = {}

    if tax is None:
        return HttpResponseNotFound('Gene not found or taxonomy information unavailable.')

    genome = tax.genomes.first() if hasattr(tax, 'genomes') else None
    formatted_species = _format_species(tax.species)
    data_dir = _species_genome_dir(formatted_species)
    _sanitize_genome_file_names(data_dir)

    protein_files = _list_sequence_files(data_dir, 'protein')
    cds_files = _list_sequence_files(data_dir, 'cds')

    annotation_source = None
    annotation_fields = []
    if annotation and isinstance(annotation, EggnogAnnotation):
        annotation_description = annotation.description or 'No description available.'
    else:
        annotation_description = annotation.get('description') if isinstance(annotation, dict) else None
        annotation_description = annotation_description or 'Loading annotation...'

    gene_family_memberships = {
        'All': [],
        'ANNEVO': [],
        'Evidence_based': [],
    }

    location = _find_gene_location(tax, gene_id, annotation_source)
    selected_protein_file, selected_cds_file = _auto_select_sequence_files(data_dir, annotation_source)
    selected_files = {}
    if selected_protein_file:
        selected_files['protein'] = os.path.join(data_dir, selected_protein_file)
    if selected_cds_file:
        selected_files['cds'] = os.path.join(data_dir, selected_cds_file)

    sequence_data = _find_gene_sequences(tax, gene_id, selected_files=selected_files, data_dir=data_dir)
    sequence_entries = []
    for seq_type in ['protein', 'cds']:
        seq_record = sequence_data.get(seq_type)
        if seq_record:
            sequence_entries.append({
                'label': 'Protein' if seq_type == 'protein' else 'CDS',
                'type': seq_type,
                'file_name': seq_record['file_name'],
                'format': seq_record['format'],
                'sequence_display': _format_sequence_lines(seq_record['sequence']),
                'header': seq_record['header'],
                'raw_sequence': seq_record['sequence'],
            })

    annotation_interpro_headers, annotation_interpro_rows = _find_interproscan_rows_for_gene(tax, gene_id)
    species_dir = _format_species(tax.species)
    jbrowse_data_path = _species_genome_url(species_dir)
    annotation_files = _find_annotation_files(_species_genome_dir(species_dir))
    jbrowse_track_ids = ','.join([f'annotationTrack{i+1}' for i in range(len(annotation_files))]) if annotation_files else 'annotationTrack'
    jbrowse_annotation = _find_jbrowse_annotation_path(tax)
    jbrowse_loc = None
    if location:
        padding = 5000
        start = max(1, location['start'] - padding)
        end = location['end'] + padding
        jbrowse_loc = f"{location['seqid']}:{start}..{end}"

    context = {
        'page_title': f'{gene_id} - Gene Detail',
        'annotation': annotation,
        'annotation_description': annotation_description,
        'tax': tax,
        'genome': genome,
        'gene_id': gene_id,
        'location': location,
        'sequence_entries': sequence_entries,
        'annotation_fields': annotation_fields,
        'annotation_interpro_headers': [],
        'annotation_interpro_rows': [],
        'species_dir': species_dir,
        'jbrowse_data_path': jbrowse_data_path,
        'jbrowse_annotation': jbrowse_annotation,
        'jbrowse_loc': jbrowse_loc,
        'jbrowse_track_ids': jbrowse_track_ids,
        'protein_sequence_files': protein_files,
        'cds_sequence_files': cds_files,
        'selected_protein_file': selected_protein_file,
        'selected_cds_file': selected_cds_file,
        'gene_detail_extras_url': reverse('gene_detail_extras'),
    }
    return render(request, 'gene_detail.html', context)


def _resolve_jbrowse_data_path(data_param):
    if not data_param:
        return None
    normalized = data_param.lstrip('/')
    full_path = os.path.normpath(os.path.join(settings.BASE_DIR, normalized))
    data_root = os.path.normpath(os.path.join(settings.BASE_DIR, 'static', 'data'))
    try:
        if os.path.commonpath([full_path, data_root]) != data_root:
            return None
    except ValueError:
        return None
    if not os.path.isdir(full_path):
        return None
    return full_path


def _find_first_file(root_dir, patterns):
    for pattern in patterns:
        for path in sorted(glob.glob(os.path.join(root_dir, pattern))):
            if os.path.isfile(path):
                return path
    return None


def _find_preferred_fasta_file(root_dir):
    patterns = [
        '*.chr.genomic.fasta.masked.gz',
        '*.genomic.fasta.masked.gz',
        '*.chr.genomic.fasta.gz',
        '*.genomic.fasta.gz',
        '*.fa.gz',
        '*.fasta.gz',
        '*.chr.genomic.fasta.masked',
        '*.genomic.fasta.masked',
        '*.chr.genomic.fasta',
        '*.genomic.fasta',
        '*.fa',
        '*.fasta',
    ]
    return _find_first_file(root_dir, patterns)


def _get_contig_lengths(data_dir):
    fasta_file = _find_preferred_fasta_file(data_dir)
    if not fasta_file:
        return {}

    fasta_indexes = _find_fasta_indexes(fasta_file)
    fasta_index = fasta_indexes.get('fai')
    lengths = {}
    if fasta_index and os.path.exists(fasta_index):
        try:
            with open(fasta_index, 'r', encoding='utf-8', errors='replace') as handle:
                for line in handle:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        lengths[parts[0]] = int(parts[1])
            if lengths:
                return lengths
        except Exception:
            pass

    try:
        current_seq = None
        current_length = 0
        open_func = gzip.open if fasta_file.lower().endswith('.gz') else open
        with open_func(fasta_file, 'rt', encoding='utf-8', errors='ignore') as handle:
            for line in handle:
                if line.startswith('>'):
                    if current_seq:
                        lengths[current_seq] = current_length
                    current_seq = line[1:].strip().split()[0]
                    current_length = 0
                else:
                    current_length += len(line.strip())
            if current_seq:
                lengths[current_seq] = current_length
    except Exception:
        pass

    return lengths


def _find_annotation_files(root_dir):
    preferred_patterns = [
        '*.chr.gtf',
        '*.chr.gtf.gz',
        '*.chr.ANNEVO.gff3',
    ]
    fallback_patterns = [
        '*ANNEVO*.gff3.gz',
        '*ANNEVO*.gff3',
        '*ANNEVO*.gtf.gz',
        '*ANNEVO*.gtf',
        '*.gff3.gz',
        '*.gff3',
        '*.gtf.gz',
        '*.gtf',
        '*.gff',
    ]
    seen = set()
    files = []

    for pattern in preferred_patterns:
        for path in sorted(glob.glob(os.path.join(root_dir, pattern))):
            if os.path.isfile(path) and path not in seen:
                seen.add(path)
                files.append(path)
    if files:
        return files

    for pattern in fallback_patterns:
        for path in sorted(glob.glob(os.path.join(root_dir, pattern))):
            if os.path.isfile(path) and path not in seen:
                seen.add(path)
                files.append(path)
    return files


def _find_preferred_annotation_file(root_dir):
    files = _find_annotation_files(root_dir)
    return files[0] if files else None


def _has_jbrowse_data_directory(species_dir):
    if not species_dir:
        return False
    return os.path.isdir(os.path.join(settings.BASE_DIR, 'static', 'data', species_dir, 'genome'))


def _has_epigenome_data_directory(species_dir):
    if not species_dir:
        return False
    return os.path.isdir(os.path.join(settings.BASE_DIR, 'static', 'data', species_dir, 'epigenome'))


def _find_epigenome_files(root_dir):
    if not root_dir or not os.path.isdir(root_dir):
        return []
    patterns = ['*.bigwig', '*.bw']
    files = []
    for pattern in patterns:
        for path in sorted(glob.glob(os.path.join(root_dir, pattern))):
            if os.path.isfile(path):
                files.append(path)
    return files


def _get_annotation_track_name(annotation_file):
    if not annotation_file:
        return 'Gene annotations'
    if '.annevo.' in os.path.basename(annotation_file).lower():
        return 'ANNEVO gene annotations'
    return 'Gene annotations'


def _is_bgzip_file(path):
    try:
        with open(path, 'rb') as f:
            header = f.read(18)
        if len(header) < 18 or header[0:2] != b'\x1f\x8b' or header[2] != 8:
            return False
        flg = header[3]
        if not (flg & 4):
            return False
        xlen = int.from_bytes(header[10:12], 'little')
        if len(header) < 12 + xlen:
            with open(path, 'rb') as f:
                header = f.read(12 + xlen)
            if len(header) < 12 + xlen:
                return False
        extra = header[12:12 + xlen]
        i = 0
        while i + 4 <= len(extra):
            si1 = extra[i]
            si2 = extra[i + 1]
            slen = int.from_bytes(extra[i + 2:i + 4], 'little')
            if i + 4 + slen > len(extra):
                break
            if si1 == 66 and si2 == 67:
                return True
            i += 4 + slen
        return False
    except Exception:
        return False


def _find_annotation_index(annotation_path):
    index_candidates = [annotation_path + '.tbi', annotation_path + '.csi']
    if annotation_path.endswith('.gz'):
        base = annotation_path[:-3]
        index_candidates.extend([base + '.tbi', base + '.csi'])

    for index_path in index_candidates:
        if os.path.exists(index_path):
            return index_path

    parent = os.path.dirname(annotation_path)
    name = os.path.basename(annotation_path)
    prefix = None
    if name.lower().endswith('.gff3.gz'):
        prefix = name[:-8]
    elif name.lower().endswith('.gtf.gz'):
        prefix = name[:-7]
    elif name.lower().endswith('.gff3'):
        prefix = name[:-5]
    elif name.lower().endswith('.gtf'):
        prefix = name[:-4]
    elif name.lower().endswith('.gff'):
        prefix = name[:-4]
    else:
        prefix = os.path.splitext(name)[0]

    for entry in os.listdir(parent):
        if entry.endswith('.tbi') or entry.endswith('.csi'):
            if entry.startswith(prefix):
                return os.path.join(parent, entry)
    return None


def _find_fasta_index(fasta_path):
    if fasta_path.endswith('.gz'):
        base = fasta_path[:-3]
        index_candidates = [fasta_path + '.gzi', fasta_path + '.gzi.gz', base + '.gzi', base + '.gzi.gz', fasta_path + '.fai', fasta_path + '.fai.gz', base + '.fai', base + '.fai.gz']
    else:
        index_candidates = [fasta_path + '.gzi', fasta_path + '.gzi.gz', fasta_path + '.fai', fasta_path + '.fai.gz']
    for index_path in index_candidates:
        if os.path.exists(index_path):
            return index_path
    return None


def _find_fasta_indexes(fasta_path):
    indexes = {'fai': None, 'gzi': None}
    if fasta_path.endswith('.gz'):
        base = fasta_path[:-3]
        fai_candidates = [fasta_path + '.fai', fasta_path + '.fai.gz', base + '.fai', base + '.fai.gz']
        gzi_candidates = [fasta_path + '.gzi', fasta_path + '.gzi.gz', base + '.gzi', base + '.gzi.gz']
    else:
        fai_candidates = [fasta_path + '.fai', fasta_path + '.fai.gz']
        gzi_candidates = [fasta_path + '.gzi', fasta_path + '.gzi.gz']

    for index_path in fai_candidates:
        if os.path.exists(index_path):
            indexes['fai'] = index_path
            break
    for index_path in gzi_candidates:
        if os.path.exists(index_path):
            indexes['gzi'] = index_path
            break
    return indexes


def _find_default_jbrowse_loc(data_dir):
    fasta_file = _find_preferred_fasta_file(data_dir)
    if not fasta_file:
        return None

    fasta_indexes = _find_fasta_indexes(fasta_file)
    fasta_index = fasta_indexes.get('fai')
    if fasta_index and os.path.exists(fasta_index):
        try:
            with open(fasta_index, 'r', encoding='utf-8', errors='replace') as handle:
                first_line = handle.readline().strip()
            if first_line:
                seqid, length_str, *_ = first_line.split('\t')
                length = int(length_str)
                window_end = min(length, 10000)
                return f"{seqid}:1..{window_end}"
        except Exception:
            pass

    try:
        open_func = gzip.open if fasta_file.lower().endswith('.gz') else open
        with open_func(fasta_file, 'rt', encoding='utf-8', errors='replace') as handle:
            for line in handle:
                if line.startswith('>'):
                    seqid = line[1:].strip().split()[0]
                    return f"{seqid}:1..10000"
    except Exception:
        pass

    return None


def _format_static_url(file_path):
    if not file_path:
        return None

    file_path = os.path.normpath(str(file_path))
    base_dir = os.path.normpath(str(settings.BASE_DIR))
    static_root = os.path.normpath(os.path.join(base_dir, 'static'))

    if os.path.isabs(file_path):
        abs_path = file_path
    else:
        abs_path = os.path.normpath(os.path.join(base_dir, file_path))

    normalized = abs_path.replace('\\', '/')
    static_index = normalized.find('/static/')
    if static_index != -1:
        return normalized[static_index:]

    if abs_path.startswith(static_root + os.path.sep) or abs_path == static_root:
        rel_path = os.path.relpath(abs_path, static_root)
        return '/static/' + rel_path.replace(os.path.sep, '/')

    return '/' + os.path.basename(abs_path)


def _resolve_indexed_fasta_path(file_path):
    if not file_path:
        return None
    path = os.path.normpath(str(file_path))
    base_dir = os.path.normpath(str(settings.BASE_DIR))

    if os.path.isabs(path):
        if os.path.isfile(path):
            return path
        # if path contains duplicate base_dir prefix, strip one copy
        if path.startswith(base_dir + os.path.sep):
            stripped = path[len(base_dir + os.path.sep):]
            candidate = os.path.join(base_dir, stripped)
            if os.path.isfile(candidate):
                return candidate
        return None

    base_name = os.path.basename(base_dir)
    if path.startswith(base_name + os.path.sep):
        stripped = path[len(base_name + os.path.sep):]
        candidate = os.path.join(base_dir, stripped)
        if os.path.isfile(candidate):
            return candidate
        path = stripped

    static_index = path.find(os.path.join('static', 'data'))
    if static_index != -1:
        candidate = os.path.normpath(path[static_index:])
        abs_candidate = os.path.join(base_dir, candidate)
        if os.path.isfile(abs_candidate):
            return abs_candidate

    candidate = os.path.join(base_dir, path)
    if os.path.isfile(candidate):
        return candidate

    # Fallback: search by basename under static/data
    basename = os.path.basename(path)
    if basename:
        data_root = os.path.join(base_dir, 'static', 'data')
        for root, _, files in os.walk(data_root):
            if basename in files:
                candidate = os.path.join(root, basename)
                if os.path.isfile(candidate):
                    return candidate
    return None


def _npm_global_root():
    npm_path = shutil.which('npm')
    if not npm_path:
        return None
    try:
        result = subprocess.run([npm_path, 'root', '-g'], capture_output=True, text=True, check=True)
        root = result.stdout.strip()
        if os.path.isdir(root):
            return root
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        pass
    return None


def _find_jbrowse_build_root():
    candidates = [
        os.path.join(settings.BASE_DIR, 'static', 'jbrowse'),
        os.path.join(settings.BASE_DIR, 'node_modules', '@jbrowse', 'web', 'dist'),
    ]
    npm_root = _npm_global_root()
    if npm_root:
        candidates.append(os.path.join(npm_root, '@jbrowse', 'web', 'dist'))
    for candidate in candidates:
        if os.path.isdir(candidate):
            return candidate
    return None


def _find_jbrowse_asset(asset_name):
    root = _find_jbrowse_build_root()
    if not root:
        return None
    asset_path = os.path.join(root, asset_name)
    return asset_path if os.path.exists(asset_path) else None


def jbrowse_asset(request, asset_name):
    allowed_assets = {
        'jbrowse-web.umd.production.min.js': 'application/javascript',
        'jbrowse-web.min.css': 'text/css',
    }
    if asset_name not in allowed_assets:
        return HttpResponseNotFound('Asset not found')
    asset_path = _find_jbrowse_asset(asset_name)
    if not asset_path:
        return HttpResponseNotFound('Asset not found')
    return FileResponse(open(asset_path, 'rb'), content_type=allowed_assets[asset_name])


def jbrowse_page(request):
    data = request.GET.get('data')
    track_param = request.GET.get('track')
    loc = request.GET.get('loc')
    data_dir = _resolve_jbrowse_data_path(data) if data else None
    fasta_url = None
    gff_url = None
    gff_index_url = None
    fai_url = None
    gzi_url = None
    fasta_adapter_type = None
    gff_adapter_type = None
    gff_file = None
    warnings = []
    species_list = []
    genome_table_data = []
    jbrowse_track = track_param if track_param else None
    jbrowse_has_epigenome = False

    if not data:
        species_qs = Taxonomy.objects.filter(genomes__isnull=False).distinct().order_by('species')
        for tax in species_qs:
            species_dir = _format_species(tax.species)
            if _has_jbrowse_data_directory(species_dir):
                species_list.append({
                    'species': tax.species,
                    'species_dir': species_dir,
                    'taxonomy_id': tax.ncbi_taxonomy_id,
                })
        genome_qs = Genome.objects.select_related('taxonomy').order_by('taxonomy__species')
        for row in genome_qs:
            if row.taxonomy:
                species_dir = _format_species(row.taxonomy.species)
                genome_table_data.append({
                    'taxonomy': row.taxonomy,
                    'assembly_level': row.assembly_level,
                    'size_mb': row.size_mb,
                    'gc_content': row.gc_content,
                    'scaffold_n50': row.scaffold_n50,
                    'completeness_genome_sd': row.completeness_genome,
                    'gene_number': row.gene_number,
                    'completeness_gene_sd': row.completeness_gene,
                    'gene_number_annevo': row.gene_number_annevo,
                    'completeness_gene_annevo_sd': row.completeness_gene_annevo,
                    'species_dir': species_dir,
                    'has_jbrowse_data': _has_jbrowse_data_directory(species_dir),
                })

    if data_dir and not loc:
        loc = _find_default_jbrowse_loc(data_dir)

    if data_dir:
        species_dir = os.path.relpath(data_dir, os.path.join(settings.BASE_DIR, 'static', 'data')).split(os.sep)[0]
        jbrowse_has_epigenome = _has_epigenome_data_directory(species_dir)

        fasta_file = _find_preferred_fasta_file(data_dir)
        gff_file = _find_preferred_annotation_file(data_dir)
        if fasta_file:
            fasta_url = _format_static_url(fasta_file)
            fasta_adapter_type = 'BgzipFastaAdapter' if fasta_file.endswith('.gz') else 'IndexedFastaAdapter'
            fasta_indexes = _find_fasta_indexes(fasta_file)
            if fasta_indexes.get('fai'):
                fai_url = _format_static_url(fasta_indexes.get('fai'))
            if fasta_indexes.get('gzi'):
                gzi_url = _format_static_url(fasta_indexes.get('gzi'))
            if not fasta_indexes.get('fai'):
                warnings.append('Genome FASTA index file (.fai) was not found. JBrowse may not be able to load the reference sequence until an index is generated.')
        else:
            warnings.append('No FASTA genome sequence was found in the requested data directory.')
        if gff_file:
            gff_url = _format_static_url(gff_file)
            if gff_file.lower().endswith('.gtf.gz'):
                gff_adapter_type = 'GtfAdapter'
            elif gff_file.lower().endswith('.gtf'):
                gff_adapter_type = 'GtfAdapter'
            elif gff_file.lower().endswith('.gff3.gz'):
                index_path = _find_annotation_index(gff_file)
                if index_path:
                    gff_adapter_type = 'Gff3TabixAdapter'
                    gff_index_url = _format_static_url(index_path)
                else:
                    gff_adapter_type = 'Gff3Adapter'
            else:
                gff_adapter_type = 'Gff3Adapter'
        else:
            warnings.append('No GFF or GTF annotation file was found in the requested data directory.')
    elif data is not None:
        warnings.append('The specified JBrowse data directory is invalid or missing.')

    asset_root = _find_jbrowse_build_root()
    jbrowse_js_url = None
    jbrowse_css_url = None
    jbrowse_app_url = None
    has_local_jbrowse = False
    if asset_root:
        js_path = os.path.join(asset_root, 'jbrowse-web.umd.production.min.js')
        css_path = os.path.join(asset_root, 'jbrowse-web.min.css')
        index_path = os.path.join(asset_root, 'index.html')
        if os.path.exists(js_path) and os.path.exists(css_path):
            jbrowse_js_url = f"/jbrowse-asset/jbrowse-web.umd.production.min.js"
            jbrowse_css_url = f"/jbrowse-asset/jbrowse-web.min.css"
            has_local_jbrowse = True
        elif os.path.exists(index_path):
            jbrowse_app_url = '/static/jbrowse/index.html'
            has_local_jbrowse = True
        else:
            warnings.append('JBrowse 2 build files were found but the required JS/CSS files are missing in the build directory.')
    else:
        warnings.append('JBrowse 2 web build directory was not found locally. Install @jbrowse/web or place JBrowse assets under static/jbrowse.')

    if not jbrowse_app_url and asset_root:
        index_path = os.path.join(asset_root, 'index.html')
        if os.path.exists(index_path):
            jbrowse_app_url = '/static/jbrowse/index.html'
            has_local_jbrowse = True

    if not jbrowse_js_url and not jbrowse_app_url:
        jbrowse_js_url = '/static/jbrowse/jbrowse-web.umd.production.min.js'
    if not jbrowse_css_url and not jbrowse_app_url:
        jbrowse_css_url = '/static/jbrowse/jbrowse-web.min.css'

    jbrowse_plugin_url = None
    plugin_path = os.path.join(settings.BASE_DIR, 'static', 'jbrowse', 'umd_plugin.js')
    if os.path.exists(plugin_path):
        jbrowse_plugin_url = '/static/jbrowse/umd_plugin.js'

    jbrowse_tax = _find_taxonomy_for_jbrowse_data(data_dir) if data_dir else None
    jbrowse_track_name = _get_annotation_track_name(gff_file) if gff_file else 'Genome annotation'

    return render(request, 'jbrowse.html', {
        'page_title': 'JBrowse 2 Viewer',
        'jbrowse_js_url': jbrowse_js_url,
        'jbrowse_css_url': jbrowse_css_url,
        'jbrowse_plugin_url': jbrowse_plugin_url,
        'jbrowse_app_url': jbrowse_app_url,
        'has_local_jbrowse': has_local_jbrowse,
        'jbrowse_data': data,
        'jbrowse_track': jbrowse_track,
        'jbrowse_loc': loc,
        'fasta_url': fasta_url,
        'gff_url': gff_url,
        'fai_url': fai_url,
        'gzi_url': gzi_url,
        'fasta_adapter_type': fasta_adapter_type,
        'gff_adapter_type': gff_adapter_type,
        'gff_index_url': gff_index_url,
        'warnings': warnings,
        'species_list': species_list,
        'genome_table_data': genome_table_data,
        'jbrowse_tax': jbrowse_tax,
        'jbrowse_has_epigenome': jbrowse_has_epigenome,
        'jbrowse_track_name': jbrowse_track_name,
    })


def jbrowse_config(request):
    data = request.GET.get('data')
    track_param = request.GET.get('track')
    loc = request.GET.get('loc')
    data_dir = _resolve_jbrowse_data_path(data)
    if not data_dir:
        return JsonResponse({
            'error': 'Invalid or missing JBrowse data directory. Ensure the requested dataset exists under static/data/<species>/genome.',
        }, status=400)

    fasta_file = _find_preferred_fasta_file(data_dir)
    annotation_files = _find_annotation_files(data_dir)
    if track_param == 'epigenome':
        if not fasta_file:
            return JsonResponse({
                'error': 'Genome FASTA file is missing in the selected data directory. Epigenome visualization requires a reference FASTA.',
            }, status=400)
    elif not (fasta_file and annotation_files):
        return JsonResponse({
            'error': 'Genome FASTA and/or annotation files are missing in the selected data directory. The dataset directory must contain both reference sequence and annotation files.',
        }, status=400)

    fasta_url = _format_static_url(fasta_file)
    fasta_adapter_type = 'BgzipFastaAdapter' if fasta_file.endswith('.gz') else 'IndexedFastaAdapter'
    track_param = request.GET.get('track')

    fasta_adapter = {
        'type': fasta_adapter_type,
        'fastaLocation': {'uri': fasta_url, 'locationType': 'UriLocation'},
    }
    fasta_indexes = _find_fasta_indexes(fasta_file)
    if fasta_adapter_type == 'IndexedFastaAdapter':
        if fasta_indexes.get('fai'):
            fasta_adapter['faiLocation'] = {'uri': _format_static_url(fasta_indexes.get('fai')), 'locationType': 'UriLocation'}
    else:
        if fasta_indexes.get('fai'):
            fasta_adapter['faiLocation'] = {'uri': _format_static_url(fasta_indexes.get('fai')), 'locationType': 'UriLocation'}
        if fasta_indexes.get('gzi'):
            fasta_adapter['gziLocation'] = {'uri': _format_static_url(fasta_indexes.get('gzi')), 'locationType': 'UriLocation'}

    def build_annotation_track(annotation_path, track_id):
        track_name = _get_annotation_track_name(annotation_path)
        uri = _format_static_url(annotation_path)
        if annotation_path.lower().endswith('.gtf.gz'):
            return {
                'type': 'FeatureTrack',
                'trackId': track_id,
                'name': track_name,
                'assemblyNames': ['assembly1'],
                'adapter': {
                    'type': 'GtfAdapter',
                    'gtfGzLocation': {'uri': uri, 'locationType': 'UriLocation'},
                },
            }
        if annotation_path.lower().endswith('.gtf'):
            return {
                'type': 'FeatureTrack',
                'trackId': track_id,
                'name': track_name,
                'assemblyNames': ['assembly1'],
                'adapter': {
                    'type': 'GtfAdapter',
                    'gtfLocation': {'uri': uri, 'locationType': 'UriLocation'},
                },
            }
        if annotation_path.lower().endswith('.gff3.gz'):
            index_path = _find_annotation_index(annotation_path)
            if index_path and _is_bgzip_file(annotation_path):
                return {
                    'type': 'FeatureTrack',
                    'trackId': track_id,
                    'name': track_name,
                    'assemblyNames': ['assembly1'],
                    'adapter': {
                        'type': 'Gff3TabixAdapter',
                        'gffGzLocation': {'uri': uri, 'locationType': 'UriLocation'},
                        'index': {'location': {'uri': _format_static_url(index_path), 'locationType': 'UriLocation'}},
                    },
                }
            return {
                'type': 'FeatureTrack',
                'trackId': track_id,
                'name': track_name,
                'assemblyNames': ['assembly1'],
                'adapter': {
                    'type': 'Gff3Adapter',
                    'gffGzLocation': {'uri': uri, 'locationType': 'UriLocation'},
                },
            }
        return {
            'type': 'FeatureTrack',
            'trackId': track_id,
            'name': track_name,
            'assemblyNames': ['assembly1'],
            'adapter': {
                'type': 'Gff3Adapter',
                'gffLocation': {'uri': uri, 'locationType': 'UriLocation'},
            },
        }

    tracks = [build_annotation_track(path, f'annotationTrack{i+1}') for i, path in enumerate(annotation_files)]

    def _find_subdirectory(dir_path, subfolder):
        candidate = os.path.join(os.path.dirname(dir_path), subfolder)
        if os.path.isdir(candidate):
            return candidate
        return None

    def _build_epigenome_track(path, track_id):
        return {
            'type': 'QuantitativeTrack',
            'trackId': track_id,
            'name': os.path.basename(path),
            'assemblyNames': ['assembly1'],
            'adapter': {
                'type': 'BigWigAdapter',
                'bigWigLocation': {'uri': _format_static_url(path), 'locationType': 'UriLocation'},
            },
        }

    if track_param == 'epigenome':
        epigenome_dir = None
        if os.path.basename(data_dir) == 'genome':
            epigenome_dir = os.path.join(os.path.dirname(data_dir), 'epigenome')
        else:
            epigenome_dir = os.path.join(data_dir, 'epigenome')
        if not os.path.isdir(epigenome_dir):
            epigenome_dir = os.path.join(os.path.dirname(data_dir), 'epigenome')
        if epigenome_dir and os.path.isdir(epigenome_dir):
            epigenome_files = _find_epigenome_files(epigenome_dir)
            if epigenome_files:
                tracks.extend([_build_epigenome_track(path, f'epigenomeTrack{i+1}') for i, path in enumerate(epigenome_files)])
            else:
                return JsonResponse({
                    'error': 'No epigenome BigWig files were found in the epigenome data directory.',
                }, status=400)
        else:
            return JsonResponse({
                'error': 'No epigenome data directory was found for the selected genome dataset.',
            }, status=400)

    config = {
        'assemblies': [
            {
                'name': 'assembly1',
                'displayName': 'Genome assembly',
                'sequence': {
                    'type': 'ReferenceSequenceTrack',
                    'trackId': 'referenceSequenceTrack',
                    'adapter': fasta_adapter,
                },
            }
        ],
        'tracks': tracks,
    }

    if loc:
        loc_match = re.match(r'^([^:]+):(\d+)\.\.(\d+)$', loc)
        if loc_match:
            ref_name, start, end = loc_match.groups()
            config['defaultSession'] = {
                'name': 'Default session',
                'views': [
                    {
                        'id': 'linearGenomeView',
                        'type': 'LinearGenomeView',
                        'assemblyName': 'assembly1',
                        'displayedRegions': [
                            {
                                'assemblyName': 'assembly1',
                                'refName': ref_name,
                                'start': int(start),
                                'end': int(end),
                            }
                        ],
                    }
                ],
            }
    return JsonResponse(config)


@csrf_exempt
def run_blast(request):
    if request.method != 'POST':
        return JsonResponse({'status': 'error', 'message': 'POST only'}, status=405)

    try:
        data = json.loads(request.body)
    except Exception:
        return JsonResponse({'status': 'error', 'message': 'Invalid JSON body.'}, status=400)

    sequence = str(data.get('sequence', '')).strip()
    if not sequence:
        return JsonResponse({'status': 'error', 'message': 'Please provide a query sequence or upload a FASTA file.'}, status=400)

    if str(data.get('sync', '')).lower() in ['true', '1', 'yes']:
        try:
            results = _perform_blast(data)
            return JsonResponse({'status': 'success', 'results': results})
        except subprocess.CalledProcessError as exc:
            message = 'BLAST failed: ' + (exc.stderr.strip() if getattr(exc, 'stderr', None) else str(exc))
            return JsonResponse({'status': 'error', 'message': message})
        except subprocess.TimeoutExpired:
            return JsonResponse({'status': 'error', 'message': f'BLAST timed out after {BLAST_CMD_TIMEOUT_SECONDS} seconds. Try a smaller query or review the database size.'})
        except Exception as exc:
            return JsonResponse({'status': 'error', 'message': str(exc)})

    job_id = uuid.uuid4().hex
    job_url = request.build_absolute_uri(reverse('blast_page') + f'?job_id={job_id}')
    email_address = str(data.get('email', '') or '').strip()
    if email_address:
        data['email_address'] = email_address
    data['job_url'] = job_url

    _write_blast_job_status(job_id, {
        'status': 'queued',
        'created_at': datetime.utcnow().isoformat() + 'Z',
        'program': data.get('program'),
        'database': data.get('database'),
        'message': 'BLAST job queued. Poll status for updates.',
        'email_address': email_address,
        'job_url': job_url,
    })

    _enqueue_job(_blast_job_worker, job_id, data)

    if email_address:
        subject = 'Your AnnelidaDB BLAST job is queued'
        message = (
            f'Your BLAST job has been queued on AnnelidaDB.\n\n'
            f'You can track its status using the following link:\n{job_url}\n\n'
            'Please keep this link to revisit your job status.'
        )
        _send_notification_email(subject, message, email_address)

    return JsonResponse({'status': 'pending', 'job_id': job_id, 'job_url': job_url})


@require_GET
def blast_status(request, job_id):
    job = _read_blast_job_status(job_id)
    if job is None:
        return JsonResponse({'status': 'error', 'message': 'BLAST job not found.'}, status=404)
    job['job_id'] = job_id
    job['job_url'] = request.build_absolute_uri(reverse('blast_page') + f'?job_id={job_id}')

    if job.get('status') in {'success', 'error'} and job.get('email_address') and not job.get('completion_email_sent'):
        subject = 'Your AnnelidaDB BLAST job completed' if job['status'] == 'success' else 'Your AnnelidaDB BLAST job failed'
        body = [
            f'Your BLAST job has {"completed successfully" if job["status"] == "success" else "finished with an error"}.',
            f'Status: {job["status"]}',
        ]
        if job['job_url']:
            body.append(f'View job status: {job["job_url"]}')

        _send_notification_email(subject, '\n\n'.join(body), job['email_address'])
        job['completion_email_sent'] = True
        _write_blast_job_status(job_id, job)

    return JsonResponse(job)


@require_GET
def gene_annotation(request):
    gene_id = str(request.GET.get('gene_id', '')).strip()
    taxonomy_id = str(request.GET.get('taxonomy_id', '')).strip()
    species = str(request.GET.get('species', '')).strip()

    if not gene_id:
        return JsonResponse({'error': 'gene_id is required'}, status=400)

    formatted_species, tax = _find_annotation_species_directory(taxonomy_id=taxonomy_id, species=species)
    if not formatted_species:
        return JsonResponse({'error': 'taxonomy_id or species is required'}, status=400)

    annotation_data = {
        'gene_id': gene_id,
        'taxonomy_id': taxonomy_id or (tax.ncbi_taxonomy_id if tax else None),
        'species': species or (tax.species if tax else None),
        'eggnog': {
            'header': [],
            'rows': [],
            'source': None,
        },
        'interpro': {
            'header': [],
            'rows': [],
            'member_dbs': [],
            'source': None,
        },
    }

    eggnog_header, eggnog_rows, eggnog_source = _find_eggnog_annotation_rows_for_gene_dir(formatted_species, gene_id)
    annotation_data['eggnog'] = {
        'header': eggnog_header,
        'rows': eggnog_rows,
        'source': eggnog_source,
    }

    interpro_header, interpro_rows, member_dbs, interpro_source = _find_interproscan_rows_for_gene_dir(formatted_species, gene_id)
    annotation_data['interpro'] = {
        'header': interpro_header,
        'rows': interpro_rows,
        'member_dbs': member_dbs,
        'source': interpro_source,
    }

    return JsonResponse(annotation_data)


@require_GET
def gene_detail_extras(request):
    gene_id = str(request.GET.get('gene_id', '')).strip()
    taxonomy_id = str(request.GET.get('taxonomy_id', '')).strip()
    species = str(request.GET.get('species', '')).strip()

    if not gene_id:
        return JsonResponse({'error': 'gene_id is required'}, status=400)

    formatted_species, tax = _find_annotation_species_directory(taxonomy_id=taxonomy_id, species=species)
    if not formatted_species:
        return JsonResponse({'error': 'taxonomy_id or species is required'}, status=400)

    annotation_fields = []
    annotation_description = 'No description available.'
    if tax:
        try:
            annotation_obj = EggnogAnnotation.objects.get(query=gene_id, taxonomy=tax)
            annotation_fields = _collect_annotation_fields(annotation_obj)
            annotation_description = annotation_obj.description or annotation_description
        except EggnogAnnotation.DoesNotExist:
            pass

    if not annotation_fields:
        header_cells, matching_rows, eggnog_source = _find_eggnog_annotation_rows_for_gene_dir(formatted_species, gene_id)
        if matching_rows:
            row = matching_rows[0]
            annotation_fields = _build_annotation_fields_from_file(header_cells, row)
            header_index = {label: idx for idx, label in enumerate(header_cells)}
            description = row[header_index['Description']] if 'Description' in header_index else None
            annotation_description = description or annotation_description

    annotation_fields = [(label, value) for label, value in annotation_fields if label != 'Preferred name']

    interpro_header, interpro_rows, member_dbs, interpro_source = _find_interproscan_rows_for_gene_dir(formatted_species, gene_id)
    gene_family_memberships = _find_gene_family_memberships_by_source(gene_id)

    return JsonResponse({
        'gene_id': gene_id,
        'taxonomy_id': taxonomy_id or (tax.ncbi_taxonomy_id if tax else None),
        'species': species or (tax.species if tax else None),
        'annotation_description': annotation_description,
        'annotation_fields': annotation_fields,
        'interpro': {
            'header': interpro_header,
            'rows': interpro_rows,
            'member_dbs': member_dbs,
            'source': interpro_source,
        },
        'gene_family_memberships': gene_family_memberships,
    })


@require_GET
def blast_job_file(request, job_id, filename):
    safe_name = os.path.basename(filename)
    if safe_name != filename:
        return HttpResponseNotFound('Invalid file requested.')
    job_dir = _blast_job_dir(job_id)
    file_path = os.path.join(job_dir, safe_name)
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return HttpResponseNotFound('Requested BLAST file not found.')
    return FileResponse(open(file_path, 'rb'), as_attachment=True, filename=safe_name)


# ──────────────────────────────────────────────
# Statistics
# ──────────────────────────────────────────────

@require_POST
def refresh_statistics(request):
    if not request.user.is_staff:
        return HttpResponseForbidden('Staff authentication required to refresh statistics.')

    stats = _get_statistics_data(force_refresh=True)
    return JsonResponse({
        'status': 'success',
        'message': 'Statistics refreshed successfully.',
        'data': stats,
    })


def get_statistics(request):
    return JsonResponse({
        'status': 'success',
        'data': _get_statistics_data()
    })


# ──────────────────────────────────────────────
# Table pages
# ──────────────────────────────────────────────

@cache_page(300)
def taxonomy_page(request):
    stats = _get_statistics_data()
    taxonomies = Taxonomy.objects.annotate(
        genome_count=Count('genomes', distinct=True),
        mitogenome_count=Count('mitogenomes', distinct=True),
        singlecell_dataset_count=Count('singlecell_datasets', distinct=True),
    ).all()
    singlecell_counts = {
        str(row['species']).strip().casefold(): row['total']
        for row in SingleCell.objects.values('species').annotate(total=Count('id'))
        if row['species']
    }
    symbiont_counts = {
        str(row['host_species']).strip().casefold(): row['total']
        for row in Symbiont.objects.values('host_species').annotate(total=Count('id'))
        if row['host_species']
    }
    epigenome_records = _load_epigenome_records()
    epigenome_counts = Counter(
        _normalize_species_name(row.get('species'))
        for row in epigenome_records
        if row.get('species')
    )
    proteome_records = _load_proteome_records()
    proteome_counts = Counter(
        _normalize_species_name(row.get('species'))
        for row in proteome_records
        if row.get('species')
    )
    ncbi_cache = _get_ncbi_omics_cache()
    counts = ncbi_cache.get('counts', {})
    cached_dna_by_taxid, cached_dna_by_species = _cached_dna_counts()
    cached_rna_counts = _cached_rna_counts_by_species()
    dna_species_count = _get_dna_species_count(ncbi_cache)
    enhanced_data = []
    for tax in taxonomies:
        ncbi_id = tax.ncbi_taxonomy_id
        entry = counts.get(str(ncbi_id), {
            'dna_runs': [],
            'rna_runs': [],
            'dna_count': 0,
            'rna_count': 0,
        })
        dna_count = int(entry.get('dna_count', 0) or 0)
        rna_count = int(entry.get('rna_count', 0) or 0)
        if dna_count == 0:
            dna_count = cached_dna_by_taxid.get(str(ncbi_id).strip(), 0)
            if dna_count == 0:
                dna_count = cached_dna_by_species.get(_normalize_species_name(tax.species), 0)
        if rna_count == 0:
            rna_count = cached_rna_counts.get(_normalize_species_name(tax.species), 0)
        has_dna = dna_count > 0 or bool(entry.get('dna_runs'))
        has_rna = rna_count > 0 or bool(entry.get('rna_runs'))
        has_genome = tax.genome_count > 0
        has_mitogenome = tax.mitogenome_count > 0
        species_name_key = _normalize_species_name(tax.species)
        proteome_csv_count = proteome_counts.get(species_name_key, 0)
        proteome_study_count = proteome_csv_count + (1 if getattr(tax, 'proteome_link', None) else 0)
        has_proteome_csv = proteome_csv_count > 0
        has_proteome = bool(getattr(tax, 'proteome_link', None) or has_proteome_csv)
        epigenome_study_count = epigenome_counts.get(species_name_key, 0)
        has_epigenome_csv = epigenome_study_count > 0
        has_epigenome = bool(getattr(tax, 'epigenome_link', None) or has_epigenome_csv)
        species_key = str(tax.species or '').strip().casefold()
        singlecell_count = singlecell_counts.get(species_key, tax.singlecell_dataset_count)
        symbiont_count = symbiont_counts.get(species_key, 0)
        enhanced_data.append({
            'ncbi_taxonomy_id': ncbi_id,
            'species': tax.species,
            'class_name': tax.class_name,
            'order': tax.order,
            'family': tax.family,
            'updated_clade1': tax.updated_clade1,
            'updated_clade2': tax.updated_clade2,
            'updated_clade3': tax.updated_clade3,
            'updated_clade4': tax.updated_clade4,
            'updated_clade5': tax.updated_clade5,
            'updated_clade6': tax.updated_clade6,
            'has_genome': has_genome,
            'has_mitogenome': has_mitogenome,
            'has_dna': has_dna,
            'has_rna': has_rna,
            'has_epigenome': has_epigenome,
            'has_epigenome_csv': has_epigenome_csv,
            'epigenome_study_count': epigenome_study_count,
            'has_proteome': has_proteome,
            'proteome_study_count': proteome_study_count,
            'singlecell_count': singlecell_count,
            'symbiont_count': symbiont_count,
            'dna_count': dna_count,
            'rna_count': rna_count,
            'dna_link': _ncbi_sra_search_url(ncbi_id, _ncbi_sra_dna_query(ncbi_id)) if has_dna else None,
            'rna_link': _ncbi_sra_search_url(ncbi_id, _ncbi_sra_rna_query(ncbi_id)) if has_rna else None,
            'epigenome_link': getattr(tax, 'epigenome_link', None),
            'proteome_link': getattr(tax, 'proteome_link', None),
        })

    taxonomy_summary = [
        {
            'title': 'Genome Species',
            'value': stats['genome_count'],
            'border_class': 'border-info',
            'subtitle': 'Species with at least one genome assembly',
        },
        {
            'title': 'SRA (RNA) Species',
            'value': stats['rna_species_count'],
            'border_class': 'border-success',
            'subtitle': 'Species with RNA/transcriptome records on NCBI',
        },
        {
            'title': 'Mito Species',
            'value': stats['mitogenome_count'],
            'border_class': 'border-warning',
            'subtitle': 'Species with mitochondrial assemblies',
        },
        {
            'title': 'SRA (DNA) Species',
            'value': _format_count_if_large(dna_species_count),
            'border_class': 'border-primary',
            'subtitle': 'Species with DNA sequence records on NCBI',
        },
        {
            'title': 'Epigenome Species',
            'value': stats['epigenome_species_count'],
            'border_class': 'border-secondary',
            'subtitle': 'Species with epigenome files or links',
        },
        {
            'title': 'Proteome Species',
            'value': stats['proteome_species_count'],
            'border_class': 'border-dark',
            'subtitle': 'Species with proteome files or links',
        },
        {
            'title': 'Single-cell Species',
            'value': stats['singlecell_species_count'],
            'border_class': 'border-info',
            'subtitle': 'Species with single-cell datasets',
        },
        {
            'title': 'Symbiont Species',
            'value': stats['symbiont_count'],
            'border_class': 'border-success',
            'subtitle': 'Species with symbiont datasets',
        },
    ]

    return render(request, 'taxonomy.html', {
        'page_title': 'Taxonomy',
        'data_type': 'taxonomy',
        'data': enhanced_data,
        'page_summary_cards': taxonomy_summary,
    })


def proteome_page(request):
    records = _load_proteome_records()
    publications = _load_proteome_publications()
    stats = _get_statistics_data()
    proteome_species_count = stats.get('proteome_species_count', 0)
    enhanced_data = []
    unique_species = set()
    unique_papers = set()

    for row in records:
        tax = None
        species_url = None
        if row['species']:
            tax = get_taxonomy_for_species_name(row['species'])
            if tax:
                species_url = reverse('species_detail', args=[tax.ncbi_taxonomy_id or tax.species])
            unique_species.add(row['species'])

        paper_id = (row.get('paper_id') or '').strip()
        if paper_id:
            unique_papers.add(paper_id)

        enhanced_data.append({
            'class_name': tax.class_name if tax else None,
            'order': row['order'],
            'family': row['family'],
            'species': row['species'],
            'species_url': species_url,
            'updated_clade6': tax.updated_clade6 if tax else None,
            'updated_clade5': tax.updated_clade5 if tax else None,
            'updated_clade4': tax.updated_clade4 if tax else None,
            'updated_clade3': tax.updated_clade3 if tax else None,
            'updated_clade2': tax.updated_clade2 if tax else None,
            'updated_clade1': tax.updated_clade1 if tax else None,
            'modality': row['modality'],
            'tissue': row['tissue'],
            'stage': row['stage'],
            'paper_id': paper_id,
        })

    proteome_summary = [
        {
            'title': 'Total Studies',
            'value': len(records),
            'border_class': 'border-warning',
            'subtitle': 'Rows in proteome.csv',
        },
        {
            'title': 'Unique Papers',
            'value': len(unique_papers),
            'border_class': 'border-primary',
            'subtitle': 'Distinct Paper_ID values',
        },
        {
            'title': 'Unique Species',
            'value': proteome_species_count,
            'border_class': 'border-info',
            'subtitle': 'Species with proteome files or links',
        },
        {
            'title': 'Connected Publications',
            'value': len(publications),
            'border_class': 'border-success',
            'subtitle': 'Publications in proteome_publication.csv',
        },
    ]

    return render(request, 'proteome.html', {
        'page_title': 'Proteome',
        'active_page': 'proteome',
        'data_type': 'proteome',
        'data': enhanced_data,
        'page_summary_cards': proteome_summary,
    })


def epigenome_page(request):
    records = _load_epigenome_records()
    publications = _load_epigenome_publications()
    enhanced_data = []
    unique_species = set()
    unique_papers = set()

    for row in records:
        tax = None
        species_url = None
        if row['species']:
            tax = get_taxonomy_for_species_name(row['species'])
            if tax:
                species_url = reverse('species_detail', args=[tax.ncbi_taxonomy_id or tax.species])
            unique_species.add(row['species'])

        paper_id = (row.get('paper_id') or '').strip()
        if paper_id:
            unique_papers.add(paper_id)

        enhanced_data.append({
            'class_name': tax.class_name if tax else None,
            'order': row['order'],
            'family': row['family'],
            'species': row['species'],
            'species_url': species_url,
            'updated_clade6': tax.updated_clade6 if tax else None,
            'updated_clade5': tax.updated_clade5 if tax else None,
            'updated_clade4': tax.updated_clade4 if tax else None,
            'updated_clade3': tax.updated_clade3 if tax else None,
            'updated_clade2': tax.updated_clade2 if tax else None,
            'updated_clade1': tax.updated_clade1 if tax else None,
            'modality': row['modality'],
            'tissue': row['tissue'],
            'stage': row['stage'],
            'paper_id': paper_id,
        })

    epigenome_summary = [
        {
            'title': 'Total Studies',
            'value': len(records),
            'border_class': 'border-warning',
            'subtitle': 'Rows in epigenome.csv',
        },
        {
            'title': 'Unique Papers',
            'value': len(unique_papers),
            'border_class': 'border-primary',
            'subtitle': 'Distinct Paper_ID values',
        },
        {
            'title': 'Unique Species',
            'value': len(unique_species),
            'border_class': 'border-info',
            'subtitle': 'Species covered by epigenome records',
        },
        {
            'title': 'Connected Publications',
            'value': len(publications),
            'border_class': 'border-success',
            'subtitle': 'Publications in epigenome_publication.csv',
        },
    ]

    return render(request, 'epigenome.html', {
        'page_title': 'Epigenome',
        'active_page': 'epigenome',
        'data_type': 'epigenome',
        'data': enhanced_data,
        'page_summary_cards': epigenome_summary,
    })


@lru_cache(maxsize=1)
def _load_epigenome_records():
    path = os.path.join(settings.BASE_DIR, 'epigenome.csv')
    records = []
    if not os.path.exists(path):
        return records
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append({
                'order': (row.get('Order') or '').strip(),
                'family': (row.get('Family') or '').strip(),
                'species': (row.get('Species') or '').strip(),
                'modality': (row.get('Modality') or '').strip(),
                'year': (row.get('Year') or '').strip(),
                'tissue': (row.get('Tissue') or '').strip(),
                'stage': (row.get('Stage') or '').strip(),
                'title': (row.get('Title') or '').strip(),
                'paper_id': (row.get('Paper_ID') or '').strip(),
            })
    return records


@lru_cache(maxsize=1)
def _load_epigenome_publications():
    path = os.path.join(settings.BASE_DIR, 'epigenome_publication.csv')
    pubs = {}
    if not os.path.exists(path):
        return pubs
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            paper_id = (row.get('Paper_ID') or '').strip()
            if not paper_id:
                continue
            pubs[paper_id] = {
                'paper_id': paper_id,
                'species_list': (row.get('Species_list') or '').strip(),
                'modality': (row.get('Modality') or '').strip(),
                'year': (row.get('Year') or '').strip(),
                'authors': (row.get('Authors') or '').strip(),
                'title': (row.get('Title') or '').strip(),
                'doi': (row.get('DOI') or '').strip(),
                'source_url': (row.get('Source_URL') or '').strip(),
                'abstract': (row.get('Abstract_Excerpt') or '').strip(),
                'accession_numbers': (row.get('Accession_numbers') or '').strip(),
            }
    return pubs


def _load_proteome_records():
    path = os.path.join(settings.BASE_DIR, 'proteome.csv')
    records = []
    if not os.path.exists(path):
        return records
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append({
                'order': (row.get('Order') or '').strip(),
                'family': (row.get('Family') or '').strip(),
                'species': (row.get('Species') or '').strip(),
                'modality': (row.get('Modality') or '').strip(),
                'year': (row.get('Year') or '').strip(),
                'tissue': (row.get('Tissue') or '').strip(),
                'stage': (row.get('Stage') or '').strip(),
                'title': (row.get('Title') or '').strip(),
                'paper_id': (row.get('Paper_ID') or '').strip(),
            })
    return records


@lru_cache(maxsize=1)
def _load_proteome_publications():
    path = os.path.join(settings.BASE_DIR, 'proteome_publication.csv')
    pubs = {}
    if not os.path.exists(path):
        return pubs
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            paper_id = (row.get('Paper_ID') or '').strip()
            if not paper_id:
                continue
            pubs[paper_id] = {
                'paper_id': paper_id,
                'species_list': (row.get('Species_list') or '').strip(),
                'modality': (row.get('Modality') or '').strip(),
                'year': (row.get('Year') or '').strip(),
                'authors': (row.get('Authors') or '').strip(),
                'title': (row.get('Title') or '').strip(),
                'doi': (row.get('DOI') or '').strip(),
                'source_url': (row.get('Source_URL') or '').strip(),
                'abstract': (row.get('Abstract_Excerpt') or '').strip(),
                'accession_numbers': (row.get('Accession_numbers') or '').strip(),
            }
    return pubs


def proteome_detail_page(request, paper_id):
    publications = _load_proteome_publications()
    publication = publications.get(paper_id)
    if not publication:
        return HttpResponseNotFound('Proteome publication not found.')
    records = [row for row in _load_proteome_records() if row['paper_id'] == paper_id]
    _attach_species_detail_links(records)
    publication_species_links = []
    for species in _split_publication_species(publication.get('species_list')):
        species_url = None
        tax = get_taxonomy_for_species_name(species)
        if tax:
            species_url = reverse('species_detail', args=[tax.ncbi_taxonomy_id or tax.species])
        publication_species_links.append({
            'name': species,
            'url': species_url,
        })
    return render(request, 'proteome_detail.html', {
        'page_title': f'Proteome Publication {paper_id}',
        'publication': publication,
        'records': records,
        'related_species': sorted(set(r['species'] for r in records if r['species'])),
        'publication_species_links': publication_species_links,
        'accession_links': _parse_accession_numbers(publication.get('accession_numbers')),
        'active_page': 'proteome',
    })


def _split_publication_species(species_list):
    if not species_list:
        return []
    return [item.strip() for item in re.split(r'\s*;\s*|\s*,\s*', species_list) if item.strip()]


def _accession_url(accession):
    accession = accession.strip()
    if not accession:
        return None
    if ':' in accession:
        prefix, value = accession.split(':', 1)
        prefix = prefix.strip().lower()
        value = value.strip()
        if prefix == 'github':
            return f'https://github.com/{value}'
        if prefix == 'zenodo':
            return f'https://doi.org/{value}'
        if prefix == 'doi':
            return f'https://doi.org/{value}'
    if re.match(r'^(SRR|ERR|DRR)\d+(?:-\d+)?$', accession, re.I):
        return f'https://www.ncbi.nlm.nih.gov/sra?term={urllib.parse.quote(accession)}'
    if re.match(r'^(PRJ[A-Z]*\d+)$', accession, re.I):
        return f'https://www.ncbi.nlm.nih.gov/bioproject/{urllib.parse.quote(accession)}'
    if re.match(r'^(ERP\d+|DRP\d+|SRP\d+)(?:-\d+)?$', accession, re.I):
        return f'https://www.ncbi.nlm.nih.gov/sra?term={urllib.parse.quote(accession)}'
    if re.match(r'^(GSE|GSM)\d+(?:-\d+)?$', accession, re.I):
        return f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={urllib.parse.quote(accession)}'
    if re.match(r'^(GCA|GCF)_\d+\.\d+$', accession, re.I):
        return f'https://www.ncbi.nlm.nih.gov/assembly/{urllib.parse.quote(accession)}'
    if re.match(r'^(ERS|ERR|DRS|SAME)\d+(?:-\d+)?$', accession, re.I):
        return f'https://www.ebi.ac.uk/ena/browser/view/{urllib.parse.quote(accession)}'
    if re.match(r'^PXD\d+', accession, re.I):
        return f'https://www.ebi.ac.uk/pride/archive/projects/{urllib.parse.quote(accession)}'
    return f'https://www.ncbi.nlm.nih.gov/search/all/?term={urllib.parse.quote(accession)}'


def _parse_accession_numbers(accession_numbers):
    accessions = []
    if not accession_numbers:
        return accessions
    for token in re.split(r'\s*;\s*|\s*,\s*', accession_numbers):
        token = token.strip()
        if not token or token.lower() == 'not reported':
            continue
        accessions.append({
            'label': token,
            'url': _accession_url(token),
        })
    return accessions


def epigenome_detail_page(request, paper_id):
    publications = _load_epigenome_publications()
    publication = publications.get(paper_id)
    if not publication:
        return HttpResponseNotFound('Epigenome publication not found.')
    records = [row for row in _load_epigenome_records() if row['paper_id'] == paper_id]
    _attach_species_detail_links(records)
    publication_species_links = []
    for species in _split_publication_species(publication.get('species_list')):
        species_url = None
        tax = get_taxonomy_for_species_name(species)
        if tax:
            species_url = reverse('species_detail', args=[tax.ncbi_taxonomy_id or tax.species])
        publication_species_links.append({
            'name': species,
            'url': species_url,
        })
    return render(request, 'epigenome_detail.html', {
        'page_title': f'Epigenome Publication {paper_id}',
        'publication': publication,
        'records': records,
        'related_species': sorted(set(r['species'] for r in records if r['species'])),
        'publication_species_links': publication_species_links,
        'accession_links': _parse_accession_numbers(publication.get('accession_numbers')),
        'active_page': 'epigenome',
    })


def _matches_species_value(source_species, target_species):
    if not source_species or not target_species:
        return False
    source_norm = _normalize_species_name(source_species)
    target_norm = _normalize_species_name(target_species)
    if not source_norm or not target_norm:
        return False
    if source_norm == target_norm:
        return True
    return source_norm in target_norm or target_norm in source_norm


def _get_epigenome_records_for_species(species):
    if not species:
        return []
    return [row for row in _load_epigenome_records() if _matches_species_value(row.get('species'), species)]


def _get_proteome_records_for_species(species):
    if not species:
        return []
    return [row for row in _load_proteome_records() if _matches_species_value(row.get('species'), species)]


def _has_epigenome_csv_data(species):
    return bool(_get_epigenome_records_for_species(species))


def _has_proteome_csv_data(species):
    return bool(_get_proteome_records_for_species(species))


_LOPHO_COMPLETENESS_MAP = None


def _get_lopho_completeness_map():
    global _LOPHO_COMPLETENESS_MAP
    if _LOPHO_COMPLETENESS_MAP is not None:
        return _LOPHO_COMPLETENESS_MAP
    path = os.path.join(settings.BASE_DIR, 'Genome_completeness.lophotrochozoa.csv')
    mapping = {}
    if os.path.exists(path):
        try:
            with open(path, 'r', encoding='utf-8-sig') as f:
                lines = [line.strip() for line in f if line.strip()]
            rows = list(csv.reader(lines))
            header_idx = -1
            for idx, r in enumerate(rows):
                if r and any('Species' in str(c) for c in r):
                    header_idx = idx
                    break
            if header_idx >= 0:
                for r in rows[header_idx + 1:]:
                    if r and r[0].strip():
                        sp_key = _normalize_species_name(' '.join(r[0].strip().split()))
                        g_raw = r[1].strip() if len(r) > 1 else ''
                        ge_raw = r[2].strip() if len(r) > 2 else ''
                        an_raw = r[3].strip() if len(r) > 3 else ''
                        mapping[sp_key] = {
                            'genome': None if g_raw.upper() in ('NA', '') else g_raw,
                            'gene': None if ge_raw.upper() in ('NA', '') else ge_raw,
                            'annevo': None if an_raw.upper() in ('NA', '') else an_raw,
                        }
        except Exception:
            pass
    _LOPHO_COMPLETENESS_MAP = mapping
    return mapping


def _enrich_genome_with_lopho_completeness(genome):
    if not genome or not getattr(genome, 'taxonomy', None):
        return genome
    if not getattr(genome, 'completeness_genome_lopho_raw', None):
        mapping = _get_lopho_completeness_map()
        sp_key = _normalize_species_name(genome.taxonomy.species)
        entry = mapping.get(sp_key)
        if entry:
            genome.completeness_genome_lopho_raw = entry['genome']
            genome.completeness_gene_lopho = entry['gene']
            genome.completeness_gene_annevo_lopho = entry['annevo']
    return genome


def _enrich_genomes_with_lopho_completeness(genomes):
    mapping = _get_lopho_completeness_map()
    for g in genomes:
        if getattr(g, 'taxonomy', None) and not getattr(g, 'completeness_genome_lopho_raw', None):
            sp_key = _normalize_species_name(g.taxonomy.species)
            entry = mapping.get(sp_key)
            if entry:
                g.completeness_genome_lopho_raw = entry['genome']
                g.completeness_gene_lopho = entry['gene']
                g.completeness_gene_annevo_lopho = entry['annevo']
    return genomes


def genome_page(request):
    data = Genome.objects.select_related('taxonomy').all()
    _enrich_genomes_with_lopho_completeness(data)
    stats = _get_statistics_data()
    total_genomes = stats['genome_count']
    total_species_genome = stats['genome_count']
    chromosome_genomes = stats['chromosome_genome_species_count']
    species_with_gene_number = stats['gene_model_species_count']

    genome_summary = [
        {
            'title': 'Total Genomes',
            'value': total_genomes,
            'border_class': 'border-warning',
        },
        {
            'title': 'Total Species (Genome)',
            'value': total_species_genome,
            'border_class': 'border-primary',
        },
        {
            'title': 'Chromosome-level Species',
            'value': chromosome_genomes,
            'border_class': 'border-info',
        },
        {
            'title': 'Species with gene model',
            'value': species_with_gene_number,
            'border_class': 'border-success',
        },
    ]

    return render(request, 'genome.html', {
        'page_title': 'Genome',
        'data_type': 'genome',
        'data': data,
        'page_summary_cards': genome_summary,
    })


def mitogenome_page(request):
    data = Mitogenome.objects.select_related('taxonomy').all()
    total_mitogenomes = data.count()
    total_species = _mitogenome_species_count()
    circular_count = data.filter(topology__icontains='Circular').count()
    length_values = [item.length_bp for item in data if item.length_bp is not None]
    average_length = int(sum(length_values) / len(length_values)) if length_values else 0

    mitogenome_summary = [
        {
            'title': 'Total Mito genomes',
            'value': total_mitogenomes,
            'border_class': 'border-info',
            'subtitle': 'Total mitogenome records',
        },
        {
            'title': 'Total Species',
            'value': total_species,
            'border_class': 'border-success',
            'subtitle': 'Distinct species represented',
        },
        {
            'title': 'Circular Topology',
            'value': circular_count,
            'border_class': 'border-warning',
            'subtitle': 'Assemblies marked as circular',
        },
        {
            'title': 'Avg. Length (bp)',
            'value': f"{average_length:,}",
            'border_class': 'border-primary',
            'subtitle': 'Average mitogenome length',
        },
    ]

    return render(request, 'mitogenome.html', {
        'page_title': 'Mitogenome',
        'data_type': 'mitogenome',
        'data': data,
        'page_summary_cards': mitogenome_summary,
    })


def symbiont_page(request):
    symbionts = Symbiont.objects.all()
    enhanced_data = []
    for sym in symbionts:
        host_tax = get_taxonomy_for_species_name(sym.host_species)
        host_species_url = None
        host_taxonomy_id = None
        if host_tax:
            host_taxonomy_id = host_tax.ncbi_taxonomy_id or host_tax.species
            host_species_url = reverse('species_detail', args=[host_taxonomy_id])
        enhanced_data.append({
            'host_species': sym.host_species,
            'sub': sym.sub,
            'host_taxonomy_id': host_taxonomy_id,
            'host_species_url': host_species_url,
            'host_class_name': host_tax.class_name if host_tax else 'Unknown',
            'host_order': host_tax.order if host_tax else 'Unknown',
            'host_family': host_tax.family if host_tax else 'Unknown',
            'symbiont_type': sym.symbiont_type,
            'size_bp': sym.size_bp,
            'gc_content': sym.gc_content,
            'contigs': sym.contigs,
            'completeness': sym.completeness,
            'contamination': sym.contamination,
            'domain': sym.domain,
            'phylum': sym.phylum,
            'class_name': sym.class_name,
            'order': sym.order,
            'family': sym.family,
            'genus': sym.genus,
            'species': sym.species,
        })
    host_species_count = symbionts.exclude(host_species__isnull=True).exclude(host_species='').values('host_species').distinct().count()
    symbiont_species_count = symbionts.count()
    complete_symbionts = symbionts.filter(contigs=1).count()
    total_genes = sum(item.genes or 0 for item in symbionts)

    symbiont_summary = [
        {
            'title': 'Host Species',
            'value': host_species_count,
            'border_class': 'border-info',
            'subtitle': 'Distinct host species represented',
        },
        {
            'title': 'Symbiont Species',
            'value': symbiont_species_count,
            'border_class': 'border-success',
            'subtitle': 'Symbiont genome records represented',
        },
        {
            'title': 'Complete Symbionts',
            'value': complete_symbionts,
            'border_class': 'border-warning',
            'subtitle': 'Symbionts assembled in a single contig',
        },
        {
            'title': 'Gene Number',
            'value': f"{total_genes:,}",
            'border_class': 'border-primary',
            'subtitle': 'Total annotated genes across symbionts',
        },
    ]

    return render(request, 'symbiont.html', {
        'page_title': 'Symbiont',
        'data_type': 'symbiont',
        'data': enhanced_data,
        'page_summary_cards': symbiont_summary,
    })


TRANSCRIPTOME_CACHE_DIR = os.path.join(settings.BASE_DIR, 'static', 'data', 'transcriptome')
TRANSCRIPTOME_CACHE_FILE = os.path.join(TRANSCRIPTOME_CACHE_DIR, 'sra_6340_runinfo.csv')
TRANSCRIPTOME_CACHE_META = os.path.join(TRANSCRIPTOME_CACHE_DIR, 'sra_6340_runinfo_meta.json')
DNA_SRA_CACHE_DIR = os.path.join(settings.BASE_DIR, 'static', 'data', 'dna_sra')
DNA_SRA_CACHE_FILE = os.path.join(DNA_SRA_CACHE_DIR, 'sra_6340_runinfo.csv')
DNA_SRA_CACHE_META = os.path.join(DNA_SRA_CACHE_DIR, 'sra_6340_runinfo_meta.json')

NCBI_OMICS_CACHE_DIR = os.path.join(settings.BASE_DIR, 'static', 'data', 'ncbi_omics')
NCBI_OMICS_CACHE_FILE = os.path.join(NCBI_OMICS_CACHE_DIR, 'ncbi_omics_counts.json')
NCBI_OMICS_CACHE_META = os.path.join(NCBI_OMICS_CACHE_DIR, 'ncbi_omics_counts_meta.json')


def _get_sample_transcriptome_rows():
    return [
        {
            'run': 'SRR000001',
            'study': 'SRP000001',
            'species': 'Capitella teleta',
            'tissue': 'Nervous system',
            'development_stage': 'Larva',
            'library_layout': 'PAIRED',
            'platform': 'Illumina HiSeq 2500',
            'spots_raw': 152000000,
            'bases_raw': 45700000000,
        },
        {
            'run': 'SRR000002',
            'study': 'SRP000002',
            'species': 'Pygospio elegans',
            'tissue': 'Whole organism',
            'development_stage': 'Juvenile',
            'library_layout': 'SINGLE',
            'platform': 'Illumina NovaSeq 6000',
            'spots_raw': 98000000,
            'bases_raw': 31200000000,
        },
        {
            'run': 'SRR000003',
            'study': 'SRP000003',
            'species': 'Platynereis dumerilii',
            'tissue': 'Larval stage',
            'development_stage': 'Embryo',
            'library_layout': 'PAIRED',
            'platform': 'Illumina HiSeq 4000',
            'spots_raw': 205000000,
            'bases_raw': 62800000000,
        },
    ]


def _save_cached_dna_sra_rows(rows):
    os.makedirs(DNA_SRA_CACHE_DIR, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    temporary_file = DNA_SRA_CACHE_FILE + '.tmp'
    with open(temporary_file, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary_file, DNA_SRA_CACHE_FILE)
    fetched_at = datetime.utcnow().isoformat() + 'Z'
    with open(DNA_SRA_CACHE_META, 'w', encoding='utf-8') as handle:
        json.dump({'fetched_at': fetched_at, 'row_count': len(rows)}, handle, indent=2)
    return fetched_at


def _load_cached_transcriptome_rows():
    if not os.path.isfile(TRANSCRIPTOME_CACHE_FILE):
        return None, None

    rows = _load_cached_transcriptome_rows_from_disk()
    if rows is not None:
        last_updated = None
        if os.path.isfile(TRANSCRIPTOME_CACHE_META):
            try:
                with open(TRANSCRIPTOME_CACHE_META, 'r', encoding='utf-8') as meta_handle:
                    meta = json.load(meta_handle)
                    last_updated = meta.get('fetched_at')
            except Exception:
                last_updated = None
        return rows, last_updated

    rows = []
    try:
        with open(TRANSCRIPTOME_CACHE_FILE, 'r', encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                raw_spots = _parse_integer_string(row.get('spots'))
                raw_bases = _parse_integer_string(row.get('bases'))
                accession_val, run_val = _get_transcriptome_accession_and_run(row)
                biosample_tissue = (row.get('BioSampleTissue') or '').strip()
                runinfo_tissue = row.get('SourceName', '').strip() or row.get('SampleName', '').strip() or row.get('Tissue', '').strip() or row.get('LibraryName', '').strip()
                rows.append({
                    'run': accession_val,
                    'study': run_val,
                    'bioproject': (row.get('BioProject') or row.get('ProjectID') or row.get('SRAStudy') or '').strip(),
                    'species': row.get('ScientificName', '').strip() or row.get('Organism', '').strip(),
                    'tissue': biosample_tissue or runinfo_tissue,
                    'development_stage': (row.get('DevelopmentStage') or '').strip(),
                    'library_layout': row.get('LibraryLayout', '').strip() or row.get('Layout', '').strip(),
                    'library_source': (row.get('LibrarySource') or row.get('Source') or '').strip(),
                    'library_strategy': (row.get('LibraryStrategy') or row.get('Strategy') or '').strip(),
                    'library_selection': (row.get('LibrarySelection') or row.get('Selection') or '').strip(),
                    'platform': row.get('Platform', '').strip() or row.get('Instrument', '').strip(),
                    'spots_raw': raw_spots,
                    'bases_raw': raw_bases,
                    'spots': _format_million_units(raw_spots),
                    'bases': _format_million_units(raw_bases),
                    'species_url': None,
                })
    except Exception:
        return None, None

    _attach_species_detail_links(rows)
    _save_cached_transcriptome_rows(rows)

    last_updated = None
    if os.path.isfile(TRANSCRIPTOME_CACHE_META):
        try:
            with open(TRANSCRIPTOME_CACHE_META, 'r', encoding='utf-8') as meta_handle:
                meta = json.load(meta_handle)
                last_updated = meta.get('fetched_at')
        except Exception:
            last_updated = None

    if not rows:
        return None, last_updated
    return rows, last_updated


def _extract_run_from_download_path(download_path):
    if not download_path:
        return None
    path = urllib.parse.urlparse(download_path).path
    accessions = re.findall(r'/(?:SRR|ERR|DRR)\d+(?=/|\.|$)', path)
    return accessions[-1].lstrip('/') if accessions else None


def _get_transcriptome_accession_and_run(row):
    accession = (row.get('Experiment') or '').strip() or (row.get('Accession') or '').strip()
    if not accession:
        accession = (row.get('Run') or '').strip()
    run = _extract_run_from_download_path((row.get('download_path') or '').strip())
    if not run:
        run = (row.get('Run') or '').strip()
    return accession, run


def _attach_species_detail_links(rows):
    if not rows:
        return

    species_names = {row.get('species') for row in rows if row.get('species')}
    species_map = {}
    for species_name in species_names:
        if not species_name:
            continue
        tax = get_taxonomy_for_species_name(species_name)
        if not tax:
            continue
        try:
            species_url = reverse('species_detail', args=[tax.ncbi_taxonomy_id or tax.species])
        except Exception:
            species_url = None
        species_map[species_name] = {
            'url': species_url,
            'taxid': str(tax.ncbi_taxonomy_id) if tax.ncbi_taxonomy_id else None,
            'class_name': tax.class_name,
            'order': tax.order,
            'family': tax.family,
            'updated_clade6': tax.updated_clade6,
            'updated_clade5': tax.updated_clade5,
            'updated_clade4': tax.updated_clade4,
            'updated_clade3': tax.updated_clade3,
            'updated_clade2': tax.updated_clade2,
            'updated_clade1': tax.updated_clade1,
        }

    for row in rows:
        species_entry = species_map.get(row.get('species'))
        row['species_url'] = species_entry.get('url') if species_entry else None
        row['taxid'] = species_entry.get('taxid') if species_entry else None
        row['class_name'] = species_entry.get('class_name') if species_entry else None
        row['order'] = species_entry.get('order') if species_entry else None
        row['family'] = species_entry.get('family') if species_entry else None
        row['updated_clade6'] = species_entry.get('updated_clade6') if species_entry else None
        row['updated_clade5'] = species_entry.get('updated_clade5') if species_entry else None
        row['updated_clade4'] = species_entry.get('updated_clade4') if species_entry else None
        row['updated_clade3'] = species_entry.get('updated_clade3') if species_entry else None
        row['updated_clade2'] = species_entry.get('updated_clade2') if species_entry else None
        row['updated_clade1'] = species_entry.get('updated_clade1') if species_entry else None


def _statistics_cache_needs_refresh():
    if not os.path.isfile(STATISTICS_CACHE_META):
        return True
    try:
        with open(STATISTICS_CACHE_META, 'r', encoding='utf-8') as handle:
            meta = json.load(handle)
            fetched_at = meta.get('fetched_at')
            if not fetched_at:
                return True
            ts = datetime.fromisoformat(fetched_at.replace('Z', '+00:00'))
            return ts.date() < datetime.utcnow().date()
    except Exception:
        return True


def _load_cached_statistics_summary(require_fresh=False):
    if not os.path.isfile(STATISTICS_CACHE_FILE):
        return None
    if require_fresh and _statistics_cache_needs_refresh():
        return None
    try:
        with open(STATISTICS_CACHE_FILE, 'r', encoding='utf-8') as handle:
            return json.load(handle)
    except Exception:
        return None


def _save_cached_statistics_summary(stats):
    os.makedirs(STATISTICS_CACHE_DIR, exist_ok=True)
    try:
        with open(STATISTICS_CACHE_FILE, 'w', encoding='utf-8') as handle:
            json.dump(stats, handle, ensure_ascii=False, indent=2)
        fetched_at = datetime.utcnow().isoformat() + 'Z'
        with open(STATISTICS_CACHE_META, 'w', encoding='utf-8') as meta_handle:
            json.dump({'fetched_at': fetched_at}, meta_handle, ensure_ascii=False, indent=2)
    except Exception:
        pass
    return stats


def _load_cached_ncbi_omics_counts():
    if not os.path.isfile(NCBI_OMICS_CACHE_FILE):
        return None, None
    try:
        with open(NCBI_OMICS_CACHE_FILE, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
    except Exception:
        return None, None

    last_updated = None
    if os.path.isfile(NCBI_OMICS_CACHE_META):
        try:
            with open(NCBI_OMICS_CACHE_META, 'r', encoding='utf-8') as meta_handle:
                meta = json.load(meta_handle)
                last_updated = meta.get('fetched_at')
        except Exception:
            last_updated = None
    return data, last_updated


def _save_cached_ncbi_omics_counts(data):
    os.makedirs(NCBI_OMICS_CACHE_DIR, exist_ok=True)
    with open(NCBI_OMICS_CACHE_FILE, 'w', encoding='utf-8') as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
    fetched_at = datetime.utcnow().isoformat() + 'Z'
    with open(NCBI_OMICS_CACHE_META, 'w', encoding='utf-8') as meta_handle:
        json.dump({'fetched_at': fetched_at}, meta_handle, ensure_ascii=False, indent=2)
    return fetched_at


def _ncbi_omics_cache_needs_refresh():
    if not os.path.isfile(NCBI_OMICS_CACHE_META):
        return True
    try:
        with open(NCBI_OMICS_CACHE_META, 'r', encoding='utf-8') as handle:
            meta = json.load(handle)
            fetched_at = meta.get('fetched_at')
            if not fetched_at:
                return True
            ts = datetime.fromisoformat(fetched_at.replace('Z', '+00:00'))
            return ts.date() < datetime.utcnow().date()
    except Exception:
        return True


def _ncbi_request_delay():
    try:
        return float(os.environ.get('NCBI_REQUEST_DELAY', 0.34))
    except (TypeError, ValueError):
        return 0.34


def _ncbi_omics_cache_matches_taxonomies(data):
    if not isinstance(data, dict) or not isinstance(data.get('counts'), dict):
        return False
    taxonomy_ids = {
        str(ncbi_id).strip()
        for ncbi_id in Taxonomy.objects.exclude(ncbi_taxonomy_id__isnull=True)
        .values_list('ncbi_taxonomy_id', flat=True)
        if str(ncbi_id).strip()
    }
    return taxonomy_ids.issubset(set(data['counts']))


def _get_ncbi_omics_cache(force_refresh=False):
    data, last_updated = _load_cached_ncbi_omics_counts()
    if not force_refresh:
        if data is not None and _ncbi_omics_cache_matches_taxonomies(data):
            return data

    if force_refresh or _ncbi_omics_cache_needs_refresh():
        try:
            return _build_ncbi_omics_cache()
        except Exception:
            pass

    if data is not None:
        return data
    return _build_ncbi_omics_cache()


def _build_ncbi_omics_cache():
    taxonomies = Taxonomy.objects.exclude(ncbi_taxonomy_id__isnull=True).values_list('ncbi_taxonomy_id', flat=True)
    counts = {}
    dna_species_count = 0
    rna_species_count = 0
    for ncbi_id in taxonomies:
        try:
            ncbi_id_val = int(str(ncbi_id).strip()) if ncbi_id else None
        except (TypeError, ValueError):
            continue
        if not ncbi_id_val:
            continue

        dna_runs = _ncbi_sra_run_accessions(_ncbi_sra_dna_query(ncbi_id_val))
        rna_runs = _ncbi_sra_run_accessions(_ncbi_sra_rna_query(ncbi_id_val))

        if dna_runs:
            dna_species_count += 1
        if rna_runs:
            rna_species_count += 1

        counts[str(ncbi_id)] = {
            'dna_runs': dna_runs,
            'rna_runs': rna_runs,
            'dna_count': len(dna_runs),
            'rna_count': len(rna_runs),
        }
    data = {
        'counts': counts,
        'dna_species_count': dna_species_count,
        'rna_species_count': rna_species_count,
    }
    _save_cached_ncbi_omics_counts(data)
    return data


def _get_ncbid_omics_for_taxid(ncbi_id):
    cache_data = _get_ncbi_omics_cache()
    return cache_data.get('counts', {}).get(str(ncbi_id), {'dna_runs': [], 'rna_runs': [], 'dna_count': 0, 'rna_count': 0})


def _has_ncbi_dna(ncbi_id, fallback=False):
    if not ncbi_id:
        return fallback
    cache_data, _ = _load_cached_ncbi_omics_counts()
    if cache_data:
        return bool(cache_data.get('counts', {}).get(str(ncbi_id), {}).get('dna_count', 0)) or fallback
    return fallback


def _has_ncbi_rna(ncbi_id, fallback=False):
    if not ncbi_id:
        return fallback
    cache_data, _ = _load_cached_ncbi_omics_counts()
    if cache_data:
        return bool(cache_data.get('counts', {}).get(str(ncbi_id), {}).get('rna_count', 0)) or fallback
    return fallback


def _parse_obis_depth_value(record):
    species_names = {row['species'] for row in rows if row.get('species')}
    if not species_names:
        return

    species_map = {}
    for species_name in species_names:
        if not species_name:
            continue
        tax = get_taxonomy_for_species_name(species_name)
        if not tax:
            continue
        try:
            species_map[species_name] = {
                'url': reverse('species_detail', args=[tax.ncbi_taxonomy_id or tax.species]),
                'taxid': str(tax.ncbi_taxonomy_id) if tax.ncbi_taxonomy_id else None,
            }
        except Exception:
            species_map[species_name] = {'url': None, 'taxid': str(tax.ncbi_taxonomy_id) if tax.ncbi_taxonomy_id else None}

    for row in rows:
        species_entry = species_map.get(row.get('species'))
        row['species_url'] = species_entry.get('url') if species_entry else None
        row['taxid'] = species_entry.get('taxid') if species_entry else None


def _normalize_species_name(name):
    if not name:
        return ''
    normalized = re.sub(r'[^a-z0-9 ]+', '', str(name).strip().lower())
    return re.sub(r'\s+', ' ', normalized).strip()


def _filter_transcriptome_rows_by_taxid(rows, taxonomy_id):
    if not taxonomy_id:
        return rows, None

    tax = get_taxonomy_by_identifier(taxonomy_id)
    if tax and tax.ncbi_taxonomy_id:
        filtered_rows = [row for row in rows if row.get('taxid') == str(tax.ncbi_taxonomy_id)]
        if filtered_rows:
            return filtered_rows, tax

    # If taxonomy identifier is a species or NCBI ID lookup failed, fall back to species matching.
    taxonomy_name = taxonomy_id.strip()
    tax_species_norm = _normalize_species_name(tax.species) if tax else _normalize_species_name(taxonomy_name)
    filtered_rows = []
    for row in rows:
        row_species = row.get('species')
        if not row_species:
            continue
        row_species_norm = _normalize_species_name(row_species)
        if row_species_norm == tax_species_norm:
            filtered_rows.append(row)
            continue
        if tax_species_norm and row_species_norm:
            if tax_species_norm in row_species_norm or row_species_norm in tax_species_norm:
                filtered_rows.append(row)
    return filtered_rows, tax


def _is_rna_transcriptome_row(row):
    if not row:
        return False
    candidate_fields = [
        row.get('library_source'),
        row.get('library_strategy'),
        row.get('library_selection'),
        row.get('library_layout'),
    ]
    for value in candidate_fields:
        if not value:
            continue
        if re.search(r'\b(rna|transcriptom|mrna|cdna)\b', str(value), re.I):
            return True
    if not any(row.get(field) for field in ('library_source', 'library_strategy', 'library_selection')):
        return True
    return False


@lru_cache(maxsize=1)
def _cached_rna_counts_by_species():
    rows, _ = _load_cached_transcriptome_rows()
    if not rows:
        return {}

    counts = Counter()
    for row in rows:
        if not _is_rna_transcriptome_row(row):
            continue
        species_key = _normalize_species_name(row.get('species'))
        run = str(row.get('run') or row.get('study') or '').strip()
        if species_key and run:
            counts[species_key] += 1
    return counts


@lru_cache(maxsize=1)
def _cached_dna_counts():
    by_taxid = Counter()
    by_species = Counter()
    if not os.path.isfile(DNA_SRA_CACHE_FILE):
        return by_taxid, by_species

    seen = set()
    try:
        with open(DNA_SRA_CACHE_FILE, 'r', encoding='utf-8', errors='replace', newline='') as handle:
            for row in csv.DictReader(handle):
                run = (row.get('Run') or row.get('run') or '').strip()
                taxid = (row.get('TaxID') or row.get('taxid') or '').strip()
                species = _normalize_species_name(row.get('ScientificName') or row.get('species'))
                key = (taxid, species, run)
                if not run or key in seen:
                    continue
                seen.add(key)
                if taxid:
                    by_taxid[taxid] += 1
                if species:
                    by_species[species] += 1
    except (OSError, csv.Error):
        return Counter(), Counter()
    return by_taxid, by_species


def _get_dna_species_count(ncbi_cache=None):
    cached_dna_by_taxid, _ = _cached_dna_counts()
    if cached_dna_by_taxid:
        return len(cached_dna_by_taxid)
    if ncbi_cache is None:
        ncbi_cache = _get_ncbi_omics_cache()
    return ncbi_cache.get('dna_species_count', 0)


def _query_transcriptome_rows(rows, search_value, study_filter, run_filter, platform_filter, layout_filter):
    if not rows:
        return []

    def norm(value):
        return str(value or '').strip().lower()

    search_value = norm(search_value)
    study_filter = norm(study_filter)
    run_filter = norm(run_filter)
    platform_filter = norm(platform_filter)
    layout_filter = norm(layout_filter)

    def matches(row):
        if study_filter and study_filter not in norm(row.get('study')):
            return False
        if run_filter and run_filter not in norm(row.get('run')):
            return False
        if platform_filter and platform_filter not in norm(row.get('platform')):
            return False
        if layout_filter and layout_filter not in norm(row.get('library_layout')):
            return False
        if search_value:
            fields = [
                row.get('run'), row.get('study'), row.get('bioproject'), row.get('species'),
                row.get('tissue'), row.get('development_stage'), row.get('library_layout'),
                row.get('library_source'), row.get('library_strategy'), row.get('library_selection'),
                row.get('platform'), row.get('class_name'), row.get('order'), row.get('family'),
                row.get('updated_clade1'), row.get('updated_clade2'), row.get('updated_clade3'),
                row.get('updated_clade4'), row.get('updated_clade5'), row.get('updated_clade6'),
            ]
            combined = ' '.join(norm(value) for value in fields)
            if search_value not in combined:
                return False
        return True

    return [row for row in rows if matches(row)]


def _map_transcriptome_order_column(index):
    mapping = {
        1: 'run',
        2: 'study',
        3: 'class_name',
        4: 'order',
        5: 'family',
        6: 'updated_clade6',
        7: 'updated_clade5',
        8: 'updated_clade4',
        9: 'updated_clade3',
        10: 'updated_clade2',
        11: 'updated_clade1',
        12: 'species',
        13: 'tissue',
        14: 'development_stage',
        15: 'platform',
        16: 'library_layout',
        17: 'spots_raw',
        18: 'bases_raw',
    }
    return mapping.get(index)


def _transcriptome_data_rows(rows):
    data = []
    for row in rows:
        data.append([
            '',
            row.get('run', ''),
            row.get('study', ''),
            row.get('class_name', '-') or '-',
            row.get('order', '-') or '-',
            row.get('family', '-') or '-',
            row.get('updated_clade6', '-') or '-',
            row.get('updated_clade5', '-') or '-',
            row.get('updated_clade4', '-') or '-',
            row.get('updated_clade3', '-') or '-',
            row.get('updated_clade2', '-') or '-',
            row.get('updated_clade1', '-') or '-',
            row.get('species', ''),
            row.get('tissue', ''),
            row.get('development_stage', ''),
            row.get('platform', ''),
            row.get('library_layout', ''),
            row.get('spots', ''),
            row.get('bases', ''),
        ])
    return data


def transcriptome_data(request):
    params = request.GET if request.method == 'GET' else request.POST
    try:
        draw = int(params.get('draw', '1') or '1')
    except ValueError:
        draw = 1
    try:
        start = int(params.get('start', '0') or '0')
    except ValueError:
        start = 0
    try:
        length = int(params.get('length', '10') or '10')
    except ValueError:
        length = 10
    if length < 0:
        length = 10

    search_value = params.get('search[value]', '').strip()
    study_filter = params.get('columns[2][search][value]', '').strip()
    run_filter = params.get('columns[1][search][value]', '').strip()
    platform_filter = params.get('columns[15][search][value]', '').strip()
    layout_filter = params.get('columns[16][search][value]', '').strip()
    selected_taxid = params.get('taxid', '').strip()

    transcriptome_rows, last_updated = _load_cached_transcriptome_rows()
    if transcriptome_rows is None:
        transcriptome_rows = _get_sample_transcriptome_rows()
        for row in transcriptome_rows:
            row['spots'] = _format_million_units(row.get('spots_raw'))
            row['bases'] = _format_million_units(row.get('bases_raw'))

    filtered_rows, _ = _filter_transcriptome_rows_by_taxid(transcriptome_rows, selected_taxid)
    filtered_rows = [row for row in filtered_rows if _is_rna_transcriptome_row(row)]

    records_total = len(filtered_rows)
    filtered_rows = _query_transcriptome_rows(filtered_rows, search_value, study_filter, run_filter, platform_filter, layout_filter)
    records_filtered = len(filtered_rows)

    order_col = params.get('order[0][column]')
    order_dir = params.get('order[0][dir]', 'asc')
    if order_col is not None:
        try:
            order_idx = int(order_col)
            order_key = _map_transcriptome_order_column(order_idx)
            if order_key:
                filtered_rows.sort(key=lambda row: str(row.get(order_key) or '').lower(), reverse=(order_dir == 'desc'))
        except ValueError:
            pass

    page_rows = filtered_rows[start:start + length]
    data = _transcriptome_data_rows(page_rows)

    return JsonResponse({
        'draw': draw,
        'recordsTotal': records_total,
        'recordsFiltered': records_filtered,
        'data': data,
    })


@cache_page(300)
def transcriptome_page(request):
    transcriptome_rows, last_updated = _load_cached_transcriptome_rows()
    source = 'cached'
    source_note = None

    if transcriptome_rows is None:
        transcriptome_rows = _get_sample_transcriptome_rows()
        for row in transcriptome_rows:
            row['spots'] = _format_million_units(row.get('spots_raw'))
            row['bases'] = _format_million_units(row.get('bases_raw'))
        last_updated = None

    selected_taxid = request.GET.get('taxid')
    filtered_rows, selected_taxonomy = _filter_transcriptome_rows_by_taxid(transcriptome_rows, selected_taxid)
    transcriptome_rows = [row for row in filtered_rows if _is_rna_transcriptome_row(row)]

    total_runs = len(transcriptome_rows)
    unique_studies = len({row['bioproject'] for row in transcriptome_rows if row.get('bioproject')})
    unique_species = len({row['species'] for row in transcriptome_rows if row.get('species')})
    unique_tissues = len({row['tissue'] for row in transcriptome_rows if row.get('tissue')})
    total_spots_raw = sum(row.get('spots_raw', 0) for row in transcriptome_rows)

    study_values = sorted({row['study'] for row in transcriptome_rows if row.get('study')})
    platform_values = sorted({row['platform'] for row in transcriptome_rows if row.get('platform')})
    layout_values = sorted({row['library_layout'] for row in transcriptome_rows if row.get('library_layout')})

    return render(request, 'transcriptome.html', {
        'page_title': 'Transcriptome',
        'active_page': 'transcriptome',
        'datasets': transcriptome_rows,
        'summary': {
            'runs': total_runs,
            'studies': unique_studies,
            'species': unique_species,
            'tissues': unique_tissues,
            'total_spots': _format_million_units(total_spots_raw),
        },
        'selected_taxonomy': selected_taxonomy,
        'selected_taxid': selected_taxid,
        'last_updated': last_updated,
        'source': source,
        'source_note': source_note,
        'study_values': study_values,
        'platform_values': platform_values,
        'layout_values': layout_values,
    })


def _parse_integer_string(value):
    if value is None:
        return 0
    try:
        return int(re.sub(r'[^0-9]', '', str(value)))
    except ValueError:
        return 0


def _format_million_units(value):
    if value is None:
        return '-'
    try:
        value = float(value)
        millions = value / 1_000_000
        if millions.is_integer():
            return f"{int(millions)} M"
        return f"{millions:.1f} M"
    except (TypeError, ValueError):
        return value or '-'


def _format_integer(value):
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return value or '-'


def _format_count_if_large(value):
    try:
        number = int(value)
        return f"{number:,}" if abs(number) >= 1000 else str(number)
    except (TypeError, ValueError):
        return value or '-'


def singlecell_page(request):
    datasets = SingleCell.objects.select_related('taxonomy').order_by('species')
    species_count = datasets.values('species').distinct().count()
    platform_count = datasets.exclude(platform__isnull=True).exclude(platform__exact='').values('platform').distinct().count()
    total_cells = sum(_parse_integer_string(item.cells) for item in datasets)

    return render(request, 'singlecell.html', {
        'page_title': 'Single-cell',
        'active_page': 'singlecell',
        'datasets': datasets,
        'summary': {
            'species_count': species_count,
            'datasets': datasets.count(),
            'platforms': platform_count,
            'total_cells': _format_integer(total_cells),
        }
    })


def _get_synteny_upload_root():
    media_root = getattr(settings, 'MEDIA_ROOT', os.path.join(settings.BASE_DIR, 'media'))
    upload_root = os.path.join(media_root, 'synteny_uploads')
    os.makedirs(upload_root, exist_ok=True)
    return upload_root


import errno

def _is_process_running(pid):
    if pid is None:
        return False
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False

    try:
        waited_pid, status = os.waitpid(pid, os.WNOHANG)
        if waited_pid == 0:
            return True
        return False
    except ChildProcessError:
        pass
    except OSError as exc:
        if exc.errno in {errno.ECHILD, errno.ESRCH}:
            return False
        if exc.errno == errno.EPERM:
            return True
        return False

    try:
        os.kill(pid, 0)
    except OSError as exc:
        if exc.errno == errno.ESRCH:
            return False
        if exc.errno == errno.EPERM:
            return True
        return False
    return True


def _find_oxford_plot_seqids_file(request_dir):
    for name in sorted(os.listdir(request_dir)):
        if name.endswith('.oxford_plot.seqids'):
            return name
    return None


def _run_synteny_upload_task(request_dir, command, log_path, status_path, species_raw, species_slug, email_address=None, page_url=None):
    email_address = str(email_address or '').strip()
    page_url = str(page_url or '').strip()
    try:
        with open(status_path, 'w', encoding='utf-8') as status_file:
            json.dump({
                'status': 'running',
                'pid': None,
                'created': datetime.utcnow().isoformat(),
                'species': species_raw,
                'species_slug': species_slug,
                'page_url': page_url,
                'email_address': email_address,
            }, status_file)

        with open(log_path, 'a', encoding='utf-8') as log_file:
            proc = subprocess.Popen(
                command,
                cwd=request_dir,
                stdout=log_file,
                stderr=log_file,
                text=True,
            )

        with open(status_path, 'w', encoding='utf-8') as status_file:
            json.dump({
                'status': 'running',
                'pid': proc.pid,
                'created': datetime.utcnow().isoformat(),
                'species': species_raw,
                'species_slug': species_slug,
                'page_url': page_url,
                'email_address': email_address,
            }, status_file)

        proc.wait()
        final_status = 'finished' if proc.returncode == 0 else 'failed'
        status_payload = {
            'status': final_status,
            'pid': proc.pid,
            'created': datetime.utcnow().isoformat(),
            'completed': datetime.utcnow().isoformat(),
            'species': species_raw,
            'species_slug': species_slug,
            'page_url': page_url,
            'email_address': email_address,
            'completion_email_sent': False,
        }
        if email_address:
            subject = 'Your AnnelidaDB synteny job completed' if final_status == 'finished' else 'Your AnnelidaDB synteny job failed'
            body_parts = [
                f'Your synteny upload for {species_raw} has {"completed successfully" if final_status == "finished" else "finished with an error"}.',
                f'Status: {final_status}',
            ]
            if page_url:
                body_parts.append(f'View job details: {page_url}')
            if final_status != 'finished':
                body_parts.append('Please check the page for logs and any error details.')
            _send_notification_email(subject, '\n\n'.join(body_parts), email_address)
            status_payload['completion_email_sent'] = True
        with open(status_path, 'w', encoding='utf-8') as status_file:
            json.dump(status_payload, status_file)

        if email_address:
            subject = 'Your AnnelidaDB synteny job completed' if final_status == 'finished' else 'Your AnnelidaDB synteny job failed'
            body_parts = [
                f'Your synteny upload for {species_raw} has {"completed successfully" if final_status == "finished" else "finished with an error"}.',
                f'Status: {final_status}',
            ]
            if page_url:
                body_parts.append(f'View job details: {page_url}')
            if final_status != 'finished':
                body_parts.append('Please check the page for logs and any error details.')
            _send_notification_email(subject, '\n\n'.join(body_parts), email_address)
    except Exception as exc:
        with open(status_path, 'w', encoding='utf-8') as status_file:
            json.dump({
                'status': 'failed',
                'pid': None,
                'created': datetime.utcnow().isoformat(),
                'completed': datetime.utcnow().isoformat(),
                'species': species_raw,
                'species_slug': species_slug,
                'error': str(exc),
                'page_url': page_url,
                'email_address': email_address,
            }, status_file)
        if email_address:
            subject = 'Your AnnelidaDB synteny job failed'
            body = (
                f'Your synteny upload for {species_raw} failed during processing.\n\n'
                f'Error: {str(exc)}\n\n'
                f'View job details: {page_url}\n'
            )
            _send_notification_email(subject, body, email_address)
            try:
                with open(status_path, 'r', encoding='utf-8') as status_file:
                    existing = json.load(status_file)
            except Exception:
                existing = {}
            existing['completion_email_sent'] = True
            try:
                with open(status_path, 'w', encoding='utf-8') as status_file:
                    json.dump(existing, status_file)
            except Exception:
                pass


def _get_synteny_upload_results(token, request):
    upload_root = _get_synteny_upload_root()
    request_dir = os.path.join(upload_root, token)
    if not os.path.isdir(request_dir):
        return None

    status_path = os.path.join(request_dir, 'status.json')
    status_data = {}
    if os.path.exists(status_path):
        try:
            with open(status_path, 'r', encoding='utf-8') as status_file:
                status_data = json.load(status_file)
        except Exception:
            status_data = {}

    status = status_data.get('status', 'unknown')
    if status == 'completed':
        status = 'finished'
        status_data['status'] = status
    pid = status_data.get('pid')
    if status == 'running' and pid:
        if not _is_process_running(pid):
            status = 'finished'
            status_data['status'] = status
            try:
                with open(status_path, 'w', encoding='utf-8') as status_file:
                    json.dump(status_data, status_file)
            except Exception:
                pass

    uploaded_files = set(status_data.get('uploaded_files', []))
    renamed_files = set(status_data.get('renamed_files', []))
    result_files = []
    log_url = None
    seqids_file_name = None
    for name in sorted(os.listdir(request_dir)):
        if name == 'status.json':
            continue
        file_url = request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{token}/{urllib.parse.quote(name)}")
        if name == 'run.log':
            log_url = file_url
            continue
        if name in {'species.pep', 'species.bed'} or name in uploaded_files or name in renamed_files:
            continue
        if name.lower().endswith('.tmp') or name.lower().endswith('.trait'):
            continue
        if name.endswith('.oxford_plot.seqids'):
            seqids_file_name = name
        result_files.append({
            'name': name,
            'url': file_url,
        })

    if status in {'finished', 'failed'} and status_data.get('email_address') and not status_data.get('completion_email_sent'):
        subject = 'Your AnnelidaDB synteny job completed' if status == 'finished' else 'Your AnnelidaDB synteny job failed'
        body_parts = [
            f'Your synteny upload for {status_data.get("species", "your species")} has {"completed successfully" if status == "finished" else "finished with an error"}.',
            f'Status: {status}',
        ]
        if status_data.get('page_url'):
            body_parts.append(f'View job details: {status_data.get("page_url")}')
        if status != 'finished':
            body_parts.append('Please check the page for logs and any error details.')
        _send_notification_email(subject, '\n\n'.join(body_parts), status_data['email_address'])
        status_data['completion_email_sent'] = True
        try:
            with open(status_path, 'w', encoding='utf-8') as status_file:
                json.dump(status_data, status_file)
        except Exception:
            pass
    return {
        'token': token,
        'base_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{token}/"),
        'page_url': request.build_absolute_uri(f"{request.path}?token={token}"),
        'files': result_files,
        'status': status,
        'log_url': log_url,
        'seqids_file_name': seqids_file_name,
        'seqids_edit_url': request.build_absolute_uri(f"/synteny-upload/edit-seqids/?token={token}") if seqids_file_name else None,
        'stdout': status_data.get('stdout'),
        'stderr': status_data.get('stderr'),
        'species': status_data.get('species'),
        'species_slug': status_data.get('species_slug'),
        'uploaded_files': status_data.get('original_uploaded_files', status_data.get('uploaded_files', [])),
        'renamed_files': status_data.get('renamed_files', []),
    }


def singlecell_detail_page(request, pk):
    dataset = get_object_or_404(SingleCell.objects.select_related('taxonomy'), pk=pk)
    marker_genes = list(dataset.marker_genes.all())

    formatted_species = _format_species(dataset.species)
    umap_datasets = _find_singlecell_umap_files(formatted_species, [dataset])
    umap_dataset = umap_datasets[0] if umap_datasets else None
    umap_file_path = _find_singlecell_umap_file(formatted_species, dataset)
    cluster_image_path = _find_singlecell_cluster_image_file(formatted_species, dataset)
    markers_file_path = _find_singlecell_markers_file(formatted_species, dataset)
    markers_file_url = _format_static_url(markers_file_path) if markers_file_path else None
    expression_file_path = _find_singlecell_expression_file(formatted_species, dataset)
    expression_file_url = _format_static_url(expression_file_path) if expression_file_path else None
    rds_files = [
        {
            'filename': os.path.basename(path),
            'url': _format_static_url(path),
        }
        for path, _ in _collect_singlecell_rds_files(dataset)
    ]
    expression_folder_url = None
    expression_folder_ready = False
    expression_generation_queued = False
    if expression_file_url:
        if expression_file_url.endswith('/Expression.csv'):
            expression_folder_url = expression_file_url[:-len('/Expression.csv')] + '/Expression/'
        else:
            expression_folder_url = os.path.dirname(expression_file_url) + '/Expression/'
        expression_folder_path = os.path.join(os.path.dirname(expression_file_path), 'Expression')
        expression_folder_ready = os.path.isdir(expression_folder_path) and any(
            os.path.isfile(os.path.join(expression_folder_path, f)) and f.lower().endswith('.csv')
            for f in os.listdir(expression_folder_path)
        )

    if expression_file_path:
        expression_folder_path = os.path.join(os.path.dirname(expression_file_path), 'Expression')
        if not os.path.isdir(expression_folder_path) or not any(
            os.path.isfile(os.path.join(expression_folder_path, f))
            for f in os.listdir(expression_folder_path)
        ):
            _enqueue_singlecell_expression_generation(expression_file_path, umap_file_path)
            expression_generation_queued = True

    return render(request, 'singlecell_detail.html', {
        'page_title': f"{dataset.species} Single-cell Detail",
        'active_page': 'singlecell',
        'dataset': dataset,
        'marker_genes': marker_genes,
        'umap_dataset': umap_dataset,
        'cluster_image_url': _format_static_url(cluster_image_path) if cluster_image_path else None,
        'cluster_image_filename': os.path.basename(cluster_image_path) if cluster_image_path else None,
        'has_umap_coordinates': bool(umap_file_path),
        'markers_file_url': markers_file_url,
        'expression_file_filename': os.path.basename(expression_file_path) if expression_file_path else None,
        'expression_file_url': expression_file_url,
        'rds_files': rds_files,
        'expression_folder_url': expression_folder_url,
        'expression_folder_ready': expression_folder_ready,
        'expression_generation_queued': expression_generation_queued,
    })


def synteny_page(request):
    synteny_dir = os.path.join(settings.BASE_DIR, 'static', 'plots', 'Synteny')
    plot_suffix = '.AnLGs.oxford_plot.png'
    fallback_suffix = '.AnLGs.oxford_plot.pdf'
    species_pairs = {}
    if os.path.exists(synteny_dir):
        for filename in sorted(os.listdir(synteny_dir)):
            if not filename.endswith(plot_suffix):
                continue
            pair_root = filename[:-len(plot_suffix)]
            if '.' not in pair_root:
                continue
            species1, species2 = pair_root.split('.', 1)
            species_pairs.setdefault(species1, []).append({
                'species2': species2,
                'filename': filename,
                'url': f'/static/plots/Synteny/{filename}'
            })
        if not species_pairs:
            for filename in sorted(os.listdir(synteny_dir)):
                if not filename.endswith(fallback_suffix):
                    continue
                pair_root = filename[:-len(fallback_suffix)]
                if '.' not in pair_root:
                    continue
                species1, species2 = pair_root.split('.', 1)
                species_pairs.setdefault(species1, []).append({
                    'species2': species2,
                    'filename': filename,
                    'url': f'/static/plots/Synteny/{filename}'
                })
    species1_options = sorted(species_pairs.keys())

    species2_unique = sorted({pair['species2'] for pairs in species_pairs.values() for pair in pairs})
    plot_file_count = sum(len(v) for v in species_pairs.values())
    unique_species_count = len(set(species1_options).union(species2_unique))
    synteny_summary = [
        {
            'title': 'Oxford Plots',
            'value': plot_file_count,
            'border_class': 'border-info',
            'subtitle': 'Available Oxford plot files for synteny comparisons',
        },
        {
            'title': 'Reference species',
            'value': len(species1_options),
            'border_class': 'border-success',
            'subtitle': 'Unique first-species references',
        },
        {
            'title': 'Query species',
            'value': len(species2_unique),
            'border_class': 'border-warning',
            'subtitle': 'Unique second-species targets',
        },
        {
            'title': 'Total species',
            'value': unique_species_count,
            'border_class': 'border-primary',
            'subtitle': 'Distinct species represented across all plots',
        },
    ]

    upload_info = None
    upload_error = None
    show_upload_info = False
    token = request.GET.get('token')
    if token and request.method != 'POST':
        upload_info = _get_synteny_upload_results(token, request)
        if upload_info is None:
            upload_error = 'The requested synteny result token was not found.'

    if request.method == 'POST':
        show_upload_info = True
        pep_file = request.FILES.get('species_pep')
        bed_file = request.FILES.get('species_bed')
        if not pep_file or not bed_file:
            upload_error = 'Please upload both species.pep and species.bed files.'
        else:
            upload_root = _get_synteny_upload_root()
            upload_token = uuid.uuid4().hex
            request_dir = os.path.join(upload_root, upload_token)
            os.makedirs(request_dir, exist_ok=True)

            pep_path = os.path.join(request_dir, 'species.pep')
            bed_path = os.path.join(request_dir, 'species.bed')
            with open(pep_path, 'wb') as f:
                for chunk in pep_file.chunks():
                    f.write(chunk)
            with open(bed_path, 'wb') as f:
                for chunk in bed_file.chunks():
                    f.write(chunk)

            script_path = '/home/yunlongli/Annelida/Synteny/Annelida-upload.py'
            command = ['python3', script_path, 'species']
            log_path = os.path.join(request_dir, 'run.log')
            status_path = os.path.join(request_dir, 'status.json')
            try:
                with open(status_path, 'w', encoding='utf-8') as status_file:
                    json.dump({'status': 'running', 'pid': None, 'created': datetime.utcnow().isoformat()}, status_file)

                with open(log_path, 'a', encoding='utf-8') as log_file:
                    proc = subprocess.Popen(
                        command,
                        cwd=request_dir,
                        stdout=log_file,
                        stderr=log_file,
                        text=True,
                    )

                with open(status_path, 'w', encoding='utf-8') as status_file:
                    json.dump({'status': 'running', 'pid': proc.pid, 'created': datetime.utcnow().isoformat()}, status_file)

                log_url = request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/run.log")
                upload_info = {
                    'token': upload_token,
                    'base_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/"),
                    'page_url': request.build_absolute_uri(f"{request.path}?token={upload_token}"),
                    'files': [],
                    'status': 'running',
                    'log_url': log_url,
                }
            except Exception as exc:
                upload_error = f'Unable to start background synteny comparison: {str(exc)}'
                upload_info = {
                    'token': upload_token,
                    'base_url': request.build_absolute_uri(f"{settings.MEDIA_URL}synteny_uploads/{upload_token}/"),
                    'page_url': request.build_absolute_uri(f"{request.path}?token={upload_token}"),
                    'files': [],
                    'status': 'failed',
                }

    selected_species1 = request.GET.get('species1')
    selected_species2 = request.GET.get('species2')

    if not selected_species1 or selected_species1 not in species_pairs:
        selected_species1 = species1_options[0] if species1_options else None

    species2_options = []
    if selected_species1:
        species2_options = [pair['species2'] for pair in species_pairs[selected_species1]]

    if selected_species2 not in species2_options:
        selected_species2 = species2_options[0] if species2_options else None

    def lookup_taxonomy(species_key):
        if not species_key:
            return None
        candidate = species_key.replace('_', ' ')
        tax = get_taxonomy_for_species_name(candidate)
        if tax:
            return tax
        alt_candidate = candidate.replace(' sp ', ' sp. ')
        return get_taxonomy_for_species_name(alt_candidate)

    selected_species1_tax = lookup_taxonomy(selected_species1)
    selected_species2_tax = lookup_taxonomy(selected_species2)

    plot_url = None
    active_filename = None
    if selected_species1 and selected_species2:
        for pair in species_pairs.get(selected_species1, []):
            if pair['species2'] == selected_species2:
                plot_url = pair['url']
                active_filename = pair['filename']
                break

    # Get Oxford plot URLs for species 1 and 2
    oxford_url1 = None
    oxford_url2 = None
    if selected_species1_tax:
        genome1 = selected_species1_tax.genomes.first()
        formatted_species1 = _format_species(selected_species1_tax.species)
        _, oxford_url1 = _get_alg_oxford_urls(genome1, formatted_species1)
    
    if selected_species2_tax:
        genome2 = selected_species2_tax.genomes.first()
        formatted_species2 = _format_species(selected_species2_tax.species)
        _, oxford_url2 = _get_alg_oxford_urls(genome2, formatted_species2)

    return render(request, 'synteny.html', {
        'page_title': 'Synteny Oxford Plot',
        'species1_options': species1_options,
        'species2_options': species2_options,
        'selected_species1': selected_species1,
        'selected_species2': selected_species2,
        'selected_species1_tax': selected_species1_tax,
        'selected_species2_tax': selected_species2_tax,
        'oxford_url1': oxford_url1,
        'oxford_url2': oxford_url2,
        'plot_url': plot_url,
        'active_filename': active_filename,
        'species_pairs_json': json.dumps({k: [pair['species2'] for pair in v] for k, v in species_pairs.items()}),
        'synteny_summary': synteny_summary,
        'upload_info': upload_info,
        'upload_error': upload_error,
        'show_upload_info': show_upload_info,
    })


def _resolve_taxonomy_identifier(value):
    if value is None:
        return None
    if str(value).isdigit():
        return Taxonomy.objects.filter(ncbi_taxonomy_id=value).first()
    return get_taxonomy_for_species_name(value)


def _expand_file_patterns(patterns):
    expanded = []
    for pattern in patterns:
        expanded.append(pattern)
        if not pattern.endswith('.gz'):
            expanded.append(f"{pattern}.gz")
    return expanded


def _collect_genome_files(tax):
    if not tax:
        return []
    formatted_species = _format_species(tax.species)
    genome_data_dir = _species_genome_dir(formatted_species)
    base_patterns = [
        '*.genomic.fasta.masked', '*.genomic.fasta', '*.fa', '*.fasta',
        '*.gff3', '*.gtf', '*.pep.fasta', '*.pep.faa', '*.CDS.fasta', '*.bed',
        '*-families.fa.gz', '*.genomic.fasta.out.gz', '*.genomic.fasta.tbl',
    ]
    patterns = _expand_file_patterns(base_patterns)

    files = []
    seen_files = set()
    if os.path.exists(genome_data_dir):
        for pattern in patterns:
            for filepath in sorted(glob.glob(os.path.join(genome_data_dir, pattern))):
                if os.path.isfile(filepath) and filepath not in seen_files:
                    seen_files.add(filepath)
                    files.append(filepath)
        if not files:
            for filepath in sorted(os.listdir(genome_data_dir)):
                abs_path = os.path.join(genome_data_dir, filepath)
                if os.path.isfile(abs_path) and abs_path not in seen_files:
                    seen_files.add(abs_path)
                    files.append(abs_path)
    base_data_dir = os.path.join(settings.BASE_DIR, 'static', 'data')
    return [(path, os.path.relpath(path, base_data_dir)) for path in files]


def _collect_mitogenome_files(tax):
    if not tax:
        return []
    formatted_species = _format_species(tax.species)
    folder_name = f"mitoz.{formatted_species}.mt.fasta.result"
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, folder_name)
    base_data_dir = os.path.join(settings.BASE_DIR, 'static', 'data')
    files = []
    if os.path.exists(base_dir):
        for root, _, filenames in os.walk(base_dir):
            for filename in sorted(filenames):
                abs_path = os.path.join(root, filename)
                if os.path.isfile(abs_path):
                    files.append((abs_path, os.path.relpath(abs_path, base_data_dir)))
    return files


def _resolve_static_data_file_path(rel_path):
    if not rel_path:
        return None
    safe_rel = os.path.normpath(rel_path).lstrip(os.sep)
    static_root = os.path.abspath(os.path.join(settings.BASE_DIR, 'static'))
    abs_path = os.path.abspath(os.path.join(static_root, safe_rel))
    if not abs_path.startswith(static_root + os.sep) and abs_path != static_root:
        return None
    return abs_path if os.path.isfile(abs_path) else None


def static_data_download(request, file_path):
    abs_path = _resolve_static_data_file_path(file_path)
    if not abs_path:
        return HttpResponseNotFound('Requested download file not found.')
    filename = os.path.basename(abs_path)
    content_type, _ = mimetypes.guess_type(filename)
    return FileResponse(
        open(abs_path, 'rb'),
        as_attachment=True,
        filename=filename,
        content_type=content_type or 'application/octet-stream',
    )


def _symbiont_directory_candidates(species, sub=None):
    if not species:
        return []
    formatted_species = _format_species(species)
    sub_path = os.path.normpath(str(sub or 'symbiont').strip())
    if sub_path in ('', '.'):
        sub_path = 'symbiont'
    if sub_path == '..' or sub_path.startswith('../') or os.path.isabs(sub_path):
        return []
    if not sub_path.startswith('symbiont'):
        sub_path = os.path.join('symbiont', sub_path)
    return [
        os.path.join(settings.BASE_DIR, 'static', 'data', candidate, sub_path)
        for candidate in (formatted_species, species)
    ]


def _split_symbiont_identifier(value):
    identifier, separator, sub = str(value or '').partition('::')
    return identifier, sub if separator else None


def _collect_symbiont_files(value, sub=None):
    value, embedded_sub = _split_symbiont_identifier(value)
    sub = sub or embedded_sub
    tax = _resolve_taxonomy_identifier(value)
    species = tax.species if tax else value
    if not species:
        return []
    files = []
    base_data_dir = os.path.join(settings.BASE_DIR, 'static', 'data')
    for base_dir in _symbiont_directory_candidates(species, sub):
        if os.path.exists(base_dir):
            for root, _, filenames in os.walk(base_dir):
                for filename in sorted(filenames):
                    abs_path = os.path.join(root, filename)
                    if os.path.isfile(abs_path):
                        files.append((abs_path, os.path.relpath(abs_path, base_data_dir)))
            break
    return files


def _collect_singlecell_rds_files(dataset):
    formatted_species = _format_species(dataset.species)
    base_data_dir = os.path.join(settings.BASE_DIR, 'static', 'data')
    candidate_dirs = _singlecell_candidate_dirs(formatted_species, dataset)
    if not candidate_dirs:
        return []

    files = []
    seen_paths = set()
    for base_dir in candidate_dirs:
        if not os.path.isdir(base_dir):
            continue
        for root, _, filenames in os.walk(base_dir):
            for filename in sorted(filenames):
                if not filename.lower().endswith('.rds'):
                    continue
                abs_path = os.path.join(root, filename)
                if abs_path in seen_paths or not os.path.isfile(abs_path):
                    continue
                seen_paths.add(abs_path)
                files.append((abs_path, os.path.relpath(abs_path, base_data_dir)))
        if getattr(dataset, 'dir', None):
            break
    return files


def _download_file_label(filename, data_type=None):
    name = filename.lower()
    if data_type == 'mitogenome':
        if name.endswith('summary.2.txt'):
            return 'Location'
        if name.endswith('cds.fasta'):
            return 'CDS'
        if name.endswith('gene.fasta'):
            return 'Gene'
        if name.endswith('cds_translation.fasta'):
            return 'Protein'
        if name.endswith('rrna.fasta'):
            return 'rRNA'
        if name.endswith('trna.fasta'):
            return 'tRNA'
        if name.endswith('.fasta') or name.endswith('.fasta.gz'):
            return 'Genome'
        return filename
    if data_type == 'singlecell':
        if name.endswith('.rds'):
            return 'RDS'
        return os.path.splitext(filename)[0]
    if data_type == 'symbiont':
        if name.endswith('.gff3') or name.endswith('.gff3.gz'):
            return 'GFF3'
        if name.endswith('.ffn') or name.endswith('.ffn.gz'):
            return 'Gene'
        if name.endswith('.faa') or name.endswith('.faa.gz'):
            return 'Protein'
        if name.endswith('.fna') or name.endswith('.fna.gz'):
            return 'Genome'
        return filename
    if data_type == 'genome':
        is_annevo = 'annevo' in name
        if 'emapper.annotations' in name or name.endswith('.annotations'):
            return 'EGGNOG Annotation (ANNEVO)' if 'annevo' in name else 'EGGNOG Annotation'
        if name.endswith('.pep.fasta.tsv') or name.endswith('.pep.fasta.tsv.gz') or name.endswith('.tsv'):
            return 'InterProScan Annotation (ANNEVO)' if 'annevo' in name else 'InterProScan Annotation'
        if name.endswith('-families.fa.gz'):
            return 'Repeat FAMILY'
        if name.endswith('.genomic.fasta.out.gz'):
            return 'Repeat LOCATION'
        if name.endswith('.genomic.fasta.tbl'):
            return 'Repeat SUMMARY'
        if name.endswith('.gff3') or name.endswith('.gff3.gz'):
            return 'Gene GFF3 (ANNEVO)' if is_annevo else 'Gene GFF3'
        if name.endswith('.gtf') or name.endswith('.gtf.gz'):
            return 'Gene GTF (ANNEVO)' if is_annevo else 'Gene GTF'
        if name.endswith('.pep.fasta') or name.endswith('.pep.fasta.gz') or name.endswith('.pep.faa') or name.endswith('.pep.faa.gz'):
            return 'Protein FASTA (ANNEVO)' if is_annevo else 'Protein FASTA'
        if name.endswith('.cds.fasta') or name.endswith('.cds.fasta.gz'):
            return 'CDS FASTA (ANNEVO)' if is_annevo else 'CDS FASTA'
        if name.endswith('.bed') or name.endswith('.bed.gz'):
            return 'Gene Bed'
        if name.endswith('.genomic.fasta.masked') or name.endswith('.genomic.fasta.masked.gz'):
            return 'Genome FASTA'
        if name.endswith('.genomic.fasta') or name.endswith('.genomic.fasta.gz') or name.endswith('.fa') or name.endswith('.fa.gz') or name.endswith('.fasta') or name.endswith('.fasta.gz'):
            return 'Genome FASTA'
        return filename
    if 'emapper.annotations' in name or name.endswith('.annotations'):
        return 'EGGNOG Annotation'
    if name.endswith('.gff3'):
        return 'Gene GFF'
    if name.endswith('.gtf'):
        return 'Gene GTF'
    if name.endswith('.bed'):
        return 'Gene Bed'
    if '.pep.' in name or name.endswith('.pep.fasta') or name.endswith('.pep.faa') or name.endswith('.pep.faa.gz'):
        return 'Protein FASTA'
    if '.cds.' in name or name.endswith('.cds.fasta'):
        return 'CDS FASTA'
    if name.endswith('.fasta') or name.endswith('.fa'):
        return 'Genome FASTA'
    if name.endswith('.tsv'):
        return 'InterProScan Annotation'
    return filename


def _download_file_type_order(data_type):
    if data_type == 'genome':
        return [
            'Genome FASTA',
            'Gene GFF3',
            'Gene GFF3 (ANNEVO)',
            'Gene GTF',
            'Gene GTF (ANNEVO)',
            'Protein FASTA',
            'Protein FASTA (ANNEVO)',
            'CDS FASTA',
            'CDS FASTA (ANNEVO)',
            'Gene Bed',
            'Repeat FAMILY',
            'Repeat LOCATION',
            'Repeat SUMMARY',
            'EGGNOG Annotation',
            'EGGNOG Annotation (ANNEVO)',
            'InterProScan Annotation',
            'InterProScan Annotation (ANNEVO)',
        ]
    if data_type == 'mitogenome':
        return [
            'CDS',
            'Gene',
            'Protein',
            'Genome',
            'rRNA',
            'tRNA',
            'Location',
        ]
    if data_type == 'symbiont':
        return [
            'Genome',
            'Protein',
            'Gene',
            'GFF3',
        ]
    return []


def _get_download_file_options(data_type, identifier):
    file_entries = []
    if data_type == 'genome':
        tax = _resolve_taxonomy_identifier(identifier)
        file_entries = _collect_genome_files(tax)
    elif data_type == 'mitogenome':
        tax = _resolve_taxonomy_identifier(identifier)
        file_entries = _collect_mitogenome_files(tax)
    elif data_type == 'symbiont':
        file_entries = _collect_symbiont_files(identifier)
    elif data_type == 'singlecell':
        try:
            dataset = SingleCell.objects.get(pk=int(identifier))
        except (TypeError, ValueError, SingleCell.DoesNotExist):
            dataset = None
        if dataset:
            file_entries = _collect_singlecell_rds_files(dataset)
    result = []
    for _, rel_path in file_entries:
        filename = os.path.basename(rel_path)
        result.append({
            'path': rel_path,
            'filename': filename,
            'label': _download_file_label(filename, data_type),
        })
    order = _download_file_type_order(data_type)
    result.sort(key=lambda entry: order.index(entry['label']) if entry['label'] in order else len(order))
    return result


def _get_symbiont_file_options(identifier, sub=None):
    result = []
    for _, rel_path in _collect_symbiont_files(identifier, sub):
        filename = os.path.basename(rel_path)
        result.append({
            'path': rel_path,
            'filename': filename,
            'label': _download_file_label(filename, 'symbiont'),
        })
    order = _download_file_type_order('symbiont')
    result.sort(key=lambda entry: order.index(entry['label']) if entry['label'] in order else len(order))
    return result


def download_file_options(request):
    data_type = request.GET.get('data_type', '')
    identifier = request.GET.get('id', '')
    if not data_type or not identifier:
        return JsonResponse({'error': 'Missing data_type or id.'}, status=400)
    sub = request.GET.get('sub') if data_type == 'symbiont' else None
    files = _get_symbiont_file_options(identifier, sub) if data_type == 'symbiont' else _get_download_file_options(data_type, identifier)
    types = sorted(set([f['label'] for f in files]))
    return JsonResponse({'files': files, 'types': types})


def _allowed_download_file_types(data_type):
    if data_type == 'mitogenome':
        return {
            'CDS', 'Gene', 'Protein', 'Genome', 'rRNA', 'tRNA', 'Location'
        }
    if data_type == 'symbiont':
        return {
            'Genome', 'Protein', 'Gene', 'GFF3'
        }
    return {
        'Genome FASTA', 'Gene GFF3', 'Gene GFF3 (ANNEVO)',
        'Gene GTF', 'Gene GTF (ANNEVO)',
        'Protein FASTA', 'Protein FASTA (ANNEVO)',
        'CDS FASTA', 'CDS FASTA (ANNEVO)',
        'Gene Bed', 'Repeat FAMILY', 'Repeat LOCATION', 'Repeat SUMMARY',
        'EGGNOG Annotation', 'EGGNOG Annotation (ANNEVO)',
        'InterProScan Annotation', 'InterProScan Annotation (ANNEVO)'
    }


def download_file_types(request):
    data_type = request.GET.get('data_type', '')
    identifiers = request.GET.getlist('ids')
    if not data_type:
        return JsonResponse({'error': 'Missing data_type.'}, status=400)
    if data_type not in ('genome', 'mitogenome', 'symbiont'):
        return JsonResponse({'error': 'Unsupported data_type.'}, status=400)
    allowed_types = _allowed_download_file_types(data_type)
    types = set()
    for identifier in identifiers:
        files = _get_download_file_options(data_type, identifier)
        for file in files:
            label = file.get('label')
            if label and label in allowed_types:
                types.add(label)
    return JsonResponse({'types': sorted(types)})


def _parse_selected_file_item(item):
    if item is None:
        return None, None
    item = str(item)
    if '|' in item:
        identifier, file_type = item.split('|', 1)
        return identifier, file_type
    return item, None


def _make_zip_response(data_type, selected_items):
    try:
        archive_filename = f"AnnelidaDB_{data_type}_download.zip"
        selection_limits = {
            'genome': 12,
            'mitogenome': 20,
            'symbiont': 20,
            'singlecell': 40,
            'transcriptome': 40,
            'taxonomy': 40,
        }
        if data_type in selection_limits and len(selected_items) > selection_limits[data_type]:
            return None, f"Too many items selected for {data_type} download. Please select up to {selection_limits[data_type]} items and try again."

        if data_type == 'genome':
            selected_files = []
            seen_paths = set()
            for item in selected_items:
                identifier, file_type = _parse_selected_file_item(item)
                tax = _resolve_taxonomy_identifier(identifier)
                files = _collect_genome_files(tax)
                if file_type:
                    files = [f for f in files if _download_file_label(os.path.basename(f[1]), 'genome') == file_type]
                for abs_path, arcname in files:
                    if abs_path not in seen_paths:
                        seen_paths.add(abs_path)
                        selected_files.append((abs_path, arcname))
        elif data_type == 'mitogenome':
            selected_files = []
            seen_paths = set()
            for item in selected_items:
                identifier, file_type = _parse_selected_file_item(item)
                tax = _resolve_taxonomy_identifier(identifier)
                files = _collect_mitogenome_files(tax)
                if file_type:
                    files = [f for f in files if _download_file_label(os.path.basename(f[1]), 'mitogenome') == file_type]
                for abs_path, arcname in files:
                    if abs_path not in seen_paths:
                        seen_paths.add(abs_path)
                        selected_files.append((abs_path, arcname))
        elif data_type == 'symbiont':
            selected_files = []
            seen_paths = set()
            for item in selected_items:
                identifier, file_type = _parse_selected_file_item(item)
                files = _collect_symbiont_files(identifier)
                if file_type:
                    files = [f for f in files if _download_file_label(os.path.basename(f[1]), 'symbiont') == file_type]
                for abs_path, arcname in files:
                    if abs_path not in seen_paths:
                        seen_paths.add(abs_path)
                        selected_files.append((abs_path, arcname))
        elif data_type == 'taxonomy':
            rows = Taxonomy.objects.filter(ncbi_taxonomy_id__in=selected_items)
            if not rows.exists():
                return None, 'No selected taxonomy entries found.'
            csv_buffer = io.StringIO()
            writer = csv.writer(csv_buffer)
            writer.writerow(['Species', 'NCBI Taxonomy ID', 'Class', 'Order', 'Family', 'Genome Available', 'Mito Available', 'DNA Available', 'RNA Available'])
            for row in rows:
                ncbi_id = row.ncbi_taxonomy_id
                writer.writerow([
                    row.species,
                    ncbi_id,
                    row.class_name or '',
                    row.order or '',
                    row.family or '',
                    'Yes' if row.genomes.exists() else 'No',
                    'Yes' if row.mitogenomes.exists() else 'No',
                    'Yes' if _has_ncbi_dna(ncbi_id, bool(getattr(row, 'dna_link', None))) else 'No',
                    'Yes' if _has_ncbi_rna(ncbi_id, bool(getattr(row, 'rna_link', None))) else 'No',
                ])
            selected_files = [('taxonomy_metadata.csv', csv_buffer.getvalue())]
        elif data_type == 'singlecell':
            rows = SingleCell.objects.select_related('taxonomy').filter(pk__in=[int(i) for i in selected_items if str(i).isdigit()])
            if not rows.exists():
                return None, 'No selected single-cell datasets found.'
            csv_buffer = io.StringIO()
            writer = csv.writer(csv_buffer)
            writer.writerow(['Species', 'Directory', 'Tissue', 'Stage', 'Clusters', 'Cells', 'Type', 'Year', 'Platform', 'BioProject', 'Journal', 'DOI', 'Title'])
            for row in rows:
                writer.writerow([
                    row.species,
                    row.dir or '',
                    row.tissue or '',
                    row.stage or '',
                    row.clusters or '',
                    row.cells or '',
                    row.type or '',
                    row.year or '',
                    row.platform or '',
                    row.bioproject or '',
                    row.journal or '',
                    row.doi or '',
                    row.title or '',
                ])
            selected_files = [('singlecell_metadata.csv', csv_buffer.getvalue())]
            for row in rows:
                selected_files.extend(_collect_singlecell_rds_files(row))
        elif data_type == 'transcriptome':
            rows = _get_sample_transcriptome_rows()
            selected_indexes = [int(i) for i in selected_items if str(i).isdigit()]
            selected_rows = [rows[i] for i in selected_indexes if 0 <= i < len(rows)]
            if not selected_rows:
                return None, 'No selected transcriptome datasets found.'
            csv_buffer = io.StringIO()
            writer = csv.writer(csv_buffer)
            writer.writerow(['Species', 'Tissue', 'Library Type', 'Reads', 'Assembly Size', 'Complete Genes', 'Status'])
            for row in selected_rows:
                writer.writerow([
                    row['species'],
                    row['tissue'],
                    row['library_type'],
                    row['reads'],
                    row['assembly_size'],
                    row['complete_genes'],
                    row['status'],
                ])
            selected_files = [('transcriptome_metadata.csv', csv_buffer.getvalue())]
        elif data_type == 'selected_files':
            selected_files = []
            allowed_base = os.path.normpath(os.path.join(settings.BASE_DIR, 'static', 'data'))
            for item in selected_items:
                if not item:
                    continue
                normalized = os.path.normpath(os.path.join(settings.BASE_DIR, item.lstrip('/\\')))
                if not normalized.startswith(allowed_base):
                    return None, 'Invalid selected file path.'
                if not os.path.isfile(normalized):
                    return None, f'File not found: {item}'
                archive_name = os.path.relpath(normalized, allowed_base)
                selected_files.append((normalized, archive_name))
        else:
            return None, 'Unsupported download type.'

        if data_type == 'singlecell':
            temp_file = tempfile.SpooledTemporaryFile(max_size=20 * 1024 * 1024)
            with zipfile.ZipFile(temp_file, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                archive.writestr(selected_files[0][0], selected_files[0][1])
                for abs_path, arcname in selected_files[1:]:
                    archive.write(abs_path, arcname)
            temp_file.seek(0)
            response = FileResponse(temp_file, as_attachment=True, filename=archive_filename)
            return response, None

        if data_type == 'transcriptome':
            temp_file = tempfile.SpooledTemporaryFile(max_size=20 * 1024 * 1024)
            with zipfile.ZipFile(temp_file, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                for filename, content in selected_files:
                    archive.writestr(filename, content)
            temp_file.seek(0)
            response = FileResponse(temp_file, as_attachment=True, filename=archive_filename)
            return response, None

        if not selected_files:
            return None, 'No downloadable files found for the selected items.'

        temp_file = tempfile.SpooledTemporaryFile(max_size=50 * 1024 * 1024)
        with zipfile.ZipFile(temp_file, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for abs_path, arcname in selected_files:
                archive.write(abs_path, arcname)
        temp_file.seek(0)
        response = FileResponse(temp_file, as_attachment=True, filename=archive_filename)
        return response, None
    except Exception as exc:
        traceback.print_exc()
        return None, f"Error creating download archive: {str(exc)}"


def download_page(request):
    return render(request, 'download.html', {
        'page_title': 'Download',
        'active_page': 'download',
        'genome_rows': Genome.objects.select_related('taxonomy').all(),
        'mitogenome_rows': Mitogenome.objects.select_related('taxonomy').all(),
        'symbiont_rows': _build_symbiont_host_rows(),
        'singlecell_rows': SingleCell.objects.select_related('taxonomy').all(),
    })


def download_batch(request):
    if request.method != 'POST':
        return redirect('download_page')

    data_type = request.POST.get('data_type', 'genome')
    selected_items = request.POST.getlist('selected')
    try:
        response, error = _make_zip_response(data_type, selected_items)
    except Exception as exc:
        traceback.print_exc()
        response, error = None, f"Download failed due to an internal error: {str(exc)}"

    if response:
        return response

    return render(request, 'download.html', {
        'page_title': 'Download',
        'active_page': 'download',
        'download_error': error,
        'genome_rows': Genome.objects.select_related('taxonomy').all(),
        'mitogenome_rows': Mitogenome.objects.select_related('taxonomy').all(),
        'symbiont_rows': _build_symbiont_host_rows(),
        'singlecell_rows': SingleCell.objects.select_related('taxonomy').all(),
    })


def _build_symbiont_host_rows():
    host_map = {}
    for sym in Symbiont.objects.order_by('host_species', 'sub'):
        host = sym.host_species or 'Unknown host'
        sub = sym.sub or 'symbiont'
        row_key = (host, sub)
        if row_key not in host_map:
            tax = get_taxonomy_for_species_name(host)
            host_id = tax.ncbi_taxonomy_id if tax and tax.ncbi_taxonomy_id else (tax.species if tax else '')
            host_map[row_key] = {
                'host_species': host,
                'sub': sub,
                'taxonomy_id': host_id,
                'detail_url_sub': sub,
                'symbiont_type': sym.symbiont_type or 'Unknown',
                'class_name': tax.class_name if tax else 'Unknown',
                'order': tax.order if tax else 'Unknown',
                'family': tax.family if tax else 'Unknown',
            }
    return list(host_map.values())


# ──────────────────────────────────────────────
# Helper: format species name for file paths
# ──────────────────────────────────────────────

def _format_species(raw_species):
    if not raw_species:
        return ''
    formatted = raw_species.replace('sp.', 'sp_')
    formatted = formatted.replace(' ', '_')
    formatted = re.sub(r'__+', '_', formatted)
    return formatted


def _species_data_dir(raw_species):
    formatted_species = _format_species(raw_species)
    if not formatted_species:
        return None
    return os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species)


def _species_genome_dir(raw_species):
    data_dir = _species_data_dir(raw_species)
    if not data_dir:
        return None
    return os.path.join(data_dir, 'genome')


def _species_genome_url(raw_species):
    formatted_species = _format_species(raw_species)
    if not formatted_species:
        return None
    return f"/static/data/{formatted_species}/genome"


def _singlecell_candidate_dirs(formatted_species, dataset):
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, 'single-cell')
    dataset_dir = (getattr(dataset, 'dir', None) or '').strip()
    if dataset_dir:
        normalized_dir = os.path.normpath(dataset_dir)
        if normalized_dir not in ('', '.') and normalized_dir != '..' and not normalized_dir.startswith('../') and not os.path.isabs(normalized_dir):
            if normalized_dir == 'single-cell':
                return [base_dir]
            if normalized_dir.startswith('single-cell' + os.sep):
                normalized_dir = normalized_dir[len('single-cell' + os.sep):]
            return [os.path.join(base_dir, normalized_dir)]
        return []

    tissue = (dataset.tissue or '').strip()
    candidate_dirs = [base_dir]
    if tissue:
        candidate_dirs.insert(0, os.path.join(base_dir, tissue.replace(' ', '_')))
        if tissue.replace(' ', '_') != tissue:
            candidate_dirs.insert(0, os.path.join(base_dir, tissue))
    return candidate_dirs


def _find_singlecell_umap_files(formatted_species, singlecell_qs):
    umap_files = []
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, 'single-cell')
    if not os.path.isdir(base_dir):
        return umap_files

    def find_in_dir(search_dir):
        if not os.path.isdir(search_dir):
            return None
        for filename in os.listdir(search_dir):
            if filename.lower().endswith('.csv') and 'umap' in filename.lower():
                return os.path.join(search_dir, filename)
        return None

    def find_any_umap():
        for root, _, files in os.walk(base_dir):
            for filename in files:
                if filename.lower().endswith('.csv') and 'umap' in filename.lower():
                    return os.path.join(root, filename)
        return None

    for dataset in singlecell_qs:
        candidate_dirs = _singlecell_candidate_dirs(formatted_species, dataset)

        umap_path = None
        for search_dir in candidate_dirs:
            umap_path = find_in_dir(search_dir)
            if umap_path:
                break
        if not umap_path and not getattr(dataset, 'dir', None):
            umap_path = find_any_umap()

        if umap_path:
            umap_files.append({
                'dataset_id': dataset.id,
                'tissue': dataset.tissue,
                'stage': dataset.stage,
                'year': dataset.year,
                'url': _format_static_url(umap_path),
                'filename': os.path.basename(umap_path),
                'label': ' '.join(filter(None, [dataset.tissue, dataset.stage, dataset.year])) or 'Single-cell UMAP',
            })

    return umap_files


def _find_singlecell_markers_file(formatted_species, dataset):
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, 'single-cell')
    candidate_dirs = _singlecell_candidate_dirs(formatted_species, dataset)

    for search_dir in candidate_dirs:
        candidate_path = os.path.join(search_dir, 'Markers.csv')
        if os.path.isfile(candidate_path):
            return candidate_path

    if getattr(dataset, 'dir', None):
        return None

    for root, _, files in os.walk(base_dir):
        for filename in files:
            if filename.lower() == 'markers.csv':
                return os.path.join(root, filename)
    return None


def _find_singlecell_expression_file(formatted_species, dataset):
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, 'single-cell')
    candidate_dirs = _singlecell_candidate_dirs(formatted_species, dataset)

    for search_dir in candidate_dirs:
        for filename in ('Expression.csv', 'Expression_with_meta.csv'):
            candidate_path = os.path.join(search_dir, filename)
            if os.path.isfile(candidate_path):
                return candidate_path

    if getattr(dataset, 'dir', None):
        return None

    for root, _, files in os.walk(base_dir):
        for filename in files:
            if filename.lower() in {'expression.csv', 'expression_with_meta.csv'}:
                return os.path.join(root, filename)
    return None


def _find_singlecell_umap_file(formatted_species, dataset):
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, 'single-cell')
    candidate_dirs = _singlecell_candidate_dirs(formatted_species, dataset)

    for search_dir in candidate_dirs:
        if not os.path.isdir(search_dir):
            continue
        for filename in os.listdir(search_dir):
            if filename.lower().endswith('.csv') and 'umap' in filename.lower():
                return os.path.join(search_dir, filename)

    if getattr(dataset, 'dir', None):
        return None

    for root, _, files in os.walk(base_dir):
        for filename in files:
            if filename.lower().endswith('.csv') and 'umap' in filename.lower():
                return os.path.join(root, filename)
    return None


def _find_singlecell_cluster_image_file(formatted_species, dataset):
    for search_dir in _singlecell_candidate_dirs(formatted_species, dataset):
        if not os.path.isdir(search_dir):
            continue
        filenames = sorted(os.listdir(search_dir))
        for filename in filenames:
            lower_name = filename.lower()
            if lower_name.endswith(('.jpg', '.jpeg', '.png', '.webp')) and ('tsne' in lower_name or 'umap' in lower_name):
                return os.path.join(search_dir, filename)
    return None


def _generate_singlecell_expression_files(expression_csv_path, umap_csv_path=None):
    script_path = os.path.join(settings.BASE_DIR, 'scripts', 'generate_expression_images.py')
    if not os.path.isfile(script_path):
        print(f'Expression generation script not found: {script_path}')
        return

    dataset_dir = os.path.dirname(expression_csv_path)
    try:
        subprocess.run(
            [sys.executable, script_path, dataset_dir, '--clean', '--overwrite', '--quiet'],
            cwd=settings.BASE_DIR,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        print('Failed to generate single-cell expression files:')
        print(exc.stdout)
        print(exc.stderr)


def _enqueue_singlecell_expression_generation(expression_csv_path, umap_csv_path=None):
    _enqueue_job(_generate_singlecell_expression_files, expression_csv_path, umap_csv_path)


def _find_taxonomy_for_jbrowse_data(data_dir):
    if not data_dir:
        return None
    data_root = os.path.normpath(os.path.join(settings.BASE_DIR, 'static', 'data'))
    try:
        rel_path = os.path.relpath(data_dir, data_root)
    except ValueError:
        return None
    species_dir = rel_path.split(os.sep)[0] if rel_path else None
    if not species_dir:
        return None
    species_name = species_dir.replace('_', ' ')
    return get_taxonomy_for_species_name(species_name)


# ──────────────────────────────────────────────
# Helper: ALG & Oxford plot URLs
# ──────────────────────────────────────────────

def _get_alg_oxford_urls(genome, formatted_species):
    alg_url = None
    oxford_url = None
    if genome and genome.assembly_level and 'chromosome' in genome.assembly_level.lower():
        alg_dir = os.path.join(settings.BASE_DIR, 'static', 'plots', 'alg')
        alg_prefix = f"Annelida.{formatted_species}.ALGs.groups"
        if os.path.exists(alg_dir):
            for file in os.listdir(alg_dir):
                if file.startswith(alg_prefix) and 'karyotype' in file:
                    alg_url = f"/static/plots/alg/{file}"
                    break
        oxford_filename = f"Annelida.{formatted_species}.oxford_plot.png"
        oxford_path = os.path.join(settings.BASE_DIR, 'static', 'plots', 'alg', oxford_filename)
        if os.path.exists(oxford_path):
            oxford_url = f"/static/plots/alg/{oxford_filename}"
        else:
            oxford_filename = f"Annelida.{formatted_species}.oxford_plot.pdf"
            oxford_path = os.path.join(settings.BASE_DIR, 'static', 'plots', 'alg', oxford_filename)
            if os.path.exists(oxford_path):
                oxford_url = f"/static/plots/alg/{oxford_filename}"
    return alg_url, oxford_url


# ──────────────────────────────────────────────
# Species Detail (Hub)
# ──────────────────────────────────────────────

def _fetch_gbif_species_data(gbif_id):
    if not gbif_id:
        return None, None

    cache_key = f"gbif_species:{gbif_id}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('data'), cached.get('error')

    url = f"https://api.gbif.org/v1/species/{gbif_id}"
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            data = json.load(response)
    except urllib.error.HTTPError as exc:
        error = f"GBIF HTTP {exc.code}"
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except urllib.error.URLError as exc:
        reason = exc.reason if hasattr(exc, 'reason') else str(exc)
        error = f"GBIF network error: {reason}"
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except Exception as exc:
        error = f"GBIF fetch error: {str(exc)}"
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error

    gbif_data = {
        'key': data.get('key'),
        'scientific_name': data.get('scientificName'),
        'canonical_name': data.get('canonicalName'),
        'authorship': data.get('authorship'),
        'rank': data.get('rank'),
        'status': data.get('taxonomicStatus'),
        'kingdom': data.get('kingdom'),
        'phylum': data.get('phylum'),
        'class_name': data.get('class'),
        'order': data.get('order'),
        'family': data.get('family'),
        'genus': data.get('genus'),
        'species': data.get('species'),
        'nub_key': data.get('nubKey'),
        'num_descendants': data.get('numDescendants'),
        'origin': data.get('origin'),
        'gbif_url': f"https://www.gbif.org/species/{gbif_id}",
    }
    cache.set(cache_key, {'data': gbif_data, 'error': None}, 3600)
    return gbif_data, None


def _flatten_aphia_classification(node):
    classification = []
    current = node
    while current:
        classification.append({
            'aphia_id': current.get('AphiaID'),
            'rank': current.get('rank'),
            'scientific_name': current.get('scientificname'),
        })
        current = current.get('child')
    return classification


def _fetch_worms_taxonomy(worms_id):
    if not worms_id:
        return None, None

    cache_key = f"worms_taxonomy:{worms_id}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('data'), cached.get('error')

    record_url = f"https://www.marinespecies.org/rest/AphiaRecordByAphiaID/{worms_id}"
    classification_url = f"https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/{worms_id}"
    worms_data = None
    error = None

    try:
        with urllib.request.urlopen(record_url, timeout=10) as response:
            record = json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            error = 'WoRMS record not found.'
        else:
            error = f'WoRMS HTTP {exc.code}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except urllib.error.URLError as exc:
        reason = exc.reason if hasattr(exc, 'reason') else str(exc)
        error = f'WoRMS network error: {reason}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except Exception as exc:
        error = f'WoRMS fetch error: {str(exc)}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error

    classification = []
    try:
        with urllib.request.urlopen(classification_url, timeout=10) as response:
            classification_payload = json.load(response)
            if classification_payload:
                classification = _flatten_aphia_classification(classification_payload)
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            error = f'WoRMS classification HTTP {exc.code}'
    except urllib.error.URLError as exc:
        reason = exc.reason if hasattr(exc, 'reason') else str(exc)
        error = f'WoRMS classification network error: {reason}'
    except Exception as exc:
        error = f'WoRMS classification fetch error: {str(exc)}'

    environment = []
    if record.get('isMarine'):
        environment.append('Marine')
    if record.get('isBrackish'):
        environment.append('Brackish')
    if record.get('isFreshwater'):
        environment.append('Freshwater')
    if record.get('isTerrestrial'):
        environment.append('Terrestrial')
    if record.get('isExtinct'):
        environment.append('Extinct')

    worms_data = {
        'aphia_id': record.get('AphiaID'),
        'scientific_name': record.get('scientificname'),
        'valid_name': record.get('valid_name'),
        'rank': record.get('rank'),
        'status': record.get('status'),
        'kingdom': record.get('kingdom'),
        'phylum': record.get('phylum'),
        'class_name': record.get('class'),
        'order': record.get('order'),
        'family': record.get('family'),
        'genus': record.get('genus'),
        'environment': environment,
        'classification': classification,
        'worms_url': record.get('url') or f'https://www.marinespecies.org/aphia.php?p=taxdetails&id={worms_id}',
    }
    cache.set(cache_key, {'data': worms_data, 'error': None}, 3600)
    return worms_data, error


def _fetch_ncbi_taxonomy(ncbi_id):
    if not ncbi_id:
        return None, None

    cache_key = f"ncbi_taxonomy:{ncbi_id}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('data'), cached.get('error')

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ncbi_id}&retmode=xml"
    try:
        with urllib.request.urlopen(url, timeout=15) as response:
            raw_xml = response.read().decode('utf-8')
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            error = 'NCBI taxonomy record not found.'
        else:
            error = f'NCBI HTTP {exc.code}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except urllib.error.URLError as exc:
        reason = exc.reason if hasattr(exc, 'reason') else str(exc)
        error = f'NCBI network error: {reason}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except Exception as exc:
        error = f'NCBI fetch error: {str(exc)}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error

    try:
        tree = ET.fromstring(raw_xml)
        taxon = tree.find('.//Taxon')
        if taxon is None:
            raise ValueError('NCBI taxonomy XML missing Taxon entry')

        classification = []
        for lineage_taxon in taxon.findall('.//LineageEx/Taxon'):
            classification.append({
                'taxid': lineage_taxon.findtext('TaxId'),
                'scientific_name': lineage_taxon.findtext('ScientificName'),
                'rank': lineage_taxon.findtext('Rank'),
            })

        ncbi_data = {
            'taxid': taxon.findtext('TaxId'),
            'scientific_name': taxon.findtext('ScientificName'),
            'rank': taxon.findtext('Rank'),
            'division': taxon.findtext('Division'),
            'lineage': taxon.findtext('Lineage'),
            'classification': classification,
            'ncbi_url': f'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={ncbi_id}',
        }
        cache.set(cache_key, {'data': ncbi_data, 'error': None}, 3600)
        return ncbi_data, None
    except Exception as exc:
        error = f'NCBI parse error: {str(exc)}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error


def _ncbi_nucleotide_search_count(ncbi_id, term):
    if not ncbi_id or not term:
        return 0
    cache_key = f"ncbi_nucleotide_count:{ncbi_id}:{hashlib.md5(term.encode('utf-8')).hexdigest()}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {
        'db': 'nuccore',
        'term': term,
        'retmode': 'json',
    }
    query_url = f"{url}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(query_url, timeout=15) as response:
            raw = response.read().decode('utf-8')
            data = json.loads(raw)
            count = int(data.get('esearchresult', {}).get('count', 0))
    except Exception:
        count = 0

    cache.set(cache_key, count, 3600)
    return count


def _ncbi_dna_query(ncbi_id):
    return f"txid{ncbi_id}[Organism:exp] AND (genomic[LibrarySource] OR Genomic[LibrarySource] OR DNA[LibrarySource] OR genomic[LibraryStrategy] OR WGS[LibraryStrategy])"


def _ncbi_rna_query(ncbi_id):
    return f"txid{ncbi_id}[Organism:exp] AND (rna[LibrarySource] OR transcriptomic[LibrarySource] OR RNA-Seq[LibraryStrategy] OR TRANSCRIPTOMIC[LibraryStrategy])"


def _ncbi_sra_dna_query(ncbi_id):
    return (
        f"txid{ncbi_id}[Organism:exp] AND (genomic[LibrarySource] OR DNA[LibrarySource] "
        f"OR genomic[LibraryStrategy] OR WGS[LibraryStrategy] OR biomol_dna[prop])"
    )


def _ncbi_sra_rna_query(ncbi_id):
    return (
        f"txid{ncbi_id}[Organism:exp] AND (rna[LibrarySource] OR transcriptomic[LibrarySource] "
        f"OR RNA-Seq[LibraryStrategy] OR RNA[LibraryStrategy] OR biomol_rna[prop])"
    )


def _ncbi_sra_ids(term, api_key=None):
    if not term:
        return []
    retmax = 500
    all_ids = []
    start = 0
    while True:
        params = {
            'db': 'sra',
            'term': term,
            'retmode': 'json',
            'retmax': str(retmax),
            'retstart': str(start),
        }
        api_key = api_key or os.environ.get('NCBI_API_KEY')
        if api_key:
            params['api_key'] = api_key
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{urllib.parse.urlencode(params)}"
        request = urllib.request.Request(url, headers={'User-Agent': 'AnnelidaDB/1.0'})
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                result = json.loads(response.read().decode('utf-8')).get('esearchresult', {})
        except (OSError, ValueError, TypeError):
            return all_ids
        page_ids = result.get('idlist', [])
        count = int(result.get('count', 0) or 0)
        if not page_ids:
            break
        all_ids.extend(page_ids)
        start += len(page_ids)
        if start >= count:
            break
        time.sleep(_ncbi_request_delay())
    return all_ids


def _ncbi_sra_run_accessions(term, api_key=None):
    ids = _ncbi_sra_ids(term, api_key=api_key)
    if not ids:
        return []
    accessions = []
    seen = set()
    chunk_size = 200
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i + chunk_size]
        params = {
            'db': 'sra',
            'id': ','.join(chunk),
            'rettype': 'runinfo',
            'retmode': 'text',
        }
        api_key = api_key or os.environ.get('NCBI_API_KEY')
        if api_key:
            params['api_key'] = api_key
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urllib.parse.urlencode(params)}"
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                text = response.read().decode('utf-8', errors='replace')
                rows = list(csv.DictReader(line for line in text.splitlines() if line.strip()))
                for row in rows:
                    run = (row.get('Run') or '').strip()
                    if run and run not in seen:
                        seen.add(run)
                        accessions.append(run)
        except Exception:
            continue
        if i + chunk_size < len(ids):
            time.sleep(_ncbi_request_delay())
    return sorted(accessions)


def _ncbi_sra_runinfo(term, api_key=None, chunk_size=200):
    ids = _ncbi_sra_ids(term, api_key=api_key)
    if not ids:
        return []

    rows = []
    api_key = api_key or os.environ.get('NCBI_API_KEY')
    for start in range(0, len(ids), chunk_size):
        chunk = ids[start:start + chunk_size]
        params = {
            'db': 'sra',
            'id': ','.join(chunk),
            'rettype': 'runinfo',
            'retmode': 'text',
        }
        if api_key:
            params['api_key'] = api_key
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urllib.parse.urlencode(params)}"
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                text = response.read().decode('utf-8', errors='replace')
            rows.extend(csv.DictReader(line for line in text.splitlines() if line.strip()))
        except (OSError, csv.Error):
            continue
        if start + chunk_size < len(ids):
            time.sleep(_ncbi_request_delay())
    return rows


def _ncbi_sra_runinfo_from_webenv(query_key, webenv, api_key=None, retmax=200, max_results=0):
    if not query_key or not webenv:
        return []
    rows = []
    start = 0
    while True:
        params = {
            'db': 'sra',
            'query_key': query_key,
            'WebEnv': webenv,
            'rettype': 'runinfo',
            'retmode': 'text',
            'retstart': str(start),
            'retmax': str(retmax),
        }
        api_key = api_key or os.environ.get('NCBI_API_KEY')
        if api_key:
            params['api_key'] = api_key
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urllib.parse.urlencode(params)}"
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                text = response.read().decode('utf-8', errors='replace')
        except Exception:
            break
        chunk_rows = list(csv.DictReader(line for line in text.splitlines() if line.strip()))
        if not chunk_rows:
            break
        rows.extend(chunk_rows)
        if max_results and len(rows) >= max_results:
            rows = rows[:max_results]
            break
        if len(chunk_rows) < retmax:
            break
        start += len(chunk_rows)
    return rows


def _ncbi_nucleotide_search_url(ncbi_id, term):
    if not ncbi_id or not term:
        return None
    base = 'https://www.ncbi.nlm.nih.gov/nucleotide/?term='
    return base + urllib.parse.quote_plus(term)


def _ncbi_sra_search_url(ncbi_id, term):
    if not ncbi_id or not term:
        return None
    base = 'https://www.ncbi.nlm.nih.gov/sra/?term='
    return base + urllib.parse.quote_plus(term)


def _parse_obis_depth_value(record):
    for key in ('depth', 'decimalDepth', 'verbatimDepth', 'bathymetry'):
        if key not in record or record[key] is None:
            continue
        try:
            value = float(record[key])
        except (TypeError, ValueError):
            continue
        if key == 'bathymetry':
            value = abs(value)
        return value
    return None


def _fetch_obis_depth_summary(worms_id, scientific_name=None, sample_size=200):
    if not worms_id or not scientific_name:
        return None, 'Missing WoRMS ID or scientific name for OBIS query.'

    cache_key = _safe_cache_key('obis_depth', f'{worms_id}:{scientific_name}')
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('data'), cached.get('error')

    query_name = urllib.parse.quote(scientific_name)
    url = f"https://api.obis.org/v3/occurrence?aphiaID={worms_id}&scientificname={query_name}&fields=depth,decimalDepth,verbatimDepth,bathymetry&size={sample_size}"

    try:
        with urllib.request.urlopen(url, timeout=15) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as exc:
        error = f'OBIS HTTP {exc.code}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except urllib.error.URLError as exc:
        reason = exc.reason if hasattr(exc, 'reason') else str(exc)
        error = f'OBIS network error: {reason}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except Exception as exc:
        error = f'OBIS fetch error: {str(exc)}'
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error

    depths = []
    for record in payload.get('results', []):
        value = _parse_obis_depth_value(record)
        if value is not None and value >= 0:
            depths.append(value)

    if not depths:
        summary = {'count': 0, 'min_depth': None, 'max_depth': None, 'mean_depth': None}
    else:
        depths.sort()
        summary = {
            'count': len(depths),
            'min_depth': depths[0],
            'max_depth': depths[-1],
            'mean_depth': sum(depths) / len(depths),
        }
    cache.set(cache_key, {'data': summary, 'error': None}, 3600)
    return summary, None


def _fetch_gbif_species_data(gbif_id):
    if not gbif_id:
        return None, None

    cache_key = f"gbif_species:{gbif_id}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached.get('data'), cached.get('error')

    url = f"https://api.gbif.org/v1/species/{gbif_id}"
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            data = json.load(response)
    except urllib.error.HTTPError as exc:
        error = f"GBIF HTTP {exc.code}"
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except urllib.error.URLError as exc:
        reason = exc.reason if hasattr(exc, 'reason') else str(exc)
        error = f"GBIF network error: {reason}"
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error
    except Exception as exc:
        error = f"GBIF fetch error: {str(exc)}"
        cache.set(cache_key, {'data': None, 'error': error}, 3600)
        return None, error

    gbif_data = {
        'key': data.get('key'),
        'scientific_name': data.get('scientificName'),
        'canonical_name': data.get('canonicalName'),
        'authorship': data.get('authorship'),
        'rank': data.get('rank'),
        'status': data.get('taxonomicStatus'),
        'kingdom': data.get('kingdom'),
        'phylum': data.get('phylum'),
        'class_name': data.get('class'),
        'order': data.get('order'),
        'family': data.get('family'),
        'genus': data.get('genus'),
        'species': data.get('species'),
        'nub_key': data.get('nubKey'),
        'num_descendants': data.get('numDescendants'),
        'origin': data.get('origin'),
        'gbif_url': f"https://www.gbif.org/species/{gbif_id}",
    }
    cache.set(cache_key, {'data': gbif_data, 'error': None}, 3600)
    return gbif_data, None


def _fetch_gbif_images(gbif_id, media_limit=12):
    images = []
    if not gbif_id:
        return images, None

    base = 'https://api.gbif.org/v1/occurrence/search'
    try:
        q = f"{base}?taxonKey={gbif_id}&mediaType=StillImage&limit={media_limit}"
        with urllib.request.urlopen(q, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            for rec in data.get('results', []):
                media = rec.get('media') or []
                for m in media:
                    url = m.get('identifier') or m.get('references') or m.get('content')
                    if not url:
                        continue
                    images.append({
                        'url': url,
                        'creator': m.get('creator'),
                        'license': m.get('license'),
                        'type': m.get('type'),
                    })
                    if len(images) >= media_limit:
                        break
                if len(images) >= media_limit:
                    break
    except Exception as e:
        return [], str(e)

    return images, None


@require_GET
def gbif_images_json(request, gbif_id):
    cached = _load_gbif_cache('images', gbif_id)
    if cached is not None:
        return JsonResponse({'images': cached, 'cached': True})

    images, error = _fetch_gbif_images(gbif_id)
    if error:
        return JsonResponse({'error': error}, status=502)

    _save_gbif_cache('images', gbif_id, images)
    return JsonResponse({'images': images, 'cached': False})


@require_GET
def gbif_species_json(request, gbif_id):
    gbif_data, error = _fetch_gbif_species_data(gbif_id)
    if error:
        return JsonResponse({'error': error}, status=502)
    return JsonResponse({'data': gbif_data, 'cached': False})


@require_GET
def worms_taxonomy_json(request, worms_id):
    worms_data, error = _fetch_worms_taxonomy(worms_id)
    if error:
        return JsonResponse({'error': error}, status=502)
    return JsonResponse({'data': worms_data, 'cached': False})


@require_GET
def ncbi_taxonomy_json(request, ncbi_id):
    ncbi_data, error = _fetch_ncbi_taxonomy(ncbi_id)
    if error:
        return JsonResponse({'error': error}, status=502)
    return JsonResponse({'data': ncbi_data, 'cached': False})


@require_GET
def obis_depth_json(request, worms_id):
    scientific_name = request.GET.get('scientific_name', '').strip()
    obis_data, error = _fetch_obis_depth_summary(worms_id, scientific_name=scientific_name)
    if error:
        return JsonResponse({'error': error}, status=502)
    return JsonResponse({'data': obis_data, 'cached': False})


@require_GET
def gbif_occurrences_json(request, gbif_id):
    """Return cached GBIF occurrence coordinates for given GBIF taxon key.
    Uses paging to collect up to a configurable maximum and caches the result.
    """
    cached = _load_gbif_cache('points', gbif_id)
    if cached is not None:
        return JsonResponse({'points': cached, 'cached': True})

    base = 'https://api.gbif.org/v1/occurrence/search'
    points = []
    limit = 300
    offset = 0
    max_total = 1200
    try:
        while True:
            q = f"{base}?taxonKey={gbif_id}&hasCoordinate=true&limit={limit}&offset={offset}"
            with urllib.request.urlopen(q, timeout=20) as resp:
                data = json.loads(resp.read().decode('utf-8'))
            results = data.get('results', [])
            if not results:
                break
            for rec in results:
                lat = rec.get('decimalLatitude')
                lon = rec.get('decimalLongitude')
                if lat is None or lon is None:
                    continue
                try:
                    latf = float(lat)
                    lonf = float(lon)
                except Exception:
                    continue
                points.append([latf, lonf])
            offset += limit
            if len(points) >= max_total:
                break
            # if results length less than limit, we've reached the end
            if len(results) < limit:
                break
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=502)

    _save_gbif_cache('points', gbif_id, points)
    return JsonResponse({'points': points, 'cached': False})


@require_GET
def backbone_distribution_json(request, taxonomy_id):
    tax = get_taxonomy_or_404(taxonomy_id)
    cache_key = _safe_cache_key('backbone_distribution', str(tax.ncbi_taxonomy_id or tax.species))
    cached = cache.get(cache_key) is not None
    records, summary = _fetch_backbone_distribution(tax)
    return JsonResponse({
        'records': records,
        'summary': summary,
        'cached': cached,
    })


def species_detail_page(request, taxonomy_id):
    tax = get_taxonomy_or_404(taxonomy_id)
    genome = tax.genomes.first()
    mitogenome = tax.mitogenomes.first()
    formatted_species = _format_species(tax.species)

    external_epigenome_link = getattr(tax, 'epigenome_link', None)
    epigenome_csv_records = _get_epigenome_records_for_species(tax.species)
    has_epigenome_csv = bool(epigenome_csv_records)
    has_epigenome = bool(external_epigenome_link or has_epigenome_csv)

    external_proteome_link = getattr(tax, 'proteome_link', None)
    proteome_csv_records = _get_proteome_records_for_species(tax.species)
    has_proteome_csv = bool(proteome_csv_records)
    has_proteome = bool(external_proteome_link or has_proteome_csv)
    proteome_link = external_proteome_link

    has_transcriptome = bool(getattr(tax, 'rna_link', None) or _species_has_transcriptome_data(tax.species))

    busco_filename = f"{formatted_species}.blobtools.snail.png"
    busco_url = f"/static/plots/busco/{busco_filename}"

    alg_url, oxford_url = _get_alg_oxford_urls(genome, formatted_species)
    singlecell_qs = SingleCell.objects.filter(species__iexact=tax.species).select_related('taxonomy')
    umap_datasets = _find_singlecell_umap_files(formatted_species, singlecell_qs)

    gbif_id = getattr(tax, 'gbif_id', None)

    updated_clades = [
        tax.updated_clade6,
        tax.updated_clade5,
        tax.updated_clade4,
        tax.updated_clade3,
        tax.updated_clade2,
        tax.updated_clade1,
    ]
    updated_taxonomy_line = '; '.join([clade for clade in updated_clades if clade])

    context = {
        'tax': tax,
        'genome': genome,
        'mitogenome': mitogenome,
        'has_symbiont': Symbiont.objects.filter(host_species__iexact=tax.species).exists(),
        'has_transcriptome': has_transcriptome,
        'has_proteome': has_proteome,
        'has_epigenome': has_epigenome,
        'has_epigenome_csv': has_epigenome_csv,
        'epigenome_csv_records': epigenome_csv_records,
        'has_proteome_csv': has_proteome_csv,
        'proteome_csv_records': proteome_csv_records,
        'external_epigenome_link': external_epigenome_link,
        'species_dir': formatted_species,
        'proteome_link': proteome_link,
        'epigenome_link': getattr(tax, 'epigenome_link', None),
        'has_singlecell': singlecell_qs.exists(),
        'singlecell_count': singlecell_qs.count(),
        'singlecell_datasets': singlecell_qs,
        'has_umap': bool(umap_datasets),
        'umap_datasets': umap_datasets,
        'umap_datasets_json': json.dumps(umap_datasets),
        'page_title': f"{tax.species} - Species Hub",
        'busco_url': busco_url,
        'alg_url': alg_url,
        'oxford_url': oxford_url,
        'publications': _find_publications_for_taxonomy(tax),
        'gbif_id': gbif_id,
        'updated_taxonomy_line': updated_taxonomy_line,
    }
    return render(request, 'species_detail.html', context)


# ──────────────────────────────────────────────
# Genome Detail
# ──────────────────────────────────────────────

def genome_detail_page(request, taxonomy_id):
    tax = get_taxonomy_or_404(taxonomy_id)
    genome = tax.genomes.first()
    _enrich_genome_with_lopho_completeness(genome)
    formatted_species = _format_species(tax.species)

    # BUSCO
    busco_filename = f"{formatted_species}.blobtools.snail.png"
    busco_url = f"/static/plots/busco/{busco_filename}"

    # ALG & Oxford
    alg_url, oxford_url = _get_alg_oxford_urls(genome, formatted_species)

    # ── Download Files (glob pattern matching) ──
    genome_data_dir = _species_genome_dir(formatted_species)
    static_base_url = _species_genome_url(formatted_species)
    _sanitize_genome_file_names(genome_data_dir)

    def _find_matching_file(patterns, exclude_annevo=False):
        for pattern in patterns:
            for filepath in sorted(glob.glob(os.path.join(genome_data_dir, pattern))):
                basename = os.path.basename(filepath)
                if exclude_annevo and 'annevo' in basename.lower():
                    continue
                return basename
        return None

    file_type_patterns = [
        ('Genome FASTA',   ['*.genomic.fasta.masked', '*.genomic.fasta', '*.genomic.fasta.masked.gz', '*.genomic.fasta.gz', '*.fa', '*.fa.gz', '*.fasta', '*.fasta.gz'], False),
        ('Gene GFF3 (ANNEVO)', ['*ANNEVO*.gff3','*ANNEVO*.gff3.gz'], False),
        ('Gene GFF3',      ['*.gff3','*.gff3.gz'], True),
        ('Gene GTF (ANNEVO)', ['*ANNEVO*.gtf','*ANNEVO*.gtf.gz'], False),
        ('Gene GTF',       ['*.gtf','*.gtf.gz'], True),
        ('Protein FASTA (ANNEVO)', ['*ANNEVO*.pep.fasta','*ANNEVO*.pep.fasta.gz','*ANNEVO*.pep.faa','*ANNEVO*.pep.faa.gz'], False),
        ('Protein FASTA',  ['*.pep.fasta','*.pep.fasta.gz','*.pep.faa','*.pep.faa.gz'], True),
        ('CDS FASTA (ANNEVO)', ['*ANNEVO*.CDS.fasta','*ANNEVO*.CDS.fasta.gz'], False),
        ('CDS FASTA',      ['*.CDS.fasta','*.CDS.fasta.gz'], True),
        ('Gene Bed',       ['*.bed','*.bed.gz'], False),
        ('EGGNOG Annotation (ANNEVO)', [
            '*ANNEVO*.1.pep.fasta.emapper.annotations.gz',
            '*ANNEVO*.1.pep.fasta.emapper.annotations',
            '*ANNEVO*.emapper.annotations.gz',
            '*ANNEVO*.emapper.annotations',
        ], False),
        ('EGGNOG Annotation', [
            '*.1.pep.fasta.emapper.annotations.gz',
            '*.1.pep.fasta.emapper.annotations',
            '*.emapper.annotations.gz',
            '*.emapper.annotations',
        ], True),
        ('InterProScan Annotation (ANNEVO)', [
            '*ANNEVO*.1.pep.fasta.tsv.gz',
            '*ANNEVO*.1.pep.fasta.tsv',
            '*ANNEVO*.pep.fasta.tsv.gz',
            '*ANNEVO*.pep.fasta.tsv',
        ], False),
        ('InterProScan Annotation', [
            '*.1.pep.fasta.tsv.gz',
            '*.1.pep.fasta.tsv',
            '*.pep.fasta.tsv.gz',
            '*.pep.fasta.tsv',
        ], True),
        ('Repeat FAMILY', ['*-families.fa.gz'], False),
        ('Repeat LOCATION', ['*.genomic.fasta.out.gz'], False),
        ('Repeat SUMMARY', ['*.genomic.fasta.tbl'], False),
    ]

    file_downloads = []
    for label, patterns, exclude_annevo in file_type_patterns:
        matched_file = _find_matching_file(patterns, exclude_annevo=exclude_annevo)
        if matched_file:
            file_downloads.append({
                'label':     label,
                'url':       reverse('static_data_download', args=[f"data/{formatted_species}/genome/{matched_file}"]),
                'path':      f"static/data/{formatted_species}/genome/{matched_file}",
                'filename':  matched_file,
                'available': True,
            })

    # ── Annotation selection and loading ──
    annotation_dir = genome_data_dir
    annotation_eggnog_file = None
    annotation_eggnog_rows = []
    annotation_eggnog_header_cells = []
    annotation_eggnog_truncated = False

    annotation_interpro_tsv_file = None
    annotation_interpro_rows = []
    annotation_interpro_headers = []
    annotation_interpro_member_databases = []
    annotation_interpro_truncated = False

    annotation_file = request.GET.get('annotation_file')
    available_annotation_files = []

    eggnog_patterns = [
        "*.1.pep.fasta.emapper.annotations",
        "*.1.pep.fasta.emapper.annotations.gz",
        "*.emapper.annotations",
        "*.emapper.annotations.gz",
    ]
    eggnog_paths = []
    for pat in eggnog_patterns:
        eggnog_paths.extend(glob.glob(os.path.join(annotation_dir, pat)))
    eggnog_paths = sorted(set(eggnog_paths))
    for path in eggnog_paths:
        basename = os.path.basename(path)
        is_annevo = '.ANNEVO.' in basename or '.annevo.' in basename
        available_annotation_files.append({
            'type': 'eggnog',
            'label': basename,
            'display_label': 'EGGNOG Annotation' + (' (ANNEVO)' if is_annevo else ''),
            'path': path,
        })

    interpro_patterns = [
        "*.1.pep.fasta.tsv",
        "*.1.pep.fasta.tsv.gz",
        "*.pep.fasta.tsv",
        "*.pep.fasta.tsv.gz",
    ]
    interpro_paths = []
    for pat in interpro_patterns:
        interpro_paths.extend(glob.glob(os.path.join(annotation_dir, pat)))
    interpro_paths = sorted(set(interpro_paths))
    for path in interpro_paths:
        basename = os.path.basename(path)
        is_annevo = '.ANNEVO.' in basename or '.annevo.' in basename
        available_annotation_files.append({
            'type': 'interproscan',
            'label': basename,
            'display_label': 'InterProScan Annotation' + (' (ANNEVO)' if is_annevo else ''),
            'path': path,
        })

    selected_annotation_type = None
    selected_annotation_file = None
    selected_annotation_path = None

    if annotation_file:
        annotation_file = os.path.basename(annotation_file)
        selected = next((item for item in available_annotation_files if item['label'] == annotation_file), None)
        if selected:
            selected_annotation_type = selected['type']
            selected_annotation_file = selected['label']
            selected_annotation_path = selected['path']

    if not selected_annotation_path and available_annotation_files:
        selected = available_annotation_files[0]
        selected_annotation_type = selected['type']
        selected_annotation_file = selected['label']
        selected_annotation_path = selected['path']

    available_annotation_types = sorted(set(item['type'] for item in available_annotation_files))

    MAX_PREVIEW_ROWS = 100
    if selected_annotation_type == 'eggnog' and selected_annotation_path and os.path.exists(selected_annotation_path):
        annotation_eggnog_file = selected_annotation_file
        annotation_eggnog_header_cells = _read_eggnog_header(selected_annotation_path)
        annotation_eggnog_rows = []
        annotation_eggnog_truncated = False

    if selected_annotation_type == 'interproscan' and selected_annotation_path and os.path.exists(selected_annotation_path):
        annotation_interpro_tsv_file = selected_annotation_file
        annotation_interpro_headers = _read_interpro_header(selected_annotation_path)
        annotation_interpro_rows = []
        annotation_interpro_member_databases = []
        annotation_interpro_truncated = False

    context = {
        'tax': tax,
        'genome': genome,
        'page_title': f"{tax.species} - Genome Detail",
        'busco_url': busco_url,
        'alg_url': alg_url,
        'oxford_url': oxford_url,
        'file_downloads': file_downloads,
        'annotation_eggnog_file': annotation_eggnog_file,
        'annotation_eggnog_header_cells': annotation_eggnog_header_cells,
        'annotation_eggnog_rows': annotation_eggnog_rows,
        'annotation_eggnog_truncated': annotation_eggnog_truncated,
        'annotation_interpro_tsv_file': annotation_interpro_tsv_file,
        'annotation_interpro_headers': annotation_interpro_headers,
        'annotation_interpro_rows': annotation_interpro_rows,
        'annotation_interpro_member_databases': annotation_interpro_member_databases,
        'annotation_interpro_truncated': annotation_interpro_truncated,
        'selected_annotation_type': selected_annotation_type,
        'selected_annotation_file': selected_annotation_file,
        'annotation_file_options': available_annotation_files,
        'available_annotation_types': available_annotation_types,
        'species_dir': formatted_species,
        'publications': _find_publications_for_taxonomy(tax),
        'repeat_percentage': _find_repeat_for_taxonomy(tax),
    }
    return render(request, 'genome_detail.html', context)


def genome_annotation_data(request, taxonomy_id):
    tax = get_taxonomy_or_404(taxonomy_id)
    params = request.POST if request.method == 'POST' else request.GET
    annotation_type = params.get('annotation')
    annotation_file = params.get('annotation_file')
    if annotation_type not in ('eggnog', 'interproscan'):
        return JsonResponse({
            'draw': int(params.get('draw', '1') or 1),
            'recordsTotal': 0,
            'recordsFiltered': 0,
            'data': [],
            'header': [],
            'file': '',
            'error': 'Invalid annotation type',
        }, status=200)

    formatted_species = _format_species(tax.species)
    annotation_dir = _species_genome_dir(formatted_species)
    search_value = str(params.get('search[value]', '') or params.get('search', '')).strip()
    member_db = str(params.get('columns[3][search][value]', '') or params.get('member_db', '')).strip()
    try:
        draw = int(params.get('draw', '1'))
    except ValueError:
        draw = 1
    try:
        start = int(params.get('start', '0'))
    except ValueError:
        start = 0
    try:
        length = int(params.get('length', '25'))
    except ValueError:
        length = 25
    if length < 0:
        length = 25

    if annotation_type == 'eggnog':
        eggnog_path = None
        eggnog_patterns = [
            '*.1.pep.fasta.emapper.annotations',
            '*.1.pep.fasta.emapper.annotations.gz',
            '*.emapper.annotations',
            '*.emapper.annotations.gz',
        ]
        if annotation_file:
            candidate = os.path.join(annotation_dir, os.path.basename(annotation_file))
            if os.path.exists(candidate):
                eggnog_path = candidate
        if not eggnog_path:
            for pat in eggnog_patterns:
                results = glob.glob(os.path.join(annotation_dir, pat))
                if results:
                    eggnog_path = results[0]
                    break
        if not eggnog_path or not os.path.exists(eggnog_path):
            return JsonResponse({
                'draw': draw,
                'recordsTotal': 0,
                'recordsFiltered': 0,
                'data': [],
                'header': [],
                'file': '',
                'error': 'EGGNOG annotation file not found',
            }, status=200)

        header_cells = _read_eggnog_header(eggnog_path)
        records_total = _count_eggnog_rows(eggnog_path)
        if search_value:
            records_filtered = _count_eggnog_rows(eggnog_path, search_query=search_value)
        else:
            records_filtered = records_total
        rows = _read_eggnog_rows(eggnog_path, start, length, search_query=search_value)
        return JsonResponse({
            'draw': draw,
            'recordsTotal': records_total,
            'recordsFiltered': records_filtered,
            'data': rows,
            'header': header_cells,
            'file': os.path.basename(eggnog_path),
        })

    interpro_path = None
    interpro_patterns = [
        '*.1.pep.fasta.tsv',
        '*.1.pep.fasta.tsv.gz',
        '*.pep.fasta.tsv',
        '*.pep.fasta.tsv.gz',
    ]
    if annotation_file:
        candidate = os.path.join(annotation_dir, os.path.basename(annotation_file))
        if os.path.exists(candidate):
            interpro_path = candidate
    if not interpro_path:
        for pat in interpro_patterns:
            results = glob.glob(os.path.join(annotation_dir, pat))
            if results:
                interpro_path = results[0]
                break
    if not interpro_path or not os.path.exists(interpro_path):
        return JsonResponse({
            'draw': draw,
            'recordsTotal': 0,
            'recordsFiltered': 0,
            'data': [],
            'header': [],
            'file': '',
            'error': 'InterProScan annotation file not found',
        }, status=200)

    header_cells = _read_interpro_header(interpro_path)
    records_total = _count_interpro_rows(interpro_path)
    if search_value or member_db:
        records_filtered = _count_interpro_rows(interpro_path, search_query=search_value, member_db=member_db)
    else:
        records_filtered = records_total
    rows = _read_interpro_rows(interpro_path, start, length, search_query=search_value, member_db=member_db)
    member_dbs = sorted(set(r[3] for r in rows if len(r) > 3 and r[3].strip()))
    return JsonResponse({
        'draw': draw,
        'recordsTotal': records_total,
        'recordsFiltered': records_filtered,
        'data': rows,
        'member_dbs': member_dbs,
        'header': header_cells,
        'file': os.path.basename(interpro_path),
    })


# ──────────────────────────────────────────────
# Mitogenome Detail
# ──────────────────────────────────────────────

def mitogenome_detail_page(request, taxonomy_id):
    tax = get_taxonomy_or_404(taxonomy_id)
    mitogenome = tax.mitogenomes.first()
    formatted_species = _format_species(tax.species)

    folder_name = f"mitoz.{formatted_species}.mt.fasta.result"
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', formatted_species, folder_name)
    static_base_url = f"/static/data/{formatted_species}/{folder_name}"
    circos_url = f"{static_base_url}/circos.svg"

    summary1_content = None
    summary2_content = None
    summary1_cards = []
    summary2_table = None
    summary2_one_line = None

    summary1_path = os.path.join(base_dir, 'summary.1.txt')
    if os.path.exists(summary1_path):
        with open(summary1_path, 'r', encoding='utf-8') as f:
            summary1_content = f.read()

    summary2_path = os.path.join(base_dir, 'summary.2.txt')
    if os.path.exists(summary2_path):
        with open(summary2_path, 'r', encoding='utf-8') as f:
            summary2_content = f.read()

    if summary1_content:
        summary1_rows = [line.strip() for line in summary1_content.splitlines() if line.strip()]
        if summary1_rows:
            headers = summary1_rows[0].replace('#', '').split('\t')
            for row in summary1_rows[1:]:
                cells = row.split('\t')
                if not cells:
                    continue
                seq_id = cells[0]
                fields = []
                for idx, cell in enumerate(cells[1:], start=1):
                    header = headers[idx] if idx < len(headers) else f'Field {idx}'
                    fields.append({'header': header, 'value': cell, 'is_first_col': False})
                summary1_cards.append({'seq_id': seq_id, 'fields': fields})

    if summary2_content:
        summary2_rows = [line.strip() for line in summary2_content.splitlines() if line.strip()]
        if summary2_rows:
            first_line = summary2_rows[0]
            if '\t' in first_line:
                summary2_table = [[cell.strip() for cell in row.split('\t')] for row in summary2_rows]
            elif ',' in first_line:
                summary2_table = [[cell.strip() for cell in row.split(',')] for row in summary2_rows]
            else:
                summary2_table = [[row] for row in summary2_rows]
        summary2_one_line = ' '.join(summary2_rows) if summary2_rows else None

    prefix = f"mitoz_{formatted_species}.mt.fasta_mitoscaf.fa.gbf"
    file_types = [
        ('CDS',     'cds.fasta'),
        ('Gene',    'gene.fasta'),
        ('Protein', 'cds_translation.fasta'),
        ('Genome',  'fasta'),
        ('rRNA',    'rrna.fasta'),
        ('tRNA',    'trna.fasta'),
        ('Location','summary.2.txt'),
    ]
    downloads = []
    for label, suffix in file_types:
        filename = f"{prefix}.{suffix}"
        filepath = os.path.join(base_dir, filename)
        if os.path.exists(filepath):
            downloads.append({
                'label':    label,
                'url':      reverse('static_data_download', args=[f"data/{formatted_species}/{folder_name}/{filename}"]),
                'path':     f"static/data/{formatted_species}/{folder_name}/{filename}",
                'filename': filename,
            })

    context = {
        'tax': tax,
        'mitogenome': mitogenome,
        'page_title': f"{tax.species} - Mitogenome Detail",
        'circos_url': circos_url,
        'summary1_cards': summary1_cards,
        'summary2_table': summary2_table,
        'summary2_one_line': summary2_one_line,
        'downloads': downloads,
    }
    return render(request, 'mitogenome_detail.html', context)


# ──────────────────────────────────────────────
# Symbiont helpers
# ──────────────────────────────────────────────

def _get_symbionts_for_taxonomy(tax, sub=None):
    symbionts = Symbiont.objects.filter(
        host_species__iexact=tax.species
    )
    if sub:
        symbionts = symbionts.filter(sub=sub)
    symbionts = symbionts.order_by('host_species', 'sub')
    unique_symbionts = []
    seen = set()
    for s in symbionts:
        key = (
            s.host_species, s.sub, s.symbiont_type, s.species, s.size_bp,
            s.contigs, s.n50_bp, s.gc_content, s.genes,
            s.completeness, s.contamination,
            s.domain, s.phylum, s.class_name, s.order, s.family, s.genus,
        )
        if key not in seen:
            seen.add(key)
            unique_symbionts.append(s)
    return unique_symbionts


# ──────────────────────────────────────────────
# Symbiont Detail
# ──────────────────────────────────────────────

def symbiont_detail_page(request, taxonomy_id):
    tax = get_taxonomy_or_404(taxonomy_id)
    requested_sub = request.GET.get('sub', '').strip()
    symbionts = _get_symbionts_for_taxonomy(tax, requested_sub or None)
    formatted_species = _format_species(tax.species)
    raw_species = tax.species

    candidate_dirs = _symbiont_directory_candidates(tax.species, requested_sub or 'symbiont')
    base_dir = None
    static_base_url = None
    download_candidate = formatted_species
    for candidate_dir in candidate_dirs:
        if os.path.exists(candidate_dir):
            base_dir = candidate_dir
            relative_dir = os.path.relpath(candidate_dir, os.path.join(settings.BASE_DIR, 'static', 'data'))
            download_candidate = relative_dir.split(os.sep, 1)[0]
            static_base_url = f"/static/data/{relative_dir.replace(os.sep, '/')}"
            break
    if base_dir is None:
        relative_dir = os.path.join(formatted_species, requested_sub or 'symbiont')
        base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', relative_dir)
        static_base_url = f"/static/data/{relative_dir.replace(os.sep, '/')}"
        download_candidate = formatted_species

    circos_filename = _find_symbiont_file(base_dir, formatted_species, 'svg')
    circos_url = f"{static_base_url}/{circos_filename}" if circos_filename else None

    summary1_filename = _find_symbiont_file(base_dir, formatted_species, 'summary.1.txt')
    summary2_filename = _find_symbiont_file(base_dir, formatted_species, 'summary.2.txt')
    summary1_content = None
    summary2_content = None
    if summary1_filename:
        summary1_path = os.path.join(base_dir, summary1_filename)
        if os.path.exists(summary1_path):
            with open(summary1_path, 'r', encoding='utf-8') as f:
                summary1_content = f.read()
    if summary2_filename:
        summary2_path = os.path.join(base_dir, summary2_filename)
        if os.path.exists(summary2_path):
            with open(summary2_path, 'r', encoding='utf-8') as f:
                summary2_content = f.read()

    summary1_rows = [line.strip() for line in summary1_content.splitlines() if line.strip()] if summary1_content else []
    allowed_fields = {
        'assembly name': 'Assembly Name',
        'publication': 'Publication',
        'doi': 'DOI',
        'data source': 'Data Source',
    }
    seen_fields = set()
    summary1_metadata = []
    for line in summary1_rows:
        if ':' not in line:
            continue
        key, value = (part.strip() for part in line.split(':', 1))
        key_lower = key.lower()
        if key_lower in allowed_fields and key_lower not in seen_fields:
            summary1_metadata.append({'key': allowed_fields[key_lower], 'value': value})
            seen_fields.add(key_lower)

    summary2_rows = [line.strip() for line in summary2_content.splitlines() if line.strip()] if summary2_content else []
    summary2_one_line = ' '.join(summary2_rows) if summary2_rows else None
    summary2_table = None
    if summary2_rows:
        first_line = summary2_rows[0]
        if '\t' in first_line:
            summary2_table = [[cell.strip() for cell in row.split('\t')] for row in summary2_rows]
        elif ',' in first_line:
            summary2_table = [[cell.strip() for cell in row.split(',')] for row in summary2_rows]
        else:
            summary2_table = [[row] for row in summary2_rows]

    file_types = [('Genome', 'fna'), ('Protein', 'faa'), ('Gene', 'ffn'), ('GFF3', 'gff3')]
    downloads = []
    for label, suffix in file_types:
        filename = _find_symbiont_file(base_dir, formatted_species, suffix)
        if filename:
            downloads.append({
                'label':    label,
                'url':      reverse('static_data_download', args=[f"data/{relative_dir.replace(os.sep, '/')}/{filename}"]),
                'path':     f"static/data/{relative_dir.replace(os.sep, '/')}/{filename}",
                'filename': filename,
            })

    context = {
        'tax': tax,
        'symbionts': symbionts,
        'page_title': f"{tax.species} - Symbiont Detail",
        'circos_url': circos_url,
        'summary1': summary1_content,
        'summary1_metadata': summary1_metadata,
        'summary2': summary2_content,
        'summary2_one_line': summary2_one_line,
        'summary2_table': summary2_table,
        'downloads': downloads,
    }
    return render(request, 'symbiont_detail.html', context)


# ──────────────────────────────────────────────
# Transcriptome Detail (placeholder)
# ──────────────────────────────────────────────

def transcriptome_detail_page(request, taxonomy_id):
    return redirect(f"{reverse('transcriptome_page')}?taxid={taxonomy_id}")

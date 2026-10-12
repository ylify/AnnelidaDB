#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict


def safe_filename(name):
    cleaned = re.sub(r'[^a-zA-Z0-9_.\-]', '_', name)
    cleaned = re.sub(r'_{2,}', '_', cleaned)
    return cleaned.strip('_.-') or 'unnamed_gene'


def resolve_file(path, expected_name):
    path = os.path.abspath(path)
    if os.path.isdir(path):
        candidate = os.path.join(path, expected_name)
        if os.path.isfile(candidate):
            return candidate
        candidate_lower = os.path.join(path, expected_name.lower())
        if os.path.isfile(candidate_lower):
            return candidate_lower
        raise FileNotFoundError(f'{expected_name} not found under directory: {path}')

    if os.path.isfile(path):
        basename = os.path.basename(path).lower()
        if basename == expected_name.lower():
            return path
        parent = os.path.dirname(path)
        candidate = os.path.join(parent, expected_name)
        if os.path.isfile(candidate):
            return candidate
        candidate_lower = os.path.join(parent, expected_name.lower())
        if os.path.isfile(candidate_lower):
            return candidate_lower
        raise FileNotFoundError(f'{expected_name} not found in the same directory as {path}')

    raise FileNotFoundError(f'Path not found: {path}')


def load_umap_coords(csv_path):
    coords = {}
    with open(csv_path, 'r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        expected = {'Cell_Barcode', 'UMAP_1', 'UMAP_2'}
        if reader.fieldnames is None or not expected.issubset({name.strip() for name in reader.fieldnames}):
            raise ValueError('UMAP.csv missing expected columns. Expected Cell_Barcode, UMAP_1, UMAP_2.')
        for row in reader:
            barcode = (row.get('Cell_Barcode') or '').strip()
            if not barcode:
                continue
            try:
                x = float((row.get('UMAP_1') or '').strip())
                y = float((row.get('UMAP_2') or '').strip())
            except ValueError:
                continue
            cell_type = (row.get('Cell_Type') or '').strip()
            coords[barcode] = (x, y, cell_type)
    if not coords:
        raise RuntimeError('No valid UMAP coordinates found in UMAP.csv.')
    return coords


def sample_gene_counts(csv_path, target_genes=None):
    counts = Counter()
    with open(csv_path, 'r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        expected = {'Gene_ID', 'Cell_Barcode', 'Expression_Value'}
        if reader.fieldnames is None or not expected.issubset({name.strip() for name in reader.fieldnames}):
            raise ValueError('Expression.csv missing expected columns. Expected Gene_ID, Cell_Barcode, Expression_Value.')
        for row in reader:
            gene = (row.get('Gene_ID') or '').strip()
            if not gene:
                continue
            if target_genes is not None and gene not in target_genes:
                continue
            counts[gene] += 1
    return counts


def write_gene_csvs(expression_csv_path, umap_coords, output_dir, target_genes=None, overwrite=False):
    active_writers = {}
    open_files = {}
    saved_files = []
    source_rows = 0
    written_rows = 0
    skipped_missing_barcode = 0
    skipped_invalid_value = 0

    def open_writer(gene):
        file_path = os.path.join(output_dir, f'{safe_filename(gene)}.csv')
        if os.path.exists(file_path) and not overwrite:
            active_writers[gene] = None
            return None
        f = open(file_path, 'w', newline='', encoding='utf-8')
        writer = csv.writer(f)
        if umap_coords is not None:
            writer.writerow(['Cell_Barcode', 'UMAP_1', 'UMAP_2', 'Expression_Value', 'Cell_Type'])
        else:
            writer.writerow(['Cell_Barcode', 'Cluster', 'Celltype', 'Stage', 'Expression_Value'])
        active_writers[gene] = writer
        open_files[gene] = f
        saved_files.append(file_path)
        return writer

    with open(expression_csv_path, 'r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        expected = {'Gene_ID', 'Cell_Barcode', 'Expression_Value'}
        if reader.fieldnames is None or not expected.issubset({name.strip() for name in reader.fieldnames}):
            raise ValueError('Expression.csv missing expected columns. Expected Gene_ID, Cell_Barcode, Expression_Value.')
        for row in reader:
            source_rows += 1
            gene = (row.get('Gene_ID') or '').strip()
            if not gene:
                continue
            if target_genes is not None and gene not in target_genes:
                continue
            barcode = (row.get('Cell_Barcode') or '').strip()
            if not barcode:
                skipped_missing_barcode += 1
                continue
            try:
                expression_value = float((row.get('Expression_Value') or '').strip())
            except ValueError:
                skipped_invalid_value += 1
                continue
            writer = active_writers.get(gene)
            if writer is None and gene not in active_writers:
                writer = open_writer(gene)
            if writer is None:
                continue
            if umap_coords is not None:
                if barcode not in umap_coords:
                    skipped_missing_barcode += 1
                    continue
                x, y, cell_type = umap_coords[barcode]
                writer.writerow([barcode, x, y, expression_value, cell_type])
            else:
                writer.writerow([
                    barcode,
                    (row.get('Cluster') or '').strip(),
                    (row.get('Celltype') or row.get('Cell_Type') or '').strip(),
                    (row.get('Stage') or '').strip(),
                    expression_value,
                ])
            written_rows += 1

    for f in open_files.values():
        f.close()
    return saved_files, {
        'source_rows': source_rows,
        'written_rows': written_rows,
        'skipped_missing_barcode': skipped_missing_barcode,
        'skipped_invalid_value': skipped_invalid_value,
    }


def parse_gene_list(gene_list_value):
    if gene_list_value is None:
        return None
    genes = set()
    for item in gene_list_value.split(','):
        item = item.strip()
        if item:
            genes.add(item)
    return genes if genes else None


def main():
    parser = argparse.ArgumentParser(description='Generate per-gene expression CSV files from Expression.csv or Expression_with_meta.csv.')
    parser.add_argument('path', help='Path to the directory containing UMAP.csv and Expression.csv, or a direct path to one of those files.')
    parser.add_argument('--umap', default=None, help='Optional path to UMAP.csv or a directory containing it.')
    parser.add_argument('--output-dir', default='Expression', help='Output folder name under the Expression.csv directory (default: Expression).')
    parser.add_argument('--genes', help='Comma-separated list of gene IDs to export. If omitted, all genes are processed unless --max-genes is used.')
    parser.add_argument('--overwrite', action='store_true', help='Regenerate existing gene CSV files.')
    parser.add_argument('--clean', action='store_true', help='Remove existing generated CSV files before rebuilding.')
    parser.add_argument('--max-genes', type=int, default=None, help='Limit the number of genes to process (for large datasets).')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress output.')
    args = parser.parse_args()

    arg_path = os.path.abspath(args.path)
    expression_path = None
    umap_path = None

    if os.path.isdir(arg_path):
        for expression_name in ('Expression.csv', 'Expression_with_meta.csv'):
            candidate = os.path.join(arg_path, expression_name)
            if os.path.isfile(candidate):
                expression_path = candidate
                break
        umap_candidate = os.path.join(arg_path, 'UMAP.csv')
        umap_path = umap_candidate if os.path.isfile(umap_candidate) else None
    elif os.path.isfile(arg_path):
        basename = os.path.basename(arg_path).lower()
        if basename in ('expression.csv', 'expression_with_meta.csv'):
            expression_path = arg_path
            umap_dir = args.umap or os.path.dirname(arg_path)
            umap_candidate = os.path.join(umap_dir, 'UMAP.csv')
            umap_path = umap_candidate if os.path.isfile(umap_candidate) else None
        elif basename == 'umap.csv':
            umap_path = arg_path
            expr_dir = args.umap or os.path.dirname(arg_path)
            expression_path = resolve_file(expr_dir, 'Expression.csv')
        else:
            parent = os.path.dirname(arg_path)
            expression_path = resolve_file(parent, 'Expression.csv')
            umap_path = resolve_file(parent, 'UMAP.csv')
    else:
        raise FileNotFoundError(f'Path not found: {args.path}')

    if args.umap:
        umap_candidate = os.path.abspath(args.umap)
        umap_path = resolve_file(umap_candidate, 'UMAP.csv')

    if not expression_path:
        raise RuntimeError('Could not locate Expression.csv or Expression_with_meta.csv.')

    output_dir = os.path.join(os.path.dirname(expression_path), args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    if args.clean:
        for filename in os.listdir(output_dir):
            if filename.lower().endswith('.csv'):
                os.remove(os.path.join(output_dir, filename))

    target_genes = parse_gene_list(args.genes)
    if args.genes and not args.quiet:
        print(f'Loading expression data from {expression_path}')
        print(f'Exporting {len(target_genes)} specified genes from {expression_path}...')
    elif not args.quiet:
        print(f'Loading UMAP coordinates from {umap_path}' if umap_path else 'No UMAP coordinates found; exporting metadata-only expression files.')
        print(f'Exporting gene CSV files from {expression_path}...')

    umap_coords = load_umap_coords(umap_path) if umap_path else None
    if target_genes is None and args.max_genes is not None:
        gene_counts = sample_gene_counts(expression_path)
        top_genes = [gene for gene, _ in gene_counts.most_common(args.max_genes)]
        if not args.quiet:
            print(f'Selected top {len(top_genes)} genes by count.')
        target_genes = set(top_genes)

    saved_files, stats = write_gene_csvs(expression_path, umap_coords, output_dir, target_genes=target_genes, overwrite=args.overwrite)
    if not args.quiet:
        print(f'Completed: {len(saved_files)} gene CSV files written to {output_dir}')
        print(f'Rows: {stats["written_rows"]}/{stats["source_rows"]} written; '
              f'{stats["skipped_missing_barcode"]} missing UMAP barcodes; '
              f'{stats["skipped_invalid_value"]} invalid expression values.')


if __name__ == '__main__':
    main()

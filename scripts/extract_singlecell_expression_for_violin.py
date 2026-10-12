#!/usr/bin/env python3
"""Extract single-cell gene expression matched to UMAP clusters for violin plotting.

Example usage:
    python scripts/extract_singlecell_expression_for_violin.py \
        --umap static/data/Paraescarpia_echinospica/single-cell/plume/UMAP.csv \
        --expression static/data/Paraescarpia_echinospica/single-cell/plume/Expression.csv \
        --markers static/data/Paraescarpia_echinospica/single-cell/plume/Markers.csv \
        --output output/Paraescarpia_echinospica_plume

Output:
    output/Paraescarpia_echinospica_plume/gene_files/<gene>.csv
    output/Paraescarpia_echinospica_plume/summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

FIELD_GENE = re.compile(r'^Gene[_ ]?ID$', re.I)
FIELD_BARCODE = re.compile(r'Cell[_ ]?Barcode|Barcode', re.I)
FIELD_UMAP1 = re.compile(r'^UMAP_?1$', re.I)
FIELD_UMAP2 = re.compile(r'^UMAP_?2$', re.I)
FIELD_CLUSTER = re.compile(r'Cell[_ ]?Type|Cluster|cell_type', re.I)
FIELD_EXPR = re.compile(r'Expression[_ ]?Value', re.I)


def find_field_index(header: List[str], patterns: Iterable[re.Pattern]) -> Optional[int]:
    for idx, col in enumerate(header):
        for pattern in patterns:
            if pattern.search(col.strip()):
                return idx
    return None


def load_marker_genes(path: Path) -> List[str]:
    if not path.exists():
        return []
    genes = []
    with path.open(newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader, None)
        if not header:
            return []
        gene_idx = None
        for idx, col in enumerate(header):
            if FIELD_GENE.search(col):
                gene_idx = idx
                break
        if gene_idx is None:
            gene_idx = 1 if len(header) > 1 else 0
        for row in reader:
            if len(row) > gene_idx:
                gene = row[gene_idx].strip()
                if gene:
                    genes.append(gene)
    return genes


def sanitize_gene_filename(gene: str) -> str:
    sanitized = re.sub(r'[^A-Za-z0-9_.-]', '_', gene)
    return sanitized or 'unnamed_gene'


def load_umap_metadata(umap_path: Path) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    if not umap_path.exists():
        raise FileNotFoundError(f'UMAP file not found: {umap_path}')
    barcode_to_meta: Dict[str, Dict[str, float]] = {}
    missing_fields = []
    with umap_path.open(newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader, None)
        if not header:
            raise ValueError('UMAP CSV is empty or malformed')
        x_idx = find_field_index(header, [FIELD_UMAP1])
        y_idx = find_field_index(header, [FIELD_UMAP2])
        barcode_idx = find_field_index(header, [FIELD_BARCODE])
        cluster_idx = find_field_index(header, [FIELD_CLUSTER])
        if x_idx is None:
            missing_fields.append('UMAP_1')
        if y_idx is None:
            missing_fields.append('UMAP_2')
        if barcode_idx is None:
            missing_fields.append('Cell_Barcode')
        if missing_fields:
            raise ValueError(f'UMAP CSV missing required columns: {", ".join(missing_fields)}')
        for row in reader:
            if len(row) <= max(x_idx, y_idx, barcode_idx):
                continue
            barcode = row[barcode_idx].strip()
            if not barcode:
                continue
            try:
                x = float(row[x_idx])
                y = float(row[y_idx])
            except ValueError:
                continue
            cluster = row[cluster_idx].strip() if cluster_idx is not None and len(row) > cluster_idx else 'Unknown'
            barcode_to_meta[barcode] = {'cluster': cluster, 'umap1': x, 'umap2': y}
    return barcode_to_meta, []


def load_expression_matrix(expression_path: Path, barcode_to_meta: Dict[str, Dict[str, float]], gene_filter: Optional[set] = None) -> Tuple[Dict[str, List[dict]], int]:
    if not expression_path.exists():
        raise FileNotFoundError(f'Expression file not found: {expression_path}')
    gene_to_rows: Dict[str, List[dict]] = defaultdict(list)
    skipped = 0
    with expression_path.open(newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader, None)
        if not header:
            raise ValueError('Expression CSV is empty or malformed')
        gene_idx = find_field_index(header, [FIELD_GENE])
        barcode_idx = find_field_index(header, [FIELD_BARCODE])
        expr_idx = find_field_index(header, [FIELD_EXPR])
        if gene_idx is None or barcode_idx is None or expr_idx is None:
            raise ValueError('Expression CSV missing expected columns: Gene_ID, Cell_Barcode, or Expression_Value')
        for row in reader:
            if len(row) <= max(gene_idx, barcode_idx, expr_idx):
                continue
            gene = row[gene_idx].strip()
            if not gene:
                continue
            if gene_filter is not None and gene not in gene_filter:
                continue
            barcode = row[barcode_idx].strip()
            if not barcode:
                continue
            if barcode not in barcode_to_meta:
                skipped += 1
                continue
            try:
                value = float(row[expr_idx])
            except ValueError:
                skipped += 1
                continue
            meta = barcode_to_meta[barcode]
            gene_to_rows[gene].append({
                'cell_barcode': barcode,
                'cluster': meta['cluster'],
                'umap1': meta['umap1'],
                'umap2': meta['umap2'],
                'expression_value': value,
            })
    return gene_to_rows, skipped


def write_gene_files(gene_to_rows: Dict[str, List[dict]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for gene, rows in gene_to_rows.items():
        filename = sanitize_gene_filename(gene) + '.csv'
        path = output_dir / filename
        with path.open('w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Gene_ID', 'Cell_Barcode', 'Cluster', 'UMAP_1', 'UMAP_2', 'Expression_Value'])
            for row in rows:
                writer.writerow([gene, row['cell_barcode'], row['cluster'], row['umap1'], row['umap2'], row['expression_value']])


def write_summary(summary: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Extract single-cell gene expression by cluster for violin plotting.')
    parser.add_argument('--umap', required=True, help='Path to UMAP.csv')
    parser.add_argument('--expression', required=True, help='Path to Expression.csv')
    parser.add_argument('--markers', help='Optional path to Markers.csv to restrict genes')
    parser.add_argument('--output', help='Optional output directory for gene files and summary. Defaults to Expression.csv parent directory.')
    parser.add_argument('--max-genes', type=int, default=0, help='Optional maximum number of genes to process (0 = all)')
    parser.add_argument('--force', action='store_true', help='Force regeneration even if output already exists.')
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    umap_path = Path(args.umap)
    expression_path = Path(args.expression)
    markers_path = Path(args.markers) if args.markers else None
    output_root = Path(args.output) if args.output else expression_path.parent / f'{expression_path.stem}_violin_output'
    gene_out_dir = output_root / 'gene_files'
    summary_path = output_root / 'summary.json'
    if not args.force and summary_path.exists() and gene_out_dir.exists() and any(gene_out_dir.iterdir()):
        print(f'Output already exists at {output_root}. Skipping extraction. Use --force to regenerate.')
        return 0
    barcode_to_meta, _ = load_umap_metadata(umap_path)
    gene_filter = None
    marker_count = 0
    if markers_path:
        marker_genes = load_marker_genes(markers_path)
        marker_count = len(marker_genes)
        gene_filter = set(marker_genes)
    gene_to_rows, skipped = load_expression_matrix(expression_path, barcode_to_meta, gene_filter)
    if args.max_genes and args.max_genes > 0:
        gene_to_rows = {gene: gene_to_rows[gene] for gene in list(sorted(gene_to_rows))[:args.max_genes]}
    write_gene_files(gene_to_rows, gene_out_dir)
    summary = {
        'umap_path': str(umap_path),
        'expression_path': str(expression_path),
        'markers_path': str(markers_path) if markers_path else None,
        'output_gene_files': str(gene_out_dir),
        'genes_processed': len(gene_to_rows),
        'samples_per_gene': {gene: len(rows) for gene, rows in gene_to_rows.items()},
        'skipped_rows_no_barcode_match': skipped,
        'marker_gene_count': marker_count,
    }
    write_summary(summary, output_root / 'summary.json')
    print(f'Wrote {len(gene_to_rows)} gene files to {gene_out_dir}')
    print(f'Skipped {skipped} expression rows without matching UMAP barcode')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

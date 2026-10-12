# import_csv.py full code
import os
import django
import csv
import re

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'db_project.settings')
django.setup()

from api.models import Taxonomy, Genome, Mitogenome, Symbiont, SingleCell, SingleCellMarkerGene

def safe_float(val):
    if not val:
        return None
    text = str(val).strip().replace(',', '')
    if not text or text.upper() == 'NA':
        return None
    match = re.search(r'([0-9]+(?:\.[0-9]+)?)', text)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _open_csv_with_fallback(path):
    encodings = ['utf-8-sig', 'utf-8', 'cp1252', 'latin-1']
    for enc in encodings:
        try:
            f = open(path, 'r', encoding=enc, newline='')
            # read a small sample to validate encoding before parsing
            for _ in range(5):
                line = f.readline()
                if not line:
                    break
            f.seek(0)
            return f
        except UnicodeDecodeError:
            try:
                f.close()
            except Exception:
                pass
            continue
    return open(path, 'r', encoding='latin-1', newline='')


def normalize_dictreader(reader):
    if reader.fieldnames:
        reader.fieldnames = [field.strip() if field else field for field in reader.fieldnames]
    return reader


def normalize_value(value):
    if value is None:
        return ''
    text = str(value).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ('"', "'"):
        text = text[1:-1].strip()
    return text


def normalize_row_dict(row):
    if not row:
        return {}
    return { (key.strip() if key else key): normalize_value(value) for key, value in row.items() }


def get_taxonomy_by_species(row_dict):
    species_name = (row_dict.get('Species') or '').strip()
    if not species_name:
        return None
    return Taxonomy.objects.filter(species__iexact=species_name).first()


def get_taxonomy_for_genome(row_dict):
    return get_taxonomy_by_species(row_dict)


def get_or_create_taxonomy_for_singlecell(row_dict):
    species_name = (row_dict.get('Species') or '').strip()
    if not species_name:
        return None
    tax_obj = Taxonomy.objects.filter(species__iexact=species_name).first()
    if tax_obj:
        return tax_obj
    return Taxonomy.objects.create(
        species=species_name,
        class_name=(row_dict.get('Class') or '').strip() or None,
        order=(row_dict.get('Order') or '').strip() or None,
        family=(row_dict.get('Family') or '').strip() or None,
        phylum=None,
        gbif_id=None,
        worms_id=None,
        genus=None,
        ncbi_taxonomy_id=None
    )


def safe_int(val):
    if not val:
        return None
    text = str(val).strip().replace(',', '')
    return int(text) if text.isdigit() else None


def import_data():
    print("Clearing old data...")
    Genome.objects.all().delete()
    Mitogenome.objects.all().delete()
    Symbiont.objects.all().delete()
    SingleCellMarkerGene.objects.all().delete()
    SingleCell.objects.all().delete()
    Taxonomy.objects.all().delete()

    print("1. Starting Taxonomy data import...")
    with _open_csv_with_fallback('taxonomy.csv') as f:
        reader = normalize_dictreader(csv.DictReader(f))
        for row in reader:
            row = normalize_row_dict(row)
            species_name = row.get('Species')
            if not species_name:
                continue
            ncbi_id = row.get('NCBI taxonomy id') or None

            Taxonomy.objects.create(
                ncbi_taxonomy_id=ncbi_id,
                species=species_name,
                genus=row.get('Geneus') or row.get('Genus') or None,
                family=row.get('Family') or None,
                order=row.get('Order') or None,
                class_name=row.get('Class') or None,
                phylum=row.get('Phylum') or None,
                worms_id=row.get('worms_id') or row.get('WoRMS_AphiaID') or None,
                gbif_id=row.get('gbif_id') or row.get('GBIF_ID') or None,
                dna_link=row.get('DNA_link') or None,
                rna_link=row.get('RNA_link') or None,
            )
    print("Taxonomy import complete.")

    print("1b. Starting Updated Taxonomy data import...")
    with _open_csv_with_fallback('taxonomy.updated.csv') as f:
        reader = normalize_dictreader(csv.DictReader(f))
        for row in reader:
            row = normalize_row_dict(row)
            species_name = row.get('Species')
            if not species_name:
                continue

            tax_obj = get_taxonomy_by_species(row)
            updated_fields = {
                'updated_clade6': row.get('taxa6') or row.get('clade6') or None,
                'updated_clade5': row.get('taxa5') or row.get('clade5') or None,
                'updated_clade4': row.get('taxa4') or row.get('clade4') or None,
                'updated_clade3': row.get('taxa3') or row.get('clade3') or None,
                'updated_clade2': row.get('taxa2') or row.get('clade2') or None,
                'updated_clade1': row.get('taxa1') or row.get('clade1') or None,
            }

            if tax_obj:
                for field, value in updated_fields.items():
                    setattr(tax_obj, field, value)

                if not tax_obj.gbif_id:
                    tax_obj.gbif_id = (row.get('gbif_id') or '').strip() or None
                if not tax_obj.dna_link:
                    tax_obj.dna_link = (row.get('DNA_link') or '').strip() or None
                if not tax_obj.rna_link:
                    tax_obj.rna_link = (row.get('RNA_link') or '').strip() or None
                if not tax_obj.species:
                    tax_obj.species = (row.get('Species') or '').strip() or None
                tax_obj.save()
            else:
                Taxonomy.objects.create(
                    ncbi_taxonomy_id=None,
                    species=(row.get('Species') or '').strip(),
                    genus=None,
                    family=None,
                    order=None,
                    class_name=None,
                    updated_clade6=updated_fields['updated_clade6'],
                    updated_clade5=updated_fields['updated_clade5'],
                    updated_clade4=updated_fields['updated_clade4'],
                    updated_clade3=updated_fields['updated_clade3'],
                    updated_clade2=updated_fields['updated_clade2'],
                    updated_clade1=updated_fields['updated_clade1'],
                    phylum=None,
                    worms_id=None,
                    gbif_id=(row.get('gbif_id') or '').strip() or None,
                    dna_link=(row.get('DNA_link') or '').strip() or None,
                    rna_link=(row.get('RNA_link') or '').strip() or None,
                )
    print("Updated taxonomy import complete.")

    print("2. Starting Genome data import...")
    lopho_busco_map = {}
    if os.path.exists('Genome_completeness.lophotrochozoa.csv'):
        with _open_csv_with_fallback('Genome_completeness.lophotrochozoa.csv') as lf:
            lines = [line.strip() for line in lf if line.strip()]
            l_reader = list(csv.reader(lines))
            h_idx = -1
            for idx, r in enumerate(l_reader):
                if r and any('Species' in str(c) for c in r):
                    h_idx = idx
                    break
            if h_idx >= 0:
                for r in l_reader[h_idx + 1:]:
                    if r and r[0].strip():
                        sp_key = ' '.join(r[0].strip().split()).lower()
                        g_raw = r[1].strip() if len(r) > 1 else ''
                        ge_raw = r[2].strip() if len(r) > 2 else ''
                        an_raw = r[3].strip() if len(r) > 3 else ''
                        lopho_busco_map[sp_key] = {
                            'genome': None if g_raw.upper() in ('NA', '') else g_raw,
                            'gene': None if ge_raw.upper() in ('NA', '') else ge_raw,
                            'annevo': None if an_raw.upper() in ('NA', '') else an_raw,
                        }

    with _open_csv_with_fallback('genome.csv') as f:
        reader = normalize_dictreader(csv.DictReader(f))
        for row in reader:
            row = normalize_row_dict(row)
            tax_obj = get_taxonomy_for_genome(row)
            if not tax_obj:
                continue
            
            sp_key = ' '.join((row.get('Species') or '').strip().split()).lower()
            lopho_entry = lopho_busco_map.get(sp_key, {})

            Genome.objects.create(
                taxonomy=tax_obj,
                size_mb=safe_float(row.get('Size (Mb)', '')),
                assembly_level=row.get('Assembly level', ''),
                gc_content=safe_float(row.get('GC (%)', '')),
                scaffold_n50=row.get('Scaffold N50', ''),
                completeness_genome=safe_float(row.get('Completeness (genome)', '')),
                completeness_genome_raw=row.get('Completeness (genome)', ''),
                gene_number=row.get('Gene number', ''),
                completeness_gene=row.get('Completeness (gene)', ''),
                gene_number_annevo=row.get('Gene number (by ANNEVO)', ''),
                completeness_gene_annevo=row.get('Completeness (% gene, by ANNEVO)', ''),
                completeness_genome_lopho=safe_float(lopho_entry.get('genome', '')),
                completeness_genome_lopho_raw=lopho_entry.get('genome'),
                completeness_gene_lopho=lopho_entry.get('gene'),
                completeness_gene_annevo_lopho=lopho_entry.get('annevo'),
            )
    print("Genome import complete.")

    print("3. Starting Mitogenome data import...")
    with _open_csv_with_fallback('mitogenome.csv') as f:
        reader = normalize_dictreader(csv.DictReader(f))
        for row in reader:
            row = normalize_row_dict(row)
            tax_obj = get_taxonomy_by_species(row)
            if not tax_obj:
                continue

            length_str = row.get('Length (bp)', '')
            length_val = int(length_str) if length_str.isdigit() else None

            Mitogenome.objects.create(
                taxonomy=tax_obj,
                length_bp=length_val,
                topology=row.get('Topology', '').strip(),
                genbank=row.get('Genbank', '').strip(),
                refseq=row.get('RefSeq', '').strip()
            )
    print("Mitogenome import complete.")
    print("4. Starting Symbiont data import...")
    with _open_csv_with_fallback('Symbiont.csv') as f:
        reader = normalize_dictreader(csv.DictReader(f))
        for row in reader:
            row = normalize_row_dict(row)
            host = row.get('Host')
            if not host:
                continue
            
            # Safely convert integers
            def safe_int(val):
                val = val.replace(',', '').strip()
                return int(val) if val.isdigit() else None

            Symbiont.objects.create(
                host_species=host,
                sub=row.get('Sub', ''),
                symbiont_type=row.get('Type', ''),
                size_bp=row.get('Size (bp)', ''),
                contigs=safe_int(row.get('Contigs', '')),
                n50_bp=row.get('N50 (bp)', ''),
                gc_content=safe_float(row.get('GC (%)', '')),
                genes=safe_int(row.get('Genes', '')),
                completeness=safe_float(row.get('Completeness (%)', '')),
                contamination=safe_float(row.get('Contamination (%)', '')),
                domain=row.get('Domain', ''),
                phylum=row.get('Phylum', ''),
                class_name=row.get('Class', ''),
                order=row.get('Order', ''),
                family=row.get('Family', ''),
                genus=row.get('Genus', ''),
                species=row.get('Species', '')
            )
    print("Symbiont import complete.")

    print("5. Starting SingleCell data import...")
    with _open_csv_with_fallback('single-cell.csv') as f:
        reader = normalize_dictreader(csv.DictReader(f))
        for row in reader:
            row = normalize_row_dict(row)
            species = row.get('Species')
            if not species:
                continue

            taxonomy = get_or_create_taxonomy_for_singlecell(row)
            SingleCell.objects.create(
                taxonomy=taxonomy,
                species=species,
                dir=(row.get('dir') or '').strip() or None,
                tissue=(row.get('Tissue') or '').strip() or None,
                stage=(row.get('Stage') or '').strip() or None,
                clusters=(row.get('cell clusters') or row.get('Cell clusters') or '').strip() or None,
                cells=(row.get('Cells') or '').strip() or None,
                study=(row.get('Study') or '').strip() or None,
                type=(row.get('Type') or '').strip() or None,
                year=(row.get('Year') or '').strip() or None,
                platform=(row.get('Platform') or '').strip() or None,
                bioproject=(row.get('Bioproject') or '').strip() or None,
                remark=(row.get('Remark') or '').strip() or None,
                journal=(row.get('Journal') or '').strip() or None,
                doi=(row.get('DOI') or '').strip() or None,
                title=(row.get('Title') or '').strip() or None,
            )
    print("SingleCell import complete.")

if __name__ == '__main__':
    import_data()
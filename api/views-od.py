# api/views.py full code
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
import json
import os
import subprocess
from django.conf import settings
from .models import Taxonomy, Genome, Mitogenome, Symbiont

def home_page(request):
    return render(request, 'index.html')

def get_species_list(request):
    taxonomies = Taxonomy.objects.all().prefetch_related('genomes', 'mitogenomes')[:50]
    data = []
    for tax in taxonomies:
        genome_info = None
        first_genome = tax.genomes.first()
        if first_genome:
            genome_info = {
                "size_mb": first_genome.size_mb,
                "gc_content": first_genome.gc_content,
            }
            
        mito_info = None
        first_mito = tax.mitogenomes.first()
        if first_mito:
            mito_info = {
                "length_bp": first_mito.length_bp,
            }

        data.append({
            "id": tax.id,
            "ncbi_id": tax.ncbi_taxonomy_id,
            "scientific_name": tax.species,
            "class_name": tax.class_name if tax.class_name else "Unassigned",
            "family": tax.family if tax.family else "Unknown Family",
            "genome": genome_info,
            "mitogenome": mito_info,
            "image_url": "https://via.placeholder.com/300x200/203a43/ffffff?text=Annelida"
        })
    return JsonResponse({"status": "success", "data": data})

def blast_page(request):
    return render(request, 'blast.html')

@csrf_exempt
def run_blast(request):
    if request.method == 'POST':
        data = json.loads(request.body)
        sequence = data.get('sequence', '').strip()
        if not sequence.startswith('>'):
            sequence = ">User_Query\n" + sequence
        query_file = "temp_query.fasta"
        with open(query_file, "w") as f:
            f.write(sequence)
        try:
            command = ['blastn', '-query', query_file, '-db', 'blast_data/annelida_db', '-outfmt', '6 sseqid pident evalue bitscore']
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            output_text = result.stdout.strip()
            real_results = []
            if output_text:
                for line in output_text.split('\n'):
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        real_results.append({
                            "subject": parts[0], "identity": parts[1] + "%",
                            "e_value": parts[2], "bit_score": parts[3]
                        })
            if os.path.exists(query_file): os.remove(query_file)
            return JsonResponse({"status": "success", "results": real_results})
        except subprocess.CalledProcessError:
            if os.path.exists(query_file): os.remove(query_file)
            return JsonResponse({"status": "error", "message": "BLAST execution failed."})
    return JsonResponse({"status": "error", "message": "POST requests only are supported."})

# --- Added statistics endpoint ---
def get_statistics(request):
    """Compute live counts for each dataset type."""
    return JsonResponse({
        "status": "success",
        "data": {
            "taxonomy_count": Taxonomy.objects.count(),
            "genome_count": Genome.objects.count(),
            "mitogenome_count": Mitogenome.objects.count(),
            "symbiont_count": Symbiont.objects.count()
        }
    })

# api/views.py appended code (added at the end of the file)

# --- Added: rendering interfaces for each data table page ---

def taxonomy_page(request):
    # Prefetch genome and mitochondrial data to improve query performance
    taxonomies = Taxonomy.objects.prefetch_related('genomes', 'mitogenomes').all()
    
    # Build an enhanced data dictionary
    enhanced_data = []
    for tax in taxonomies:
        enhanced_data.append({
            'ncbi_taxonomy_id': tax.ncbi_taxonomy_id,
            'species': tax.species,
            'class_name': tax.class_name,
            'order': tax.order,
            'family': tax.family,
            
            # Determine if genome and mitogenome records exist (boolean True/False)
            'has_genome': tax.genomes.exists(),
            'has_mitogenome': tax.mitogenomes.exists(),
            
            # Determine whether the original CSV-provided DNA/RNA links exist
            'has_dna': bool(tax.dna_link),
            'has_rna': bool(tax.rna_link),
            
            # Pass the original links to the frontend for direct access
            'dna_link': tax.dna_link,
            'rna_link': tax.rna_link
        })
        
    return render(request, 'table_page.html', {
        'page_title': 'Taxonomy',
        'data_type': 'taxonomy',
        'data': enhanced_data
    })
def genome_page(request):
    # Get genome data and include the related Taxonomy species name (select_related improves query speed)
    data = Genome.objects.select_related('taxonomy').all()
    return render(request, 'table_page.html', {
        'page_title': 'Genome',
        'data_type': 'genome',
        'data': data
    })

def mitogenome_page(request):
    data = Mitogenome.objects.select_related('taxonomy').all()
    return render(request, 'table_page.html', {
        'page_title': 'Mitogenome',
        'data_type': 'mitogenome',
        'data': data
    })

def symbiont_page(request):
    data = Symbiont.objects.all()
    return render(request, 'table_page.html', {
        'page_title': 'Symbiont',
        'data_type': 'symbiont',
        'data': data
    })

from django.shortcuts import get_object_or_404

# api/views.py (replace species_detail_page function)
from django.shortcuts import get_object_or_404
from django.conf import settings
import os

# Ensure Symbiont is imported at the top
from .models import Taxonomy, Genome, Mitogenome, Symbiont

# Modifications to species_detail_page in api/views.py

def species_detail_page(request, taxonomy_id):
    tax = get_object_or_404(Taxonomy, ncbi_taxonomy_id=taxonomy_id)
    genome = tax.genomes.first()
    mitogenome = tax.mitogenomes.first()
    
    # Determine whether other omics data are available
    has_symbiont = Symbiont.objects.filter(host_species=tax.species).exists()
    has_transcriptome = bool(tax.rna_link)
    
    # Keep the chart scanning logic unchanged...
    raw_species = tax.species
    formatted_species = raw_species.replace("sp.", "sp_").replace(" ", "_")
    
    busco_filename = f"{formatted_species}.blobtools.snail.png"
    busco_url = f"/static/plots/busco/{busco_filename}"
    
    alg_url = None
    oxford_url = None
    if genome and genome.assembly_level and 'chromosome' in genome.assembly_level.lower():
        alg_dir = os.path.join(settings.BASE_DIR, 'static', 'plots', 'alg')
        alg_prefix = f"Annelida.{formatted_species}.ALGs.groups"
        if os.path.exists(alg_dir):
            for file in os.listdir(alg_dir):
                if file.startswith(alg_prefix) and "karyotype" in file:
                    alg_url = f"/static/plots/alg/{file}"
                    break
                    
        oxford_filename = f"Annelida.{formatted_species}.oxford_plot.pdf"
        oxford_path = os.path.join(settings.BASE_DIR, 'static', 'plots', 'alg', oxford_filename)
        if os.path.exists(oxford_path):
            oxford_url = f"/static/plots/alg/{oxford_filename}"

    context = {
        'tax': tax,
        'genome': genome,
        'mitogenome': mitogenome,
        'has_symbiont': has_symbiont,
        'has_transcriptome': has_transcriptome,
        'page_title': f"{tax.species} - Species Hub",
        'busco_url': busco_url,
        'alg_url': alg_url,
        'oxford_url': oxford_url
    }
    
    return render(request, 'species_detail.html', context)

import os
from django.conf import settings
from django.shortcuts import render, get_object_or_404
from .models import Taxonomy

def mitogenome_detail_page(request, taxonomy_id):
    tax = get_object_or_404(Taxonomy, ncbi_taxonomy_id=taxonomy_id)
    mitogenome = tax.mitogenomes.first()
    
    # Format species names: replace "sp." with "sp_" and spaces with "_"
    raw_species = tax.species
    formatted_species = raw_species.replace("sp.", "sp_").replace(" ", "_")
    
    # Define the static file directory name and base URL
    folder_name = f"mitoz.{formatted_species}.mt.fasta.result"
    base_dir = os.path.join(settings.BASE_DIR, 'static', 'data', folder_name)
    static_base_url = f"/static/data/{folder_name}/"
    
    # Build the chart URL
    circos_url = f"{static_base_url}circos.png"
    
    # Read the summary file text content (if it exists)
    summary1_content = None
    summary2_content = None
    
    sum1_path = os.path.join(base_dir, "summary.1.txt")
    if os.path.exists(sum1_path):
        with open(sum1_path, 'r', encoding='utf-8') as f:
            summary1_content = f.read()
            
    sum2_path = os.path.join(base_dir, "summary.2.txt")
    if os.path.exists(sum2_path):
        with open(sum2_path, 'r', encoding='utf-8') as f:
            summary2_content = f.read()
            
    # Build the download link list
    prefix = f"mitozoa_{formatted_species}.mt.fasta_mitoscaf.fa"
    file_types = [
        ("CDS", "cds.fasta"),
        ("Gene", "gene.fasta"),
        ("Protein", "cds_translation.fasta"),
        ("Genome", "fasta"),
        ("rRNA", "rrna.fasta"),
        ("tRNA", "trna.fasta"),
    ]
    
    downloads = []
    for label, suffix in file_types:
        file_name = f"{prefix}.{suffix}"
        downloads.append({
            "label": label,
            "url": f"{static_base_url}{file_name}",
            "filename": file_name
        })
    
    context = {
        'tax': tax,
        'mitogenome': mitogenome,
        'page_title': f"{tax.species} - Mitogenome Detail",
        'circos_url': circos_url,
        'summary1': summary1_content,
        'summary2': summary2_content,
        'downloads': downloads,
    }
    return render(request, 'mitogenome_detail.html', context)
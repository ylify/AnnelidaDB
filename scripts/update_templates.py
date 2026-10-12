from pathlib import Path

TEMPLATES = {
    'index.html': '''{% extends 'base.html' %}
{% block title %}AnnelidaDB - The Annelid Multi-Omics Database{% endblock %}
{% block content %}
<div class="container mt-5">
    <div class="row mb-4 align-items-stretch">
        <div class="col-lg-6 mb-4 mb-lg-0">
            <div class="dashboard-section h-100">
                <h3 class="fw-bold mb-4 border-bottom pb-2 border-2 border-info"><i class="fas fa-book-open text-info me-2"></i>About Annelida</h3>
                <p class="text-secondary" style="font-size: 1.05rem; line-height: 1.8; text-align: justify;">
                    The phylum Annelida, also known as segmented worms, is a major animal group with over 22,000 living species, including earthworms, polychaetes, and leeches. They inhabit a vast range of ecosystems, from deep-sea hydrothermal vents and intertidal zones to freshwater and moist terrestrial soils.
                    <br><br>
                    <strong>AnnelidaDB</strong> is dedicated to collecting and integrating multi-omics data for annelids, providing high-quality genome, mitochondrial, and transcriptome resources for evolutionary biology, developmental biology (Evo-Devo), and ecotoxicology research.
                </p>
            </div>
        </div>
        <div class="col-lg-6">
            <div class="dashboard-section h-100 d-flex flex-column justify-content-center">
                <h3 class="fw-bold mb-4 text-center"><i class="fas fa-chart-pie text-warning me-2"></i>Data Statistics</h3>
                <div class="row text-center mt-3">
                    <div class="col-4">
                        <div class="stat-number text-primary" id="stat-taxonomy">
                            <div class="spinner-border spinner-border-sm" role="status"></div>
                        </div>
                        <div class="stat-label mt-2">Species</div>
                    </div>
                    <div class="col-4 border-start border-end">
                        <div class="stat-number text-success" id="stat-genome">
                            <div class="spinner-border spinner-border-sm" role="status"></div>
                        </div>
                        <div class="stat-label mt-2">Genomes</div>
                    </div>
                    <div class="col-4">
                        <div class="stat-number text-warning" id="stat-mitogenome">
                            <div class="spinner-border spinner-border-sm" role="status"></div>
                        </div>
                        <div class="stat-label mt-2">Mitogenomes</div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <div class="row">
        <div class="col-lg-8 mb-4 mb-lg-0">
            <div class="dashboard-section h-100 p-0 overflow-hidden">
                <div class="bg-light p-3 border-bottom">
                    <h4 class="fw-bold mb-0 m-2"><i class="fas fa-project-diagram text-success me-2"></i>Phylogenetic Overview</h4>
                </div>
                <div class="p-4 text-center d-flex flex-column align-items-center justify-content-center" style="min-height: 350px;">
                    <p class="text-muted small mb-3">(Illustration only: the phylogenetic relationship between Polychaeta, Oligochaeta, and Hirudinea. An interactive JS tree viewer will be added in a future version.)</p>
                    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/6f/Annelida_phylogeny.svg/800px-Annelida_phylogeny.svg.png" alt="Annelida Phylogenetic Tree" class="img-fluid" style="max-height: 280px; opacity: 0.85;">
                </div>
            </div>
        </div>
        <div class="col-lg-4">
            <div class="dashboard-section h-100 p-0">
                <div class="bg-light p-3 border-bottom">
                    <h4 class="fw-bold mb-0 m-2"><i class="fas fa-rocket text-danger me-2"></i>Quick Tools</h4>
                </div>
                <div class="p-4">
                    <a href="/blast/" class="btn btn-outline-primary quick-btn w-100 mb-3 text-start p-3 rounded-3 d-flex align-items-center">
                        <i class="fas fa-align-left fa-2x ms-2 me-4"></i>
                        <div><h5 class="mb-1 fw-bold">BLAST</h5></div>
                    </a>
                    <a href="#" class="btn btn-outline-success quick-btn w-100 mb-3 text-start p-3 rounded-3 d-flex align-items-center">
                        <i class="fas fa-search fa-2x ms-2 me-4"></i>
                        <div><h5 class="mb-1 fw-bold">Gene Search</h5></div>
                    </a>
                    <a href="#" class="btn btn-outline-info quick-btn w-100 text-start p-3 rounded-3 d-flex align-items-center">
                        <i class="fas fa-download fa-2x ms-2 me-4"></i>
                        <div><h5 class="mb-1 fw-bold">Data Download</h5></div>
                    </a>
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}
{% block extra_scripts %}
<script>
    document.addEventListener("DOMContentLoaded", function() {
        fetch('/api/statistics/')
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    document.getElementById('stat-taxonomy').innerText = data.data.taxonomy_count;
                    document.getElementById('stat-genome').innerText = data.data.genome_count;
                    document.getElementById('stat-mitogenome').innerText = data.data.mitogenome_count;
                }
            })
            .catch(error => console.error('Failed to fetch statistics:', error));
    });
</script>
{% endblock %}
''',
    'blast.html': '''{% extends 'base.html' %}
{% block title %}BLAST - AnnelidaDB{% endblock %}
{% block content %}
<div class="hero-section">
    <div class="container text-center">
        <h2 class="fw-bold">Sequence Alignment (BLAST)</h2>
        <p class="lead text-white-50">Align nucleotide or protein sequences against the local annelid multi-omics database</p>
    </div>
</div>
<div class="container my-5">
    <div class="row justify-content-center">
        <div class="col-md-8">
            <div class="dashboard-section shadow-sm p-4">
                <form id="blastForm">
                    <div class="mb-3">
                        <label class="form-label fw-bold">Enter FASTA-formatted sequence (DNA/RNA/Protein):</label>
                        <textarea class="form-control" id="sequenceInput" rows="6" placeholder=">Sequence_1\nATGCGTACGTAGCTAGCT..."></textarea>
                    </div>
                    <div class="d-grid">
                        <button type="button" class="btn btn-primary btn-lg" id="runBlastBtn">
                            <span id="btnText">Run BLAST</span>
                            <div class="spinner-border spinner-border-sm text-light ms-2" id="loadingSpinner" role="status" style="display:none;"></div>
                        </button>
                    </div>
                </form>
            </div>
            <div class="dashboard-section shadow-sm mt-4" id="resultSection" style="display:none;">
                <h4 class="mb-3">Alignment results:</h4>
                <div class="table-responsive">
                    <table class="table table-hover table-bordered">
                        <thead class="table-dark">
                            <tr>
                                <th>Subject</th>
                                <th>Identity</th>
                                <th>E-value</th>
                                <th>Bit Score</th>
                            </tr>
                        </thead>
                        <tbody id="resultBody"></tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}
{% block extra_scripts %}
<script>
    document.getElementById('runBlastBtn').addEventListener('click', function() {
        const seq = document.getElementById('sequenceInput').value;
        if (!seq.trim()) {
            alert('Please enter a sequence.');
            return;
        }
        document.getElementById('btnText').innerText = 'BLAST in progress, please wait...';
        document.getElementById('loadingSpinner').style.display = 'inline-block';
        document.getElementById('resultSection').style.display = 'none';
        fetch('/api/run_blast/', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sequence: seq })
        })
        .then(response => response.json())
        .then(data => {
            document.getElementById('btnText').innerText = 'Run BLAST';
            document.getElementById('loadingSpinner').style.display = 'none';
            if (data.status === 'success') {
                const tbody = document.getElementById('resultBody');
                tbody.innerHTML = '';
                data.results.forEach(item => {
                    tbody.innerHTML += `\n<tr>\n<td><a href="#">${item.subject}</a></td>\n<td><span class="badge bg-success">${item.identity}</span></td>\n<td>${item.evalue}</td>\n<td>${item.bitscore}</td>\n</tr>`;
                });
                document.getElementById('resultSection').style.display = 'block';
            }
        })
        .catch(error => {
            alert('A network error occurred.');
            document.getElementById('btnText').innerText = 'Run BLAST';
            document.getElementById('loadingSpinner').style.display = 'none';
        });
    });
</script>
{% endblock %}
''',
    'table_page.html': '''{% extends 'base.html' %}
{% block title %}{{ page_title }} - AnnelidaDB{% endblock %}
{% block content %}
<div class="hero-section pb-5">
    <div class="container text-center">
        <h2 class="fw-bold">{{ page_title }}</h2>
        <p class="lead text-white-50">Browse, search, and sort complete data sets.</p>
    </div>
</div>
<div class="container mb-5">
    <div class="dashboard-section p-4">
        <div class="table-responsive">
            <table id="omicsTable" class="table table-striped table-hover align-middle w-100">
                <thead class="table-dark">
                    <tr>
                        {% if data_type == 'taxonomy' %}
                            <th>Class</th><th>Order</th><th>Family</th><th>Species</th><th>NCBI ID</th><th class="text-center" title="Whole Genome Assembly">Genome</th><th class="text-center" title="Mitochondrial Genome">Mito</th><th class="text-center" title="DNA Sequence (SRA)">DNA</th><th class="text-center" title="RNA/Transcriptome (SRA)">RNA</th>
                        {% elif data_type == 'genome' %}
                            <th>Species</th><th>Size (Mb)</th><th>GC Content (%)</th><th>Scaffold N50</th><th>Genome Completeness</th><th>Gene Number</th><th>Gene Completeness</th>
                        {% elif data_type == 'mitogenome' %}
                            <th>Species</th><th>Length (bp)</th><th>Topology</th><th>Genbank</th><th>RefSeq</th>
                        {% elif data_type == 'symbiont' %}
                            <th>Host Species</th><th>Symbiont Type</th><th>Genome Size</th><th>GC %</th><th>Completeness</th><th>Symbiont Taxonomy</th>
                        {% endif %}
                    </tr>
                </thead>
                <tbody>
                    {% for row in data %}
                    <tr>
                        {% if data_type == 'taxonomy' %}
                            <td>{{ row.class_name|default:'-' }}</td>
                            <td>{{ row.order|default:'-' }}</td>
                            <td>{{ row.family|default:'-' }}</td>
                            <td><a href="/species/{{ row.ncbi_taxonomy_id }}/" class="text-decoration-none fw-bold"><i>{{ row.species }}</i></a></td>
                            <td><a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={{ row.ncbi_taxonomy_id }}" target="_blank" class="fw-bold badge bg-dark text-decoration-none">{{ row.ncbi_taxonomy_id }}</a></td>
                            <td class="text-center">{% if row.has_genome %}<i class="fas fa-check-circle text-success fs-5" title="Available"></i>{% else %}<i class="fas fa-times text-black-50" title="Not Available"></i>{% endif %}</td>
                            <td class="text-center">{% if row.has_mitogenome %}<i class="fas fa-check-circle text-success fs-5" title="Available"></i>{% else %}<i class="fas fa-times text-black-50" title="Not Available"></i>{% endif %}</td>
                            <td class="text-center">{% if row.has_dna %}<a href="{{ row.dna_link }}" target="_blank"><i class="fas fa-external-link-alt text-primary fs-5"></i></a>{% else %}<i class="fas fa-times text-black-50"></i>{% endif %}</td>
                            <td class="text-center">{% if row.has_rna %}<a href="{{ row.rna_link }}" target="_blank"><i class="fas fa-external-link-alt text-primary fs-5"></i></a>{% else %}<i class="fas fa-times text-black-50"></i>{% endif %}</td>
                        {% elif data_type == 'genome' %}
                            <td><a href="/species/{{ row.taxonomy.ncbi_taxonomy_id }}/" class="text-decoration-none fw-bold"><i>{{ row.taxonomy.species }}</i></a></td>
                            <td><span class="badge bg-success">{{ row.size_mb|default:'-' }} Mb</span></td>
                            <td>{{ row.gc_content|default:'-' }} %</td>
                            <td>{{ row.scaffold_n50|default:'-' }}</td>
                            <td>{{ row.completeness_genome|default:'-' }} %</td>
                            <td>{{ row.gene_number|default:'-' }}</td>
                            <td>{{ row.completeness_gene|default:'-' }} %</td>
                        {% elif data_type == 'mitogenome' %}
                            <td><a href="/species/{{ row.taxonomy.ncbi_taxonomy_id }}/" class="text-decoration-none fw-bold"><i>{{ row.taxonomy.species }}</i></a></td>
                            <td><span class="badge bg-warning text-dark">{{ row.length_bp }} bp</span></td>
                            <td>{{ row.topology }}</td>
                            <td>{{ row.genbank }}</td>
                            <td>{{ row.refseq }}</td>
                        {% elif data_type == 'symbiont' %}
                            <td class="fw-bold"><a href="{% url 'symbiont_detail' row.ncbi_taxonomy_id %}" class="text-decoration-none text-primary fw-bold"><i>{{ row.host_species }}</i><i class="fas fa-arrow-right ms-1 small"></i></a></td>
                            <td><span class="badge bg-info text-dark">{{ row.symbiont_type|default:'-' }}</span></td>
                            <td>{{ row.size_bp|default:'-' }} bp</td>
                            <td>{{ row.gc_content|default:'-' }} %</td>
                            <td>{% if row.completeness >= 90 %}<span class="badge bg-success">{{ row.completeness }}%</span>{% else %}<span class="badge bg-warning text-dark">{{ row.completeness }}%</span>{% endif %}</td>
                            <td style="max-width: 250px;"><small class="text-muted"><b>{{ row.phylum|default:'-' }}</b> > {{ row.class_name|default:'-' }} > <i>{{ row.species|default:'Unknown' }}</i></small></td>
                        {% endif %}
                    </tr>
                    {% empty %}
                    <tr>
                        <td colspan="9" class="text-center text-muted py-4">
                            <i class="fas fa-inbox fa-3x mb-3 text-light"></i><br>
                            There are currently no records in the database, or this module is under construction.
                        </td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
    </div>
</div>
{% endblock %}
{% block extra_scripts %}
<script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.5/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/1.13.5/js/dataTables.bootstrap5.min.js"></script>
<script>
    $(document).ready(function() {
        $('#omicsTable').DataTable({
            pageLength: 25,
            order: [[1, 'asc']],
            language: {
                search: 'Search:',
                lengthMenu: 'Rows per page _MENU_',
                info: ' _START_ to _END_ rows，in _TOTAL_ rows'
            }
        });
    });
</script>
{% endblock %}
''',
    'synteny.html': '''{% extends 'base.html' %}
{% block title %}{{ page_title }} - AnnelidaDB{% endblock %}
{% block content %}
<div class="container">
    <div class="page-title text-center py-5">
        <h1 class="display-5 fw-bold">Synteny Oxford Plot</h1>
        <p class="text-secondary">Choose two species to automatically match existing Oxford plot files under /static/plots/Synteny/.</p>
    </div>
    <div class="row gy-4">
        <div class="col-lg-4">
            <div class="dashboard-section selector-card shadow-sm border-0">
                <div class="card-body">
                    <h5 class="card-title mb-4"><i class="fas fa-random text-secondary me-2"></i>Select species pair</h5>
                    {% if species1_options %}
                    <form id="syntenyForm" method="get" action="/synteny/">
                        <div class="mb-4">
                            <label for="species1" class="form-label">Species 1</label>
                            <select id="species1" name="species1" class="form-select" required>
                                {% for species in species1_options %}
                                <option value="{{ species }}" {% if species == selected_species1 %}selected{% endif %}>{{ species }}</option>
                                {% endfor %}
                            </select>
                        </div>
                        <div class="mb-4">
                            <label for="species2" class="form-label">Species 2</label>
                            <select id="species2" name="species2" class="form-select" required>
                                {% for species in species2_options %}
                                <option value="{{ species }}" {% if species == selected_species2 %}selected{% endif %}>{{ species }}</option>
                                {% endfor %}
                            </select>
                        </div>
                        <div class="mb-4">
                            <button type="submit" class="btn btn-primary w-100"><i class="fas fa-play me-2"></i>View matching files</button>
                        </div>
                        <div class="alert alert-info small">
                            Only combinations with existing files are shown. If you select a non-matching pair, the page will automatically revert to an available pair.
                        </div>
                    </form>
                    {% else %}
                    <div class="alert alert-warning">
                        No Oxford plot files were detected under /static/plots/Synteny/.
                    </div>
                    {% endif %}
                </div>
            </div>
        </div>
        <div class="col-lg-8">
            <div class="dashboard-section shadow-sm border-0">
                <div class="card-body">
                    <div class="d-flex justify-content-between align-items-center mb-4">
                        <div>
                            <h5 class="card-title mb-1">Oxford Plot Preview</h5>
                            <p class="mb-0 text-muted">File naming format: Species1.Species2.AnLGs.oxford_plot.pdf</p>
                        </div>
                        {% if active_filename %}
                        <span class="badge bg-success">{{ active_filename }}</span>
                        {% endif %}
                    </div>
                    {% if plot_url %}
                    <div class="border rounded overflow-hidden plot-frame" style="min-height: 760px; height: 76vh;">
                        <object data="{{ plot_url }}" type="application/pdf" width="100%" height="100%">
                            <p class="p-4 text-center text-muted">Unable to display the PDF. Please <a href="{{ plot_url }}" target="_blank">click here to download</a>.</p>
                        </object>
                    </div>
                    {% elif species1_options %}
                    <div class="alert alert-secondary text-center py-5">
                        <i class="fas fa-file-excel fa-2x mb-3 text-muted"></i>
                        <div>The selected species pair currently has no available Oxford plot file.</div>
                        <div class="mt-2 text-muted">Please choose another available pair from the left side.</div>
                    </div>
                    {% else %}
                    <div class="alert alert-warning text-center py-5">
                        <i class="fas fa-folder-open fa-2x mb-3"></i>
                        <div>No Synteny Oxford plot files were found.</div>
                    </div>
                    {% endif %}
                </div>
            </div>
        </div>
    </div>
</div>
{% endblock %}
{% block extra_scripts %}
<script>
    const speciesPairs = {{ species_pairs_json|safe }};
    const species1Select = document.getElementById('species1');
    const species2Select = document.getElementById('species2');
    const syntenyForm = document.getElementById('syntenyForm');
    function updateSpecies2Options() {
        const species1 = species1Select.value;
        const options = speciesPairs[species1] || [];
        species2Select.innerHTML = '';
        options.forEach(species2 => {
            const option = document.createElement('option');
            option.value = species2;
            option.textContent = species2;
            species2Select.appendChild(option);
        });
    }
    if (species1Select && species2Select && syntenyForm) {
        species1Select.addEventListener('change', () => {
            updateSpecies2Options();
            species2Select.selectedIndex = 0;
            syntenyForm.submit();
        });
        species2Select.addEventListener('change', () => {
            syntenyForm.submit();
        });
    }
</script>
{% endblock %}
''',
    'species_detail.html': '''{% extends 'base.html' %}
{% block title %}{{ page_title|default:"Species Hub" }} - AnnelidaDB{% endblock %}
{% block content %}
<div class="container mt-5 mb-5">
    <div class="page-title text-center mb-5">
        <h1 class="display-5 fw-bold">{{ tax.species }}</h1>
        <p class="text-secondary">The Species Hub page now follows the homepage Bootstrap style and displays multi-omics data associated with this species.</p>
    </div>
    <div class="row gx-4 gy-4 mb-4">
        <div class="col-lg-8">
            <div class="dashboard-section p-4">
                <div class="d-flex justify-content-between align-items-start mb-4">
                    <div>
                        <h2 class="h4 fw-bold mb-2">Species Overview</h2>
                        <p class="text-secondary mb-0">NCBI Taxonomy ID: <strong>{{ tax.ncbi_taxonomy_id }}</strong></p>
                    </div>
                    <a href="javascript:history.back()" class="btn btn-outline-info btn-sm"><i class="fas fa-arrow-left me-2"></i>Back</a>
                </div>
                <div class="row g-3">
                    <div class="col-sm-6"><div class="p-3 rounded-4 bg-light"><h5 class="mb-1">Phylum</h5><p class="mb-0 text-secondary">{{ tax.phylum|default:"Unassigned" }}</p></div></div>
                    <div class="col-sm-6"><div class="p-3 rounded-4 bg-light"><h5 class="mb-1">Class</h5><p class="mb-0 text-secondary">{{ tax.class_name|default:"Unassigned" }}</p></div></div>
                    <div class="col-sm-6"><div class="p-3 rounded-4 bg-light"><h5 class="mb-1">Order</h5><p class="mb-0 text-secondary">{{ tax.order|default:"Unassigned" }}</p></div></div>
                    <div class="col-sm-6"><div class="p-3 rounded-4 bg-light"><h5 class="mb-1">Family</h5><p class="mb-0 text-secondary">{{ tax.family|default:"Unassigned" }}</p></div></div>
                    <div class="col-12"><div class="p-3 rounded-4 bg-light"><h5 class="mb-1">Genus</h5><p class="mb-0 text-secondary">{{ tax.genus|default:"Unassigned" }}</p></div></div>
                </div>
            </div>
        </div>
        <div class="col-lg-4">
            <div class="dashboard-section p-4 h-100">
                <h2 class="h5 fw-bold mb-4">External Links</h2>
                <div class="list-group list-group-flush">
                    <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={{ tax.ncbi_taxonomy_id }}" target="_blank" class="list-group-item list-group-item-action px-3 py-3"><i class="fas fa-database me-2 text-info"></i>NCBI Taxonomy</a>
                    {% if tax.worms_id %}
                    <a href="https://www.marinespecies.org/aphia.php?p=taxdetails&id={{ tax.worms_id }}" target="_blank" class="list-group-item list-group-item-action px-3 py-3"><i class="fas fa-water me-2 text-primary"></i>WoRMS</a>
                    {% else %}
                    <div class="list-group-item px-3 py-3 text-muted"><i class="fas fa-ban me-2"></i>No WoRMS ID</div>
                    {% endif %}
                    {% if tax.gbif_id %}
                    <a href="https://www.gbif.org/species/{{ tax.key_field }}" target="_blank" class="list-group-item list-group-item-action px-3 py-3"><i class="fas fa-leaf me-2 text-success"></i>GBIF</a>
                    {% else %}
                    <div class="list-group-item px-3 py-3 text-muted"><i class="fas fa-ban me-2"></i>No GBIF Key</div>
                    {% endif %}
                </div>
            </div>
        </div>
    </div>
    <div class="dashboard-section p-4 mb-5">
        <div class="d-flex flex-column flex-md-row justify-content-between align-items-start gap-3 mb-4">
            <div><h2 class="h4 fw-bold mb-1">Multi-Omics Dashboard</h2><p class="text-secondary mb-0">Click the buttons below to quickly switch between omics detail views.</p></div>
            <div class="btn-group tab-btn-group" role="group" aria-label="Tabs">
                <button type="button" class="btn btn-outline-secondary active" id="tab-genome" onclick="switchTab('genome')">Genome</button>
                <button type="button" class="btn btn-outline-secondary" id="tab-mitogenome" onclick="switchTab('mitogenome')">Mitogenome</button>
                <button type="button" class="btn btn-outline-secondary" id="tab-transcriptome" onclick="switchTab('transcriptome')">Transcriptome</button>
                <button type="button" class="btn btn-outline-secondary" id="tab-singlecell" onclick="switchTab('singlecell')">Single-cell</button>
                <button type="button" class="btn btn-outline-secondary" id="tab-symbiont" onclick="switchTab('symbiont')">Symbiont</button>
            </div>
        </div>
        <div id="content-genome" class="tab-pane active">
            {% if genome %}
            <div class="text-center py-6">
                <div class="mb-4 rounded-circle bg-info bg-opacity-10 d-inline-flex align-items-center justify-content-center" style="width:88px; height:88px;">
                    <i class="fas fa-dna text-info fa-2x"></i>
                </div>
                <h3 class="h5 fw-bold mb-2">Genome Assembly Available</h3>
                <p class="text-secondary mb-4">Assembly size <strong>{{ genome.size_mb|default:'-' }} Mb</strong> · GC content <strong>{{ genome.gc_content|default:'-' }}%</strong>.</p>
                <a href="{% url 'genome_detail' tax.ncbi_taxonomy_id %}" class="btn btn-primary btn-lg px-5">View Genome Details</a>
            </div>
            {% else %}
            <div class="text-center py-6 text-secondary">
                <i class="fas fa-dna fa-3x mb-3"></i>
                <p>No Genome data available for this species.</p>
            </div>
            {% endif %}
        </div>
        <div id="content-mitogenome" class="tab-pane d-none">
            {% if mitogenome %}
            <div class="text-center py-6">
                <div class="mb-4 rounded-circle bg-info bg-opacity-10 d-inline-flex align-items-center justify-content-center" style="width:88px; height:88px;">
                    <i class="fas fa-circle-notch text-info fa-2x"></i>
                </div>
                <h3 class="h5 fw-bold mb-2">Mitogenome Records Available</h3>
                <p class="text-secondary mb-4">Length <strong>{{ mitogenome.length_bp|default:'-' }} bp</strong> · Topology <strong>{{ mitogenome.topology|default:'-' }}</strong>.</p>
                <a href="{% url 'mitogenome_detail' tax.ncbi_taxonomy_id %}" class="btn btn-primary btn-lg px-5">View Mitogenome Details</a>
            </div>
            {% else %}
            <div class="text-center py-6 text-secondary">
                <i class="fas fa-circle-notch fa-3x mb-3"></i>
                <p>No Mitogenome data available for this species.</p>
            </div>
            {% endif %}
        </div>
        <div id="content-transcriptome" class="tab-pane d-none">
            {% if has_transcriptome %}
            <div class="text-center py-6">
                <div class="mb-4 rounded-circle bg-info bg-opacity-10 d-inline-flex align-items-center justify-content-center" style="width:88px; height:88px;">
                    <i class="fas fa-wave-square text-info fa-2x"></i>
                </div>
                <h3 class="h5 fw-bold mb-2">Transcriptome / RNA Data</h3>
                <p class="text-secondary mb-4">External transcriptome data is available for this species.</p>
                <a href="{% url 'transcriptome_page' %}?taxid={{ tax.ncbi_taxonomy_id }}" class="btn btn-primary btn-lg px-5">Access Transcriptome Data</a>
            </div>
            {% else %}
            <div class="text-center py-6 text-secondary">
                <i class="fas fa-wave-square fa-3x mb-3"></i>
                <p>No Transcriptome data available for this species.</p>
            </div>
            {% endif %}
        </div>
        <div id="content-singlecell" class="tab-pane d-none">
            <div class="text-center py-6 text-secondary">
                <i class="fas fa-microscope fa-3x mb-3"></i>
                <p>Single-cell data module is currently under development.</p>
            </div>
        </div>
        <div id="content-symbiont" class="tab-pane d-none">
            {% if has_symbiont %}
            <div class="text-center py-6">
                <div class="mb-4 rounded-circle bg-info bg-opacity-10 d-inline-flex align-items-center justify-content-center" style="width:88px; height:88px;">
                    <i class="fas fa-bacterium text-info fa-2x"></i>
                </div>
                <h3 class="h5 fw-bold mb-2">Symbiont Records Found</h3>
                <p class="text-secondary mb-4">Associated symbionts are available for <strong>{{ tax.species }}</strong>.</p>
                <a href="{% url 'symbiont_detail' tax.ncbi_taxonomy_id %}" class="btn btn-primary btn-lg px-5">View Symbiont Details</a>
            </div>
            {% else %}
            <div class="text-center py-6 text-secondary">
                <i class="fas fa-bacterium fa-3x mb-3"></i>
                <p>No Symbiont records found for this host.</p>
            </div>
            {% endif %}
        </div>
    </div>
</div>
{% endblock %}
{% block extra_scripts %}
<script>
    function switchTab(tab) {
        document.querySelectorAll('.tab-pane').forEach(el => el.classList.add('d-none'));
        document.querySelectorAll('.tab-btn-group .btn').forEach(el => el.classList.remove('active'));
        document.getElementById('content-' + tab).classList.remove('d-none');
        document.getElementById('tab-' + tab).classList.add('active');
    }
</script>
{% endblock %}
''',
    'genome_detail.html': '''{% extends 'base.html' %}
{% block title %}{{ page_title|default:"Genome Detail" }} - AnnelidaDB{% endblock %}
{% block content %}
<div class="container mb-5">
    <div class="row mb-4 align-items-stretch">
        <div class="col-12">
            <div class="dashboard-section p-4 mb-4">
                <div class="d-flex flex-column flex-md-row justify-content-between align-items-start gap-3">
                    <div>
                        <h1 class="h3 fw-bold mb-2"><i>{{ tax.species }}</i></h1>
                        <p class="text-secondary mb-0">NCBI Taxonomy ID: <strong>{{ tax.ncbi_taxonomy_id }}</strong></p>
                    </div>
                    <a href="/genome/" class="btn btn-outline-info btn-sm"><i class="fas fa-arrow-left me-2"></i>Back to Genomes</a>
                </div>
            </div>
        </div>
        <div class="col-lg-4 d-flex flex-column gap-4">
            <div class="dashboard-section p-4 border-start border-info border-4">
                <h2 class="h5 fw-bold mb-3">Taxonomy</h2>
                <ul class="list-group list-group-flush">
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Phylum</span><span>{{ tax.phylum|default:'Unassigned' }}</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Class</span><span>{{ tax.class_name|default:'Unassigned' }}</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Order</span><span>{{ tax.order|default:'Unassigned' }}</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Family</span><span>{{ tax.family|default:'Unassigned' }}</span></li>
                </ul>
            </div>
            <div class="dashboard-section p-4 border-start border-success border-4 flex-grow-1">
                <h2 class="h5 fw-bold mb-3">Assembly Statistics</h2>
                {% if genome %}
                <ul class="list-group list-group-flush">
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Size</span><span>{{ genome.size_mb|default:'-' }} Mb</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Assembly Level</span><span>{{ genome.assembly_level|default:'-' }}</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">GC Content</span><span>{{ genome.gc_content|default:'-' }} %</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Scaffold N50</span><span>{{ genome.scaffold_n50|default:'-' }}</span></li>
                    <li class="list-group-item px-0 py-2 d-flex justify-content-between align-items-center"><span class="text-secondary">Gene Count</span><span>{{ genome.gene_number|default:'-' }}</span></li>
                </ul>
                <div class="mt-3 pt-3 border-top">
                    <div class="d-flex justify-content-between align-items-center mb-1"><span class="text-secondary">Genome</span><span class="fw-bold text-success">{{ genome.completeness_genome|default:'-' }}%</span></div>
                    <div class="d-flex justify-content-between align-items-center"><span class="text-secondary">Gene</span><span class="fw-bold text-info">{{ genome.completeness_gene|default:'-' }}%</span></div>
                </div>
                {% else %}
                <p class="text-muted">No genome data available.</p>
                {% endif %}
            </div>
        </div>
        <div class="col-lg-8">
            <div class="dashboard-section p-4">
                <h2 class="h5 fw-bold mb-3">BUSCO Completeness</h2>
                <div class="bg-light rounded-4 p-4 text-center">
                    {% if busco_url %}
                    <img src="{{ busco_url }}" alt="BUSCO plot" class="img-fluid rounded shadow-sm" style="max-height:550px; object-fit:contain;" onerror="this.onerror=null;this.style.display='none';">
                    {% else %}
                    <div class="py-5 text-muted fst-italic">BUSCO plot not available.</div>
                    {% endif %}
                </div>
            </div>
            <div class="row g-4">
                {% if alg_url %}
                <div class="col-md-6">
                    <div class="dashboard-section p-4 border-start border-danger border-4">
                        <h2 class="h6 fw-bold mb-3">Ancestral Linkage Groups</h2>
                        <div class="bg-light rounded-4 p-3 text-center">
                            <img src="{{ alg_url }}" alt="ALG Karyotype" class="img-fluid rounded shadow-sm" style="max-height: 400px; object-fit: contain;" onerror="this.onerror=null;this.style.display='none';">
                        </div>
                    </div>
                </div>
                {% endif %}
                {% if oxford_url %}
                <div class="col-md-6">
                    <div class="dashboard-section p-4 border-start border-secondary border-4">
                        <h2 class="h6 fw-bold mb-3">Oxford Plot Synteny</h2>
                        <div class="bg-light rounded-4 p-3 text-center">
                            <img src="{{ oxford_url }}" alt="Oxford Plot" class="img-fluid rounded shadow-sm" style="max-height: 400px; object-fit: contain;" onerror="this.onerror=null;this.style.display='none';">
                        </div>
                    </div>
                </div>
                {% endif %}
            </div>
            <div class="dashboard-section p-4 border-start border-warning border-4 mt-4">
                <h2 class="h5 fw-bold mb-3">Download Files</h2>
                <div class="row g-3">
                    {% for item in file_downloads %}
                    <div class="col-12 col-sm-6 col-xl-4">
                        {% if item.available %}
                        <a href="{{ item.url }}" download="{{ item.filename }}" class="btn btn-outline-warning w-100 text-start py-3 shadow-sm bg-white">
                            <i class="fas fa-file-arrow-down me-2"></i>{{ item.label }}
                        </a>
                        {% else %}
                        <button class="btn btn-outline-secondary w-100 text-start py-3 shadow-sm" disabled>
                            <i class="fas fa-file-arrow-down me-2"></i>{{ item.label }}
                        </button>
                        {% endif %}
                    </div>
                    {% endfor %}
                </div>
            </div>
            <div class="dashboard-section p-4 border-start border-info border-4 mt-4">
                <h2 class="h5 fw-bold mb-3">Annotation</h2>
                {% if annotation_eggnog_file or annotation_interpro_tsv_file %}
                <div class="row row-cols-1 row-cols-md-4 g-3 mb-4">
                    <div class="col">
                        <label class="form-label fw-semibold">Choose File</label>
                        <select id="annotation-file-select" class="form-select">
                            {% if annotation_eggnog_file %}<option value="eggnog">EGGNOG Annotation</option>{% endif %}
                            {% if annotation_interpro_tsv_file %}<option value="interproscan" {% if not annotation_eggnog_file %}selected{% endif %}>InterProScan TSV</option>{% endif %}
                        </select>
                    </div>
                    <div class="col" id="member-db-filter-wrapper" style="display:none;">
                        <label class="form-label fw-semibold">Member database</label>
                        <select id="member-db-filter" class="form-select"></select>
                    </div>
                </div>
                <div id="annotation-results"></div>
                {% else %}
                <div class="alert alert-secondary">No annotation files were located.</div>
                {% endif %}
            </div>
        </div>
    </div>
</div>
{% endblock %}
{% block extra_scripts %}
<script>
    const annotationData = {
        eggnog: [],
        interpro: []
    };
    document.getElementById('annotation-file-select')?.addEventListener('change', function() {
        document.getElementById('member-db-filter-wrapper').style.display = this.value === 'interproscan' ? 'block' : 'none';
    });
</script>
{% endblock %}
''',
    'mitogenome_detail.html': '''{% extends 'base.html' %}
{% block title %}{{ page_title|default:"Mitogenome Detail" }} - AnnelidaDB{% endblock %}
{% block content %}
<div class="container mt-5 mb-5">
    <div class="dashboard-section p-4 mb-4">
        <div class="d-flex flex-column flex-md-row justify-content-between align-items-start gap-3">
            <div>
                <h1 class="h3 fw-bold mb-2"><i>{{ tax.species }}</i></h1>
                <p class="text-secondary mb-0">NCBI Taxonomy ID: <strong>{{ tax.ncbi_taxonomy_id }}</strong></p>
            </div>
            <a href="/mitochondria/" class="btn btn-outline-info btn-sm"><i class="fas fa-arrow-left me-2"></i>Back to Mitogenomes</a>
        </div>
    </div>
    <div class="row align-items-stretch g-4">
        <div class="col-lg-6">
            <div class="dashboard-section p-4">
                <h2 class="h5 fw-bold mb-3">Circos Plot</h2>
                <div class="bg-light rounded-4 p-4 text-center">
                    {% if circos_url %}
                    <img src="{{ circos_url }}" alt="Circos Plot" class="img-fluid rounded shadow-sm" style="max-height: 520px; object-fit: contain;">
                    {% else %}
                    <div class="py-5 text-muted fst-italic">Circos plot file not found.</div>
                    {% endif %}
                </div>
            </div>
        </div>
        <div class="col-lg-6 d-flex flex-column gap-4">
            <div class="dashboard-section p-4 border-start border-primary border-4">
                <h2 class="h5 fw-bold mb-3">Download Fasta Files</h2>
                {% if downloads %}
                <div class="row g-2">
                    {% for item in downloads %}
                    <div class="col-6">
                        <a href="{{ item.url }}" download="{{ item.filename }}" class="btn btn-outline-primary w-100 bg-white shadow-sm fw-bold">{{ item.label }}</a>
                    </div>
                    {% endfor %}
                </div>
                {% else %}
                <div class="alert alert-secondary mb-0">No downloadable sequence files were found.</div>
                {% endif %}
            </div>
            <div class="dashboard-section p-4 border-start border-info border-4 flex-grow-1 d-flex flex-column">
                <div class="d-flex justify-content-between align-items-center mb-3">
                    <h2 class="h5 fw-bold mb-0">Gene Features Summary</h2>
                    <span class="badge bg-light text-dark border"><i class="fas fa-file-alt me-1"></i>summary.1.txt</span>
                </div>
                <div class="bg-light rounded-4 p-3 flex-grow-1 overflow-auto" style="max-height: 540px;">
                    {% if summary1_cards %}
                        {% for card in summary1_cards %}
                        <div class="gene-card mb-3">
                            <div class="gene-header d-flex justify-content-between align-items-center">
                                <span class="text-primary fw-bold">{{ card.seq_id }}</span>
                                <a href="https://www.ncbi.nlm.nih.gov/nuccore/{{ card.seq_id }}" target="_blank" class="btn btn-sm btn-outline-primary py-0 px-2">NCBI</a>
                            </div>
                            {% for field in card.fields %}
                            <div class="gene-row">
                                <div class="gene-label">{{ field.header }}</div>
                                <div class="gene-value">{% if field.is_first_col %}<a href="https://www.ncbi.nlm.nih.gov/nuccore/{{ field.value }}" target="_blank" class="seq-link">{{ field.value }}</a>{% else %}{{ field.value }}{% endif %}</div>
                            </div>
                            {% endfor %}
                        </div>
                        {% endfor %}
                    {% else %}
                    <div class="alert alert-secondary mb-0">Could not locate summary.1.txt in the directory.</div>
                    {% endif %}
                </div>
            </div>
        </div>
    </div>
    <div class="dashboard-section p-4 border-start border-warning border-4 mt-4">
        <h2 class="h5 fw-bold mb-3">Codon Usage Stats</h2>
        {% if summary2_table %}
        <div class="table-responsive bg-white rounded shadow-sm" style="max-height: 420px; overflow-y: auto;">
            <table class="table table-hover table-striped text-center align-middle mb-0" style="font-size: 0.85rem; white-space: nowrap;">
                <thead class="table-dark position-sticky top-0">
                    <tr>
                        {% for header in summary2_table.0 %}
                        <th class="fw-bold px-3 py-3 border-0">{{ header }}</th>
                        {% endfor %}
                    </tr>
                </thead>
                <tbody>
                    {% for row in summary2_table|slice:'1:' %}
                    <tr>
                        {% for cell in row %}
                        <td class="px-3 py-2">{{ cell }}</td>
                        {% endfor %}
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
        {% elif summary2_one_line %}
        <div class="alert alert-light">{{ summary2_one_line }}</div>
        {% else %}
        <div class="alert alert-secondary">Could not locate summary.2.txt.</div>
        {% endif %}
    </div>
</div>
{% endblock %}
''',
    'symbiont_detail.html': '''{% extends 'base.html' %}
{% block title %}{{ page_title|default:"Symbiont Detail" }} - AnnelidaDB{% endblock %}
{% block content %}
<div class="container mt-5 mb-5">
    <div class="dashboard-section p-4 mb-4">
        <div class="d-flex flex-column flex-md-row justify-content-between align-items-start gap-3">
            <div>
                <h1 class="h3 fw-bold mb-2">{{ tax.species }}</h1>
                <p class="text-secondary mb-0">Host NCBI Taxonomy ID: <strong>{{ tax.ncbi_taxonomy_id }}</strong></p>
            </div>
            <a href="/symbiont/" class="btn btn-outline-info btn-sm"><i class="fas fa-arrow-left me-2"></i>Back to Symbionts</a>
        </div>
    </div>
    <div class="row g-4 align-items-stretch">
        <div class="col-lg-6">
            <div class="dashboard-section p-4 border-start border-info border-4 flex-column h-100">
                <h2 class="h5 fw-bold mb-3">General Information</h2>
                <div class="table-responsive">
                    {% if symbionts %}
                        {% for sym in symbionts %}
                        <table class="table table-sm table-hover table-custom mb-4 rounded shadow-sm bg-white">
                            <tbody>
                                <tr class="table-light"><td colspan="2" class="fw-bold text-primary">{{ sym.species|default:'Unknown Symbiont' }} — <span class="text-muted fst-italic">{{ sym.symbiont_type|default:'Unknown Type' }}</span></td></tr>
                                <tr><td class="text-muted w-50">Completeness</td><td>{{ sym.completeness|default:'-' }}%</td></tr>
                                <tr><td class="text-muted">Contamination</td><td class="text-danger">{{ sym.contamination|default:'-' }}%</td></tr>
                                <tr><td class="text-muted">Size</td><td>{{ sym.size_bp|default:'-' }} bp</td></tr>
                                <tr><td class="text-muted">Contigs</td><td>{{ sym.contigs|default:'-' }}</td></tr>
                                <tr><td class="text-muted">N50</td><td>{{ sym.n50_bp|default:'-' }} bp</td></tr>
                                <tr><td class="text-muted">Genes</td><td>{{ sym.genes|default:'-' }}</td></tr>
                                <tr><td class="text-muted">Domain</td><td>{{ sym.domain|default:'-' }}</td></tr>
                                <tr><td class="text-muted">Phylum</td><td>{{ sym.phylum|default:'-' }}</td></tr>
                                <tr><td class="text-muted">Class</td><td>{{ sym.class_name|default:'-' }}</td></tr>
                                <tr><td class="text-muted">Order</td><td>{{ sym.order|default:'-' }}</td></tr>
                                <tr><td class="text-muted">Family</td><td>{{ sym.family|default:'-' }}</td></tr>
                                <tr><td class="text-muted">Genus</td><td class="fst-italic">{{ sym.genus|default:'-' }}</td></tr>
                            </tbody>
                        </table>
                        {% endfor %}
                    {% else %}
                    <div class="text-center py-5 text-muted">No General Information found for this host.</div>
                    {% endif %}
                </div>
            </div>
        </div>
        <div class="col-lg-6">
            <div class="dashboard-section p-4 border-start border-success border-4 h-100">
                <div class="d-flex justify-content-between align-items-center mb-3">
                    <h2 class="h5 fw-bold mb-0">Circos Plot</h2>
                    <span class="badge bg-light text-dark border"><i class="fas fa-chart-pie me-1"></i>BakTa</span>
                </div>
                <div class="bg-light rounded-4 p-3 text-center">
                    {% if circos_url %}
                    <img src="{{ circos_url }}" alt="Circos Plot" class="img-fluid rounded shadow-sm" style="max-height: 520px; object-fit: contain;">
                    {% else %}
                    <div class="py-5 text-muted fst-italic">Circos plot file was not found.</div>
                    {% endif %}
                </div>
            </div>
            <div class="dashboard-section p-4 border-start border-warning border-4 mt-4">
                <h2 class="h5 fw-bold mb-3">Annotation Summary</h2>
                {% if summary2_table %}
                <div class="table-responsive bg-white rounded shadow-sm" style="max-height: 400px; overflow-y: auto;">
                    <table class="table table-sm table-hover table-striped table-bordered text-center align-middle mb-0" style="font-size: 0.85rem; white-space: nowrap;">
                        <thead class="table-dark position-sticky top-0"><tr>{% for header in summary2_table.0 %}<th class="fw-bold px-3 py-2">{{ header }}</th>{% endfor %}</tr></thead>
                        <tbody>{% for row in summary2_table|slice:'1:' %}<tr>{% for cell in row %}<td class="px-3 py-2">{{ cell }}</td>{% endfor %}</tr>{% endfor %}</tbody>
                    </table>
                </div>
                {% elif summary2_one_line %}
                <div class="alert alert-light">{{ summary2_one_line }}</div>
                {% else %}
                <div class="alert alert-secondary">summary.2.txt file was not found.</div>
                {% endif %}
            </div>
            <div class="dashboard-section p-4 border-start border-primary border-4 mt-4">
                <h2 class="h5 fw-bold mb-3">Downloads</h2>
                {% if downloads %}
                <div class="row g-3">
                    {% for item in downloads %}
                    <div class="col-md-4 col-sm-6">
                        <a href="{{ item.url }}" download="{{ item.filename }}" class="btn btn-outline-primary w-100 bg-white py-3 shadow-sm d-flex justify-content-between align-items-center">
                            <span>{{ item.label }}</span><i class="fas fa-file-download"></i>
                        </a>
                    </div>
                    {% endfor %}
                </div>
                {% else %}
                <div class="alert alert-secondary">No SOB download files found.</div>
                {% endif %}
            </div>
        </div>
    </div>
</div>
{% endblock %}
'''
}

root = Path(__file__).resolve().parent.parent
templates_dir = root / 'templates'

for name, content in TEMPLATES.items():
    target = templates_dir / name
    target.write_text(content, encoding='utf-8')
    print(f'Wrote {target}')

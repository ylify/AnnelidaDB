# api/models.py full code
import re
import uuid
from django.db import models

# 1. Taxonomy (taxonomy information table)
class Taxonomy(models.Model):
    ncbi_taxonomy_id = models.CharField(max_length=50, null=True, unique=True, verbose_name="NCBI Taxonomy ID")
    species = models.CharField(max_length=100, verbose_name="Species")
    gbif_id = models.CharField(max_length=50, verbose_name="GBIF ID", null=True, blank=True)
    worms_id = models.CharField(max_length=50, verbose_name="WoRMS ID", null=True, blank=True)
    genus = models.CharField(max_length=100, verbose_name="Genus", null=True, blank=True)
    family = models.CharField(max_length=100, verbose_name="Family", null=True, blank=True)
    order = models.CharField(max_length=100, verbose_name="Order", null=True, blank=True)
    class_name = models.CharField(max_length=100, verbose_name="Class", null=True, blank=True)
    updated_clade6 = models.CharField(max_length=100, verbose_name="Updated Clade 6", null=True, blank=True)
    updated_clade5 = models.CharField(max_length=100, verbose_name="Updated Clade 5", null=True, blank=True)
    updated_clade4 = models.CharField(max_length=100, verbose_name="Updated Clade 4", null=True, blank=True)
    updated_clade3 = models.CharField(max_length=100, verbose_name="Updated Clade 3", null=True, blank=True)
    updated_clade2 = models.CharField(max_length=100, verbose_name="Updated Clade 2", null=True, blank=True)
    updated_clade1 = models.CharField(max_length=100, verbose_name="Updated Clade 1", null=True, blank=True)
    phylum = models.CharField(max_length=100, verbose_name="Phylum", null=True, blank=True)
    
    # Additional link information
    dna_link = models.URLField(max_length=500, null=True, blank=True, verbose_name="DNA")
    rna_link = models.URLField(max_length=500, null=True, blank=True, verbose_name="RNA")
    proteome_link = models.URLField(max_length=500, null=True, blank=True, verbose_name="Proteome")
    epigenome_link = models.URLField(max_length=500, null=True, blank=True, verbose_name="Epigenome")

    def __str__(self):
        return self.species

# 2. Genome (genome data table)
class Genome(models.Model):
    # Linked to the taxonomy table; one species can have one or more assemblies
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='genomes')
    
    size_mb = models.FloatField(null=True, blank=True, verbose_name="Size (Mb)")
    assembly_level = models.CharField(max_length=100, null=True, blank=True, verbose_name="Assembly Level")
    gc_content = models.FloatField(null=True, blank=True, verbose_name="GC Content (%)")
    scaffold_n50 = models.CharField(max_length=100, null=True, blank=True, verbose_name="Scaffold N50")
    completeness_genome = models.FloatField(null=True, blank=True, verbose_name="Completeness (Genome) %")
    completeness_genome_raw = models.CharField(max_length=200, null=True, blank=True, verbose_name="Completeness (Genome) raw")
    gene_number = models.CharField(max_length=50, null=True, blank=True, verbose_name="Gene Number")
    completeness_gene = models.CharField(max_length=200, null=True, blank=True, verbose_name="Completeness (Gene) %")
    gene_number_annevo = models.CharField(max_length=50, null=True, blank=True, verbose_name="Gene Number (by ANNEVO)")
    completeness_gene_annevo = models.CharField(max_length=200, null=True, blank=True, verbose_name="Completeness (% gene, by ANNEVO)")

    completeness_genome_lopho = models.FloatField(null=True, blank=True, verbose_name="Completeness (Genome, Lophotrochozoa) %")
    completeness_genome_lopho_raw = models.CharField(max_length=200, null=True, blank=True, verbose_name="Completeness (Genome, Lophotrochozoa) raw")
    completeness_gene_lopho = models.CharField(max_length=200, null=True, blank=True, verbose_name="Completeness (Gene, Lophotrochozoa)")
    completeness_gene_annevo_lopho = models.CharField(max_length=200, null=True, blank=True, verbose_name="Completeness (% gene, by ANNEVO, Lophotrochozoa)")

    def __str__(self):
        return f"{self.taxonomy.species} - Genome"

    @staticmethod
    def _extract_busco_s_d(raw_value):
        if not raw_value:
            return None
        s_match = re.search(r'S:([\d\.]+%)', raw_value)
        d_match = re.search(r'D:([\d\.]+%)', raw_value)
        if s_match and d_match:
            return f"S:{s_match.group(1)},D:{d_match.group(1)}"
        return raw_value.strip() or None

    @property
    def completeness_genome_sd(self):
        return self._extract_busco_s_d(self.completeness_genome_raw)

    @property
    def completeness_gene_sd(self):
        return self._extract_busco_s_d(self.completeness_gene)

    @property
    def completeness_gene_annevo_sd(self):
        return self._extract_busco_s_d(self.completeness_gene_annevo)

    @property
    def completeness_genome_lopho_sd(self):
        return self._extract_busco_s_d(self.completeness_genome_lopho_raw)

    @property
    def completeness_gene_lopho_sd(self):
        return self._extract_busco_s_d(self.completeness_gene_lopho)

    @property
    def completeness_gene_annevo_lopho_sd(self):
        return self._extract_busco_s_d(self.completeness_gene_annevo_lopho)

# 3. Mitogenome (mitochondrial genome table)
class Mitogenome(models.Model):
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='mitogenomes')
    
    length_bp = models.IntegerField(null=True, blank=True, verbose_name="Length (bp)")
    topology = models.CharField(max_length=50, null=True, blank=True, verbose_name="Topology")
    genbank = models.CharField(max_length=100, null=True, blank=True, verbose_name="Genbank")
    refseq = models.CharField(max_length=100, null=True, blank=True, verbose_name="RefSeq")

    def __str__(self):
        return f"{self.taxonomy.species} - Mitogenome"

# 4. Symbiont (symbiont data table)
class Symbiont(models.Model):
    # Here we use the host name as the association; if unsure whether it exists in the Taxonomy table, a CharField can be used
    host_species = models.CharField(max_length=150, verbose_name="Host Species")
    sub = models.CharField(max_length=150, null=True, blank=True, verbose_name="Sub")
    symbiont_type = models.CharField(max_length=50, null=True, blank=True, verbose_name="Type")
    size_bp = models.CharField(max_length=50, null=True, blank=True, verbose_name="Size (bp)")
    contigs = models.IntegerField(null=True, blank=True, verbose_name="Contigs")
    n50_bp = models.CharField(max_length=50, null=True, blank=True, verbose_name="N50 (bp)")
    gc_content = models.FloatField(null=True, blank=True, verbose_name="GC %")
    genes = models.IntegerField(null=True, blank=True, verbose_name="Genes")
    completeness = models.FloatField(null=True, blank=True, verbose_name="Completeness %")
    contamination = models.FloatField(null=True, blank=True, verbose_name="Contamination %")
    
    # Taxonomic information for the symbiont
    domain = models.CharField(max_length=100, null=True, blank=True, verbose_name="Domain")
    phylum = models.CharField(max_length=100, null=True, blank=True, verbose_name="Phylum")
    class_name = models.CharField(max_length=100, null=True, blank=True, verbose_name="Class")
    order = models.CharField(max_length=100, null=True, blank=True, verbose_name="Order")
    family = models.CharField(max_length=100, null=True, blank=True, verbose_name="Family")
    genus = models.CharField(max_length=100, null=True, blank=True, verbose_name="Genus")
    species = models.CharField(max_length=100, null=True, blank=True, verbose_name="Symbiont Species")

    def __str__(self):
        return f"{self.host_species} - {self.symbiont_type}"


class SingleCell(models.Model):
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.SET_NULL, related_name='singlecell_datasets', null=True, blank=True)
    species = models.CharField(max_length=200, verbose_name="Species")
    dir = models.CharField(max_length=250, null=True, blank=True, verbose_name="Data directory")
    tissue = models.CharField(max_length=200, null=True, blank=True, verbose_name="Tissue")
    stage = models.CharField(max_length=200, null=True, blank=True, verbose_name="Stage")
    clusters = models.CharField(max_length=100, null=True, blank=True, verbose_name="Cell clusters")
    cells = models.CharField(max_length=100, null=True, blank=True, verbose_name="Cells")
    study = models.CharField(max_length=250, null=True, blank=True, verbose_name="Study")
    type = models.CharField(max_length=100, null=True, blank=True, verbose_name="Type")
    year = models.CharField(max_length=25, null=True, blank=True, verbose_name="Year")
    platform = models.CharField(max_length=150, null=True, blank=True, verbose_name="Platform")
    bioproject = models.CharField(max_length=150, null=True, blank=True, verbose_name="BioProject")
    remark = models.TextField(null=True, blank=True, verbose_name="Remark")
    journal = models.CharField(max_length=250, null=True, blank=True, verbose_name="Journal")
    doi = models.CharField(max_length=250, null=True, blank=True, verbose_name="DOI")
    title = models.TextField(null=True, blank=True, verbose_name="Title")

    @property
    def class_name(self):
        return self.taxonomy.class_name if self.taxonomy else None

    @property
    def order(self):
        return self.taxonomy.order if self.taxonomy else None

    @property
    def family(self):
        return self.taxonomy.family if self.taxonomy else None

    def __str__(self):
        return f"{self.species} ({self.year or 'unknown'})"


class SingleCellMarkerGene(models.Model):
    singlecell = models.ForeignKey(SingleCell, on_delete=models.CASCADE, related_name='marker_genes')
    cluster = models.CharField(max_length=100, verbose_name="Cluster")
    genes = models.TextField(verbose_name="Genes")

    def __str__(self):
        return f"{self.singlecell.species} - {self.cluster}"


class ContigMapping(models.Model):
    contig_id = models.CharField(max_length=100, verbose_name="Contig ID")
    species = models.CharField(max_length=200, verbose_name="Species")
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='contig_mappings')

    class Meta:
        unique_together = ('contig_id', 'taxonomy')

    def __str__(self):
        return f"{self.contig_id} - {self.species}"


class SequenceFileIndex(models.Model):
    SEQ_TYPE_CHOICES = [
        ('protein', 'Protein'),
        ('cds', 'CDS'),
    ]

    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='sequence_file_indexes')
    seq_type = models.CharField(max_length=10, choices=SEQ_TYPE_CHOICES, verbose_name='Sequence Type')
    seq_id = models.CharField(max_length=200, verbose_name='Sequence ID', db_index=True)
    file_name = models.CharField(max_length=300, verbose_name='Source FASTA Filename')
    file_path = models.CharField(max_length=500, verbose_name='Indexed FASTA Path')

    class Meta:
        unique_together = ('taxonomy', 'seq_type', 'seq_id')
        indexes = [
            models.Index(fields=['seq_type', 'seq_id']),
        ]

    def __str__(self):
        return f"{self.seq_type}:{self.seq_id} -> {self.file_path}"


class EggnogAnnotation(models.Model):
    query = models.CharField(max_length=100, verbose_name="Query ID", db_index=True)
    seed_ortholog = models.CharField(max_length=200, null=True, blank=True, verbose_name="Seed Ortholog")
    evalue = models.FloatField(null=True, blank=True, verbose_name="E-value")
    score = models.FloatField(null=True, blank=True, verbose_name="Score")
    eggnog_ogs = models.TextField(null=True, blank=True, verbose_name="EggNOG OGs")
    max_annot_lvl = models.CharField(max_length=50, null=True, blank=True, verbose_name="Max Annotation Level")
    cog_category = models.CharField(max_length=10, null=True, blank=True, verbose_name="COG Category")
    description = models.TextField(null=True, blank=True, verbose_name="Description", db_index=True)
    preferred_name = models.CharField(max_length=100, null=True, blank=True, verbose_name="Preferred Name", db_index=True)
    gos = models.TextField(null=True, blank=True, verbose_name="GOs")
    ec = models.CharField(max_length=20, null=True, blank=True, verbose_name="EC")
    kegg_ko = models.TextField(null=True, blank=True, verbose_name="KEGG KO")
    kegg_pathway = models.TextField(null=True, blank=True, verbose_name="KEGG Pathway")
    kegg_module = models.TextField(null=True, blank=True, verbose_name="KEGG Module")
    kegg_reaction = models.TextField(null=True, blank=True, verbose_name="KEGG Reaction")
    kegg_rclass = models.TextField(null=True, blank=True, verbose_name="KEGG RClass")
    brite = models.TextField(null=True, blank=True, verbose_name="BRITE")
    kegg_tc = models.TextField(null=True, blank=True, verbose_name="KEGG TC")
    cazys = models.TextField(null=True, blank=True, verbose_name="CAZys")
    bigg_reaction = models.TextField(null=True, blank=True, verbose_name="BiGG Reaction")
    pfams = models.TextField(null=True, blank=True, verbose_name="PFAMs")
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='eggnog_annotations')

    def __str__(self):
        return f"{self.query} - {self.description}"


class InterproScanAnnotation(models.Model):
    query = models.CharField(max_length=200, verbose_name='Protein ID', db_index=True)
    md5 = models.CharField(max_length=64, null=True, blank=True, verbose_name='MD5')
    length = models.IntegerField(null=True, blank=True, verbose_name='Length')
    analysis = models.CharField(max_length=100, null=True, blank=True, verbose_name='Analysis')
    signature_accession = models.CharField(max_length=100, null=True, blank=True, verbose_name='Signature Accession')
    signature_description = models.TextField(null=True, blank=True, verbose_name='Signature Description')
    start = models.IntegerField(null=True, blank=True, verbose_name='Start')
    end = models.IntegerField(null=True, blank=True, verbose_name='End')
    score = models.FloatField(null=True, blank=True, verbose_name='Score')
    status = models.CharField(max_length=50, null=True, blank=True, verbose_name='Status')
    date = models.CharField(max_length=100, null=True, blank=True, verbose_name='Date')
    ipr_accession = models.CharField(max_length=100, null=True, blank=True, verbose_name='IPR Accession')
    ipr_description = models.TextField(null=True, blank=True, verbose_name='IPR Description')
    gos = models.TextField(null=True, blank=True, verbose_name='GO')
    pathway = models.TextField(null=True, blank=True, verbose_name='Pathway')
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='interpro_annotations')

    def __str__(self):
        return f"{self.query} - {self.analysis} ({self.signature_accession})"


class SubmissionRequest(models.Model):
    STATUS_PENDING = 'pending'
    STATUS_REVIEWED = 'reviewed'
    STATUS_REJECTED = 'rejected'
    STATUS_CHOICES = [
        (STATUS_PENDING, 'Pending Review'),
        (STATUS_REVIEWED, 'Reviewed'),
        (STATUS_REJECTED, 'Rejected'),
    ]

    DATA_TYPE_CHOICES = [
        ('taxonomy', 'Taxonomy'),
        ('genome', 'Genome'),
        ('proteome', 'Proteome'),
        ('epigenome', 'Epigenome'),
        ('mitogenome', 'Mitochondria'),
        ('transcriptome', 'Transcriptome'),
        ('singlecell', 'Single-cell'),
        ('symbiont', 'Symbiont'),
        ('synteny', 'Synteny'),
        ('other', 'Other'),
    ]

    submitted_at = models.DateTimeField(auto_now_add=True, verbose_name='Submitted At')
    data_type = models.CharField(max_length=50, choices=DATA_TYPE_CHOICES, default='other', verbose_name='Data Type')
    title = models.CharField(max_length=255, verbose_name='Title')
    submitter_name = models.CharField(max_length=120, blank=True, verbose_name='Submitter Name')
    submitter_email = models.EmailField(blank=True, verbose_name='Submitter Email')
    source_url = models.CharField(max_length=512, blank=True, verbose_name='Source URL')
    details = models.TextField(verbose_name='Details')
    access_token = models.CharField(max_length=64, unique=True, editable=False, default=uuid.uuid4, verbose_name='Access Token')
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default=STATUS_PENDING, verbose_name='Status')
    admin_note = models.TextField(blank=True, verbose_name='Admin Note')
    uploaded_files = models.JSONField(blank=True, default=list, verbose_name='Uploaded Files')

    class Meta:
        ordering = ['-submitted_at']
        verbose_name = 'Submission Request'
        verbose_name_plural = 'Submission Requests'

    def __str__(self):
        return f"{self.title} ({self.get_data_type_display()})"

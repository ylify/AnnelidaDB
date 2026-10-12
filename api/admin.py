# api/admin.py full code
from django.contrib import admin
from .models import Taxonomy, Genome, Mitogenome, SubmissionRequest

class GenomeInline(admin.StackedInline):
    model = Genome
    extra = 0

class MitogenomeInline(admin.StackedInline):
    model = Mitogenome
    extra = 0

class TaxonomyAdmin(admin.ModelAdmin):
    list_display = ('species', 'class_name', 'order', 'family', 'updated_clade1', 'updated_clade6', 'ncbi_taxonomy_id')
    search_fields = ('species', 'ncbi_taxonomy_id')
    inlines = [GenomeInline, MitogenomeInline]

class SubmissionRequestAdmin(admin.ModelAdmin):
    list_display = ('title', 'data_type', 'status', 'submitted_at', 'submitter_email')
    list_filter = ('data_type', 'status', 'submitted_at')
    search_fields = ('title', 'details', 'submitter_name', 'submitter_email')
    readonly_fields = ('submitted_at',)

admin.site.register(Taxonomy, TaxonomyAdmin)
admin.site.register(Genome)
admin.site.register(Mitogenome)
admin.site.register(SubmissionRequest, SubmissionRequestAdmin)
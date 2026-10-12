"""
URL configuration for db_project project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
# db_project/urls.py full code
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path
from api import views

urlpatterns = [
    # ==========================
    # 1. Home and API endpoints
    # ==========================
    path('', views.home_page, name='home'),
    path('api/species/', views.get_species_list, name='get_species_list'),
    path('api/species/query/', views.query_species_by_name, name='query_species_by_name'),
    path('api/species/search/', views.search_species, name='search_species'),
    path('api/phylogeny/', views.get_phylogeny_tree, name='get_phylogeny_tree'),
    path('api/statistics/', views.get_statistics, name='get_statistics'),
    path('api/statistics/refresh/', views.refresh_statistics, name='refresh_statistics'),
    path('api/transcriptome/data/', views.transcriptome_data, name='transcriptome_data'),
    
    # ==========================
    # 2. BLAST-related pages
    # ==========================
    path('blast/', views.blast_page, name='blast_page'),
    path('tools/', views.tools_page, name='tools_page'),
    path('gene-family/', views.gene_family_page, name='gene_family_page'),
    path('gene-search/', views.gene_search_page, name='gene_search_page'),
    path('species/<str:taxonomy_id>/gene/<str:gene_id>/', views.gene_detail_page, name='species_gene_detail'),
    path('gene/<str:gene_id>/', views.gene_detail_page, name='gene_detail'),
    path('jbrowse/', views.jbrowse_page, name='jbrowse_page'),
    path('jbrowse-config/', views.jbrowse_config, name='jbrowse_config'),
    path('jbrowse-asset/<str:asset_name>/', views.jbrowse_asset, name='jbrowse_asset'),
    path('api/run_blast/', views.run_blast, name='run_blast'),
    path('api/blast_status/<str:job_id>/', views.blast_status, name='blast_status'),
    path('api/gene_annotation/', views.gene_annotation, name='gene_annotation'),
    path('api/gene_detail_extras/', views.gene_detail_extras, name='gene_detail_extras'),
    path('api/gene_family/', views.gene_family_api, name='gene_family_api'),
    path('gene-family/<str:family_id>/', views.gene_family_detail_page, name='gene_family_detail_page'),
    path('api/blast_job/<str:job_id>/<str:filename>/', views.blast_job_file, name='blast_job_file'),
    path('api/gene_search/', views.gene_search_api, name='gene_search_api'),

    # ==========================
    # 3. Data table overview pages
    # ==========================
    path('taxonomy/', views.taxonomy_page, name='taxonomy_page'),
    path('genome/', views.genome_page, name='genome_page'),
    path('proteome/', views.proteome_page, name='proteome_page'),
    path('epigenome/', views.epigenome_page, name='epigenome_page'),
    path('epigenome/<str:paper_id>/', views.epigenome_detail_page, name='epigenome_detail_page'),
    path('proteome/<str:paper_id>/', views.proteome_detail_page, name='proteome_detail_page'),
    path('mitochondria/', views.mitogenome_page, name='mitogenome_page'),
    path('transcriptome/', views.transcriptome_page, name='transcriptome_page'),
    path('singlecell/', views.singlecell_page, name='singlecell_page'),
    path('singlecell/detail/<int:pk>/', views.singlecell_detail_page, name='singlecell_detail'),
    path('symbiont/', views.symbiont_page, name='symbiont_page'),
    path('download/', views.download_page, name='download_page'),
    path('download/batch/', views.download_batch, name='download_batch'),
    path('download/batch-job/', views.download_batch_job, name='download_batch_job'),
    path('download/batch-job-status/<str:job_id>/', views.download_job_status, name='download_job_status'),
    path('download/batch-job-file/<str:job_id>/<str:filename>', views.download_job_file, name='download_job_file'),
    path('download/file-options/', views.download_file_options, name='download_file_options'),
    path('download/file-types/', views.download_file_types, name='download_file_types'),
    path('download/static-data/<path:file_path>/', views.static_data_download, name='static_data_download'),
    path('submit/', views.submit_data_page, name='submit_page'),
    path('submit/view/<str:token>/', views.submission_request_detail, name='submission_request_detail'),
    path('about/', views.about_page, name='about_page'),
    path('user-guide/', views.user_guide, name='user_guide'),
    path('user-guide/pdf/', views.user_guide_pdf, name='user_guide_pdf'),
    path('about/file-preview/<str:file_key>/', views.about_file_preview, name='about_file_preview'),
    path('about/file-download/<str:file_key>/', views.about_file_download, name='about_file_download'),
    path('synteny/', views.synteny_page, name='synteny_page'),
    path('synteny-upload/', views.synteny_upload_page, name='synteny_upload_page'),
    path('synteny-upload/download-selected/', views.synteny_upload_download_selected, name='synteny_upload_download_selected'),
    path('synteny-upload/edit-seqids/', views.synteny_upload_edit_seqids_page, name='synteny_upload_edit_seqids_page'),

    # ==========================
    # 4. Hub pages (previously species detail pages)
    # ==========================
    path('species/<str:taxonomy_id>/', views.species_detail_page, name='species_detail'),
    path('api/backbone/distribution/<str:taxonomy_id>/', views.backbone_distribution_json, name='backbone_distribution_json'),
    path('api/worms/taxonomy/<int:worms_id>/', views.worms_taxonomy_json, name='worms_taxonomy_json'),
    path('api/ncbi/taxonomy/<int:ncbi_id>/', views.ncbi_taxonomy_json, name='ncbi_taxonomy_json'),
    path('api/obis/depth/<int:worms_id>/', views.obis_depth_json, name='obis_depth_json'),
    path('api/gbif/species/<int:gbif_id>/', views.gbif_species_json, name='gbif_species_json'),
    path('api/gbif/occurrences/<int:gbif_id>/', views.gbif_occurrences_json, name='gbif_occurrences_json'),
    path('api/gbif/images/<int:gbif_id>/', views.gbif_images_json, name='gbif_images_json'),

    # ==========================
    # 5. Four independent omics detail pages
    # ==========================
    path('genome/detail/<str:taxonomy_id>/', views.genome_detail_page, name='genome_detail'),
    path('genome/detail/<str:taxonomy_id>/annotation-data/', views.genome_annotation_data, name='genome_annotation_data'),
    path('mitogenome/detail/<str:taxonomy_id>/', views.mitogenome_detail_page, name='mitogenome_detail'),
    path('symbiont/detail/<str:taxonomy_id>/', views.symbiont_detail_page, name='symbiont_detail'),
    path('transcriptome/detail/<str:taxonomy_id>/', views.transcriptome_detail_page, name='transcriptome_detail'),
]

if settings.DEBUG:
    static_doc_root = settings.STATICFILES_DIRS[0] if getattr(settings, 'STATICFILES_DIRS', None) else settings.STATIC_ROOT
    urlpatterns += static(settings.STATIC_URL, document_root=static_doc_root)
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

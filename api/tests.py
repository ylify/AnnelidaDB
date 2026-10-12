from django.test import TestCase, RequestFactory
from unittest.mock import patch
import os
import tempfile

from .models import Mitogenome, Taxonomy
from .views import _build_ncbi_omics_cache, _download_file_label, _download_file_type_order, _get_download_job, _get_dna_species_count, _mitogenome_species_count, _ncbi_sra_ids, _set_download_job, species_detail_page


class DownloadJobStateTests(TestCase):
    def test_download_job_reads_latest_state_after_worker_cache(self):
        job_id = 'test-download-job'
        with tempfile.TemporaryDirectory() as temp_dir, patch('api.views.DOWNLOAD_JOB_DIR', temp_dir):
            _set_download_job(job_id, status='pending', filename='', file_path='')
            self.assertEqual(_get_download_job(job_id)['status'], 'pending')

            zip_path = os.path.join(temp_dir, 'download.zip')
            _set_download_job(job_id, status='completed', filename='download.zip', file_path=zip_path)
            _get_download_job(job_id)['status'] = 'pending'

            job = _get_download_job(job_id)
            self.assertEqual(job['status'], 'completed')
            self.assertEqual(job['filename'], 'download.zip')


class DownloadFileLabelTests(TestCase):
    def test_genome_labels_match_detail_page_types(self):
        self.assertEqual(_download_file_label('species.genomic.fasta', 'genome'), 'Genome FASTA')
        self.assertEqual(_download_file_label('species.ANNEVO.gff3', 'genome'), 'Gene GFF3 (ANNEVO)')
        self.assertEqual(_download_file_label('species.gff3.gz', 'genome'), 'Gene GFF3')
        self.assertEqual(_download_file_label('species.ANNEVO.gtf.gz', 'genome'), 'Gene GTF (ANNEVO)')
        self.assertEqual(_download_file_label('species.pep.fasta.gz', 'genome'), 'Protein FASTA')
        self.assertEqual(_download_file_label('species.ANNEVO.CDS.fasta.gz', 'genome'), 'CDS FASTA (ANNEVO)')
        self.assertEqual(_download_file_label('species.emapper.annotations.gz', 'genome'), 'EGGNOG Annotation')
        self.assertEqual(_download_file_label('species.pep.fasta.tsv.gz', 'genome'), 'InterProScan Annotation')
        self.assertEqual(_download_file_label('species-families.fa.gz', 'genome'), 'Repeat FAMILY')
        self.assertEqual(_download_file_label('species.genomic.fasta.out.gz', 'genome'), 'Repeat LOCATION')
        self.assertEqual(_download_file_label('species.genomic.fasta.tbl', 'genome'), 'Repeat SUMMARY')

    def test_download_file_type_order(self):
        order = _download_file_type_order('genome')
        self.assertTrue(order.index('Genome FASTA') < order.index('Gene GFF3'))
        self.assertTrue(order.index('Repeat FAMILY') < order.index('Repeat LOCATION'))
        self.assertTrue(order.index('Repeat LOCATION') < order.index('Repeat SUMMARY'))

    def test_symbiont_labels_match_detail_page_types(self):
        self.assertEqual(_download_file_label('species.fna', 'symbiont'), 'Genome')
        self.assertEqual(_download_file_label('species.fna.gz', 'symbiont'), 'Genome')
        self.assertEqual(_download_file_label('species.faa', 'symbiont'), 'Protein')
        self.assertEqual(_download_file_label('species.ffn', 'symbiont'), 'Gene')
        self.assertEqual(_download_file_label('species.gff3', 'symbiont'), 'GFF3')

    def test_mitogenome_labels_match_detail_page_types(self):
        self.assertEqual(_download_file_label('mitoz_species.mt.fasta.cds.fasta', 'mitogenome'), 'CDS')
        self.assertEqual(_download_file_label('mitoz_species.mt.fasta.gene.fasta', 'mitogenome'), 'Gene')
        self.assertEqual(_download_file_label('mitoz_species.mt.fasta.cds_translation.fasta', 'mitogenome'), 'Protein')
        self.assertEqual(_download_file_label('mitoz_species.mt.fasta.rrna.fasta', 'mitogenome'), 'rRNA')
        self.assertEqual(_download_file_label('mitoz_species.mt.fasta.trna.fasta', 'mitogenome'), 'tRNA')
        self.assertEqual(_download_file_label('mitoz_species.mt.fasta.summary.2.txt', 'mitogenome'), 'Location')


class NCBIOmicsCacheTests(TestCase):
    def setUp(self):
        self.taxonomy = Taxonomy.objects.create(
            ncbi_taxonomy_id='123',
            species='Testus specius',
        )

    @patch('api.views._save_cached_ncbi_omics_counts')
    @patch('api.views._ncbi_sra_run_accessions')
    def test_build_ncbi_omics_cache_records_dna_and_rna_counts(self, mock_run_accessions, mock_save_cache):
        def side_effect(term, api_key=None):
            if 'biomol_dna' in term or 'genomic[LibrarySource]' in term:
                return ['SRR000001', 'SRR000002']
            if 'biomol_rna' in term or 'RNA-Seq[LibraryStrategy]' in term:
                return ['SRR000010']
            return []

        mock_run_accessions.side_effect = side_effect

        data = _build_ncbi_omics_cache()
        self.assertIn('123', data['counts'])
        self.assertEqual(data['counts']['123']['dna_count'], 2)
        self.assertEqual(data['counts']['123']['rna_count'], 1)
        self.assertEqual(data['dna_species_count'], 1)
        self.assertEqual(data['rna_species_count'], 1)

    @patch('api.views.urllib.request.urlopen')
    def test_ncbi_sra_ids_uses_stable_pagination(self, mock_urlopen):
        class Response:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc_value, traceback):
                return False

            def read(self):
                return b'{"esearchresult":{"count":"1","idlist":["7591051"]}}'

        mock_urlopen.return_value = Response()

        ids = _ncbi_sra_ids('txid60939[Organism:exp]')

        self.assertEqual(ids, ['7591051'])
        self.assertIn('retmax=500', mock_urlopen.call_args.args[0].full_url)

    @patch('api.views._cached_dna_counts')
    def test_dna_species_count_prefers_local_taxid_cache(self, mock_cached_dna):
        mock_cached_dna.return_value = ({'123': 2, '456': 1}, {})

        self.assertEqual(_get_dna_species_count({'dna_species_count': 99}), 2)

    @patch('api.views._cached_dna_counts')
    def test_dna_species_count_falls_back_to_legacy_cache(self, mock_cached_dna):
        mock_cached_dna.return_value = ({}, {})

        self.assertEqual(_get_dna_species_count({'dna_species_count': 7}), 7)


    @patch('api.views._get_ncbi_omics_cache')
    @patch('api.views._load_cached_transcriptome_rows')
    @patch('api.views._cached_dna_counts')
    def test_taxonomy_page_shows_dna_and_cached_rna_counts(self, mock_cached_dna, mock_load_transcriptome, mock_get_ncbi_cache):
        mock_cached_dna.return_value = ({'123': 4}, {})
        mock_load_transcriptome.return_value = ([{
            'run': 'SRR000010',
            'species': 'Testus specius',
            'library_source': 'TRANSCRIPTOMIC',
            'library_strategy': 'RNA-Seq',
        }], '2026-08-07T00:00:00Z')
        mock_get_ncbi_cache.return_value = {
            'counts': {
                '123': {
                    'dna_runs': ['SRR000001', 'SRR000002', 'SRR000003'],
                    'rna_runs': [],
                    'dna_count': 0,
                    'rna_count': 0,
                }
            },
            'dna_species_count': 1,
            'rna_species_count': 0,
        }
        response = self.client.get('/taxonomy/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'SRA (DNA)')
        self.assertContains(response, '<span class="badge bg-primary">4</span>', html=False)
        self.assertContains(response, 'SRA (RNA)')
        self.assertContains(response, '<span class="badge bg-success">1</span>', html=False)
        self.assertContains(response, 'https://www.ncbi.nlm.nih.gov/sra/?term=')


class MitogenomeStatisticsTests(TestCase):
    def test_mitogenome_species_count_uses_distinct_taxonomy_species(self):
        first_taxonomy = Taxonomy.objects.create(species='Species one')
        second_taxonomy = Taxonomy.objects.create(species='Species two')
        Mitogenome.objects.create(taxonomy=first_taxonomy)
        Mitogenome.objects.create(taxonomy=first_taxonomy)
        Mitogenome.objects.create(taxonomy=second_taxonomy)

        self.assertEqual(_mitogenome_species_count(), 2)


class SpeciesDetailViewTests(TestCase):
    def setUp(self):
        self.taxonomy = Taxonomy.objects.create(
            ncbi_taxonomy_id='123',
            species='Testus specius',
        )

    def test_species_detail_does_not_show_sra_run_counts(self):
        response = self.client.get(f'/species/{self.taxonomy.ncbi_taxonomy_id}/')
        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, 'SRA DNA runs')
        self.assertNotContains(response, 'SRA RNA runs')


class TranscriptomeViewTests(TestCase):
    @patch('api.views._load_cached_transcriptome_rows')
    def test_transcriptome_page_loads_rows_through_data_api(self, mock_load_rows):
        mock_load_rows.return_value = ([{
            'run': 'SRR000001',
            'study': 'SRX000001',
            'species': 'Testus specius',
            'tissue': 'whole organism',
            'development_stage': 'adult',
            'platform': 'ILLUMINA',
            'library_layout': 'PAIRED',
            'spots': '1 M',
            'bases': '2 M',
            'bioproject': 'PRJ000001',
            'species_url': None,
            'class_name': 'Clitellata',
            'order': 'Hirudinida',
            'family': 'Hirudinidae',
            'updated_clade1': 'Clitellata',
            'updated_clade2': '',
            'updated_clade3': '',
            'updated_clade4': '',
            'updated_clade5': '',
            'updated_clade6': 'Annelida',
        }], '2026-08-07T00:00:00Z')

        response = self.client.get('/transcriptome/')

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, '/api/transcriptome/data/')

        data_response = self.client.get('/api/transcriptome/data/', {
            'draw': '1',
            'start': '0',
            'length': '10',
        })
        self.assertEqual(data_response.status_code, 200)
        self.assertEqual(data_response.json()['data'][0][1], 'SRR000001')
        self.assertIn('Testus specius', data_response.json()['data'][0])

        post_response = self.client.post('/api/transcriptome/data/', {
            'draw': '2',
            'start': '0',
            'length': '10',
        })
        self.assertEqual(post_response.status_code, 200)
        self.assertEqual(post_response.json()['draw'], 2)

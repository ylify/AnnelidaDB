from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0021_add_eggnog_query_index'),
    ]

    operations = [
        migrations.CreateModel(
            name='InterproScanAnnotation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('query', models.CharField(db_index=True, max_length=200, verbose_name='Protein ID')),
                ('md5', models.CharField(blank=True, max_length=64, null=True, verbose_name='MD5')),
                ('length', models.IntegerField(blank=True, null=True, verbose_name='Length')),
                ('analysis', models.CharField(blank=True, max_length=100, null=True, verbose_name='Analysis')),
                ('signature_accession', models.CharField(blank=True, max_length=100, null=True, verbose_name='Signature Accession')),
                ('signature_description', models.TextField(blank=True, null=True, verbose_name='Signature Description')),
                ('start', models.IntegerField(blank=True, null=True, verbose_name='Start')),
                ('end', models.IntegerField(blank=True, null=True, verbose_name='End')),
                ('score', models.FloatField(blank=True, null=True, verbose_name='Score')),
                ('status', models.CharField(blank=True, max_length=50, null=True, verbose_name='Status')),
                ('date', models.CharField(blank=True, max_length=100, null=True, verbose_name='Date')),
                ('ipr_accession', models.CharField(blank=True, max_length=100, null=True, verbose_name='IPR Accession')),
                ('ipr_description', models.TextField(blank=True, null=True, verbose_name='IPR Description')),
                ('gos', models.TextField(blank=True, null=True, verbose_name='GO')),
                ('pathway', models.TextField(blank=True, null=True, verbose_name='Pathway')),
                ('taxonomy', models.ForeignKey(on_delete=models.deletion.CASCADE, related_name='interpro_annotations', to='api.taxonomy')),
            ],
        ),
    ]

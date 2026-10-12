from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0019_submissionrequest_uploaded_files'),
    ]

    operations = [
        migrations.CreateModel(
            name='SequenceFileIndex',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('seq_type', models.CharField(choices=[('protein', 'Protein'), ('cds', 'CDS')], max_length=10, verbose_name='Sequence Type')),
                ('seq_id', models.CharField(db_index=True, max_length=200, verbose_name='Sequence ID')),
                ('file_name', models.CharField(max_length=300, verbose_name='Source FASTA Filename')),
                ('file_path', models.CharField(max_length=500, verbose_name='Indexed FASTA Path')),
                ('taxonomy', models.ForeignKey(on_delete=models.deletion.CASCADE, related_name='sequence_file_indexes', to='api.taxonomy')),
            ],
            options={
                'unique_together': {('taxonomy', 'seq_type', 'seq_id')},
            },
        ),
        migrations.AddIndex(
            model_name='sequencefileindex',
            index=models.Index(fields=['seq_type', 'seq_id'], name='api_sequence_type_seq_id_idx'),
        ),
    ]

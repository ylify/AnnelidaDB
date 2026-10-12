from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0020_sequencefileindex'),
    ]

    operations = [
        migrations.AlterField(
            model_name='eggnogannotation',
            name='query',
            field=models.CharField(max_length=100, verbose_name='Query ID', db_index=True),
        ),
    ]

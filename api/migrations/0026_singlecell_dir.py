from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0025_symbiont_sub'),
    ]

    operations = [
        migrations.AddField(
            model_name='singlecell',
            name='dir',
            field=models.CharField(blank=True, max_length=250, null=True, verbose_name='Data directory'),
        ),
    ]
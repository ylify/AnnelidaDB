from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0024_taxonomy_updated_clade1_taxonomy_updated_clade2_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='symbiont',
            name='sub',
            field=models.CharField(blank=True, max_length=150, null=True, verbose_name='Sub'),
        ),
    ]
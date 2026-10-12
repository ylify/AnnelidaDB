from django.db import migrations, models
import uuid


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0015_submissionrequest'),
    ]

    operations = [
        migrations.AddField(
            model_name='submissionrequest',
            name='access_token',
            field=models.CharField(
                default=uuid.uuid4,
                editable=False,
                max_length=64,
                unique=True,
                verbose_name='Access Token',
            ),
            preserve_default=False,
        ),
    ]

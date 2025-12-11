# Generated manually for API preferences field

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemtools', '0003_chemjob_progress_message_chemjob_progress_percent'),
    ]

    operations = [
        migrations.AddField(
            model_name='chemjob',
            name='api_preferences',
            field=models.JSONField(blank=True, default=dict),
        ),
    ]

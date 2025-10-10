from __future__ import annotations

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("simulator", "0001_initial"),
    ]

    operations = [
        migrations.AddField(
            model_name="simulationattempt",
            name="parameters",
            field=models.JSONField(blank=True, default=dict),
        ),
        migrations.AddField(
            model_name="simulationattempt",
            name="results",
            field=models.JSONField(blank=True, default=dict),
        ),
        migrations.AddField(
            model_name="simulationattempt",
            name="results_summary",
            field=models.JSONField(blank=True, default=dict),
        ),
    ]

from __future__ import annotations

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("simulator", "0002_simulationattempt_parameters"),
    ]

    operations = [
        migrations.AddField(
            model_name="simulationattempt",
            name="artifacts",
            field=models.JSONField(blank=True, default=dict),
        ),
        migrations.AddField(
            model_name="simulationattempt",
            name="cohort_size",
            field=models.PositiveIntegerField(default=1),
        ),
        migrations.AddField(
            model_name="simulationattempt",
            name="seed",
            field=models.IntegerField(blank=True, null=True),
        ),
    ]

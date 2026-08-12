from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        (
            "twin_engine",
            "0002_observationresidual_stage_longitudinallabresult_and_more",
        ),
    ]

    operations = [
        migrations.AddField(
            model_name="patienttwinstate",
            name="lineage",
            field=models.JSONField(
                blank=True,
                default=dict,
            ),
        ),
    ]

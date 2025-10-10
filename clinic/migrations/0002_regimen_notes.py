from __future__ import annotations

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("clinic", "0001_initial"),
    ]

    operations = [
        migrations.AddField(
            model_name="regimen",
            name="notes",
            field=models.TextField(blank=True),
        ),
    ]

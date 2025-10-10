from __future__ import annotations

from rest_framework import serializers

from . import models


class ChemJobSerializer(serializers.ModelSerializer):
    out_html_url = serializers.SerializerMethodField()
    out_csv_url = serializers.SerializerMethodField()
    status_label = serializers.SerializerMethodField()
    thumbnail_url = serializers.SerializerMethodField()

    class Meta:
        model = models.ChemJob
        fields = [
            "id",
            "created",
            "kind",
            "input_a",
            "input_b",
            "out_html_url",
            "out_csv_url",
            "status_label",
            "thumbnail_url",
        ]

    def get_out_html_url(self, obj: models.ChemJob) -> str | None:
        return obj.out_html.url if obj.out_html else None

    def get_out_csv_url(self, obj: models.ChemJob) -> str | None:
        return obj.out_csv.url if obj.out_csv else None

    def get_status_label(self, obj: models.ChemJob) -> str:
        text, _variant = obj.status_label()
        return text

    def get_thumbnail_url(self, obj: models.ChemJob) -> str | None:
        return obj.thumbnail_url()

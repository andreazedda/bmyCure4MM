from __future__ import annotations

from rest_framework import mixins, viewsets

from . import models, serializers


class ChemJobViewSet(mixins.ListModelMixin, mixins.RetrieveModelMixin, viewsets.GenericViewSet):
    """Read-only viewset exposing ChemJob metadata."""

    queryset = models.ChemJob.objects.select_related("user").all()
    serializer_class = serializers.ChemJobSerializer
    ordering = ["-created"]
    http_method_names = ["get", "head", "options"]

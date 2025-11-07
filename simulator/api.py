from __future__ import annotations

from rest_framework.decorators import api_view
from rest_framework.response import Response

from . import explain
from .pharmaco import registry as pharmaco_registry


@api_view(["GET"])
def glossary(_request):
    """Return bilingual KPI glossary."""
    return Response(explain.KPI)


@api_view(["GET"])
def drugs(_request):
    """Return PK/PD preset registry for UI reuse."""
    return Response(pharmaco_registry.list_profiles())

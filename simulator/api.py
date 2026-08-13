from __future__ import annotations

from rest_framework.decorators import api_view
from rest_framework.response import Response

from mmportal.governance import EpistemicLabel, governance_metadata

from . import explain
from .pharmaco import registry as pharmaco_registry


@api_view(["GET"])
def glossary(_request):
    """Return bilingual KPI glossary."""
    return Response({
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.DERIVED,
            output_kind="model_kpi_glossary",
        ),
        "kpi": explain.KPI,
    })


@api_view(["GET"])
def drugs(_request):
    """Return hypothetical PK/PD model-input profiles for UI reuse."""
    return Response({
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.HYPOTHETICAL,
            output_kind="pkpd_model_input_profiles",
        ),
        "profiles": pharmaco_registry.list_profiles(),
    })

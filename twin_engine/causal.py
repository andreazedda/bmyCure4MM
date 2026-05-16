from __future__ import annotations

from typing import Any

from django.core.exceptions import ValidationError

from .models import CausalAssumptionSet


def build_causal_assumption_set(
    *,
    name: str,
    graph_definition: dict[str, Any],
    variables: list[str],
    intervention: dict[str, Any],
    outcome: dict[str, Any],
    adjustment_set: list[str] | None = None,
    assumptions: dict[str, Any] | None = None,
    identification_status: str | None = None,
    notes: str = "",
    patient=None,
    created_by=None,
) -> CausalAssumptionSet:
    validate_causal_graph_definition(graph_definition)
    assumption_set = CausalAssumptionSet.objects.create(
        patient=patient,
        name=name,
        graph_definition=graph_definition,
        variables=variables,
        intervention=intervention,
        outcome=outcome,
        adjustment_set=adjustment_set or [],
        assumptions=assumptions or {},
        identification_status=identification_status or CausalAssumptionSet.IDENT_MECHANISTIC_ONLY,
        notes=notes,
        created_by=created_by,
    )
    return assumption_set


def validate_causal_graph_definition(graph_definition: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(graph_definition, dict):
        raise ValidationError("graph_definition must be a JSON object.")
    nodes = graph_definition.get("nodes", [])
    edges = graph_definition.get("edges", [])
    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise ValidationError("graph_definition requires list-valued nodes and edges.")
    for edge in edges:
        if not isinstance(edge, (list, tuple)) or len(edge) != 2:
            raise ValidationError("Each edge must contain exactly two nodes.")
    return graph_definition


def classify_estimand(
    *,
    graph_definition: dict[str, Any] | None = None,
    intervention: dict[str, Any] | None = None,
    outcome: dict[str, Any] | None = None,
    adjustment_set: list[str] | None = None,
    identification_status: str | None = None,
) -> str:
    if identification_status == CausalAssumptionSet.IDENT_MECHANISTIC_ONLY:
        return "mechanistic_simulation_only"
    if not graph_definition or not intervention or not outcome:
        return "insufficient_data"
    if identification_status == CausalAssumptionSet.IDENT_IDENTIFIED_UNDER_ASSUMPTIONS:
        return "causal_estimand_identified_under_assumptions"
    if identification_status in {
        CausalAssumptionSet.IDENT_NOT_IDENTIFIED,
        CausalAssumptionSet.IDENT_PARTIALLY_IDENTIFIED,
    }:
        return "causal_estimand_defined_not_identified"
    if adjustment_set:
        return "causal_estimand_defined_not_identified"
    return "insufficient_data"


def distinguish_mechanistic_counterfactual_vs_causal_estimand(
    *,
    graph_definition: dict[str, Any] | None = None,
    intervention: dict[str, Any] | None = None,
    outcome: dict[str, Any] | None = None,
    adjustment_set: list[str] | None = None,
    identification_status: str | None = None,
) -> dict[str, Any]:
    classification = classify_estimand(
        graph_definition=graph_definition,
        intervention=intervention,
        outcome=outcome,
        adjustment_set=adjustment_set,
        identification_status=identification_status,
    )
    return {
        "classification": classification,
        "mechanistic_label": "mechanistic model counterfactual",
        "causal_label": "causal inference estimand",
        "exploratory_label": "unvalidated exploratory hypothesis",
        "warning": (
            "Single-patient simulation does not establish causal proof."
            if classification != "causal_estimand_identified_under_assumptions"
            else "Interpretation remains conditional on the declared assumptions."
        ),
    }

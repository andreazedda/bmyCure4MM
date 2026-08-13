"""
API endpoint for exposing the decision algorithm in machine-readable format.

This makes the algorithm fully transparent and auditable.
"""

from __future__ import annotations

import json

from django.http import JsonResponse, HttpRequest
from django.views.decorators.http import require_GET

from .decision_algorithm import get_algorithm, evaluate_metrics, get_applicable_rules
from mmportal.governance import EpistemicLabel, governance_metadata


@require_GET
def decision_algorithm_api(request: HttpRequest) -> JsonResponse:
    """
    Return the complete decision algorithm as JSON.
    
    Query parameters:
    - format: 'full' (default) or 'summary'
    - lang: 'en' (default) or 'it'
    - evaluate: if present with tumor_reduction, healthy_loss, time_to_recurrence,
                returns applicable rules for those values
    
    Examples:
        GET /simulator/api/decision-algorithm/
        GET /simulator/api/decision-algorithm/?format=summary
        GET /simulator/api/decision-algorithm/?evaluate=1&tumor_reduction=0.45&healthy_loss=0.25
    """
    fmt = request.GET.get("format", "full")
    lang = request.GET.get("lang", "en")
    
    algorithm = get_algorithm()
    
    # If evaluation requested, compute applicable rules
    if request.GET.get("evaluate"):
        try:
            tumor_reduction = float(request.GET.get("tumor_reduction", ""))
        except (ValueError, TypeError):
            tumor_reduction = None
        
        try:
            healthy_loss = float(request.GET.get("healthy_loss", ""))
        except (ValueError, TypeError):
            healthy_loss = None
        
        try:
            time_to_recurrence = float(request.GET.get("time_to_recurrence", ""))
        except (ValueError, TypeError):
            time_to_recurrence = None
        
        try:
            time_horizon = int(request.GET.get("time_horizon", "168"))
        except (ValueError, TypeError):
            time_horizon = 168
        
        # Evaluate metrics
        evaluation = evaluate_metrics(tumor_reduction, healthy_loss, time_to_recurrence)
        
        # Get applicable rules
        applicable_rules = get_applicable_rules(
            tumor_reduction, healthy_loss, time_to_recurrence, time_horizon
        )
        
        return JsonResponse({
            "governance": governance_metadata(
                epistemic_label=EpistemicLabel.HEURISTIC,
                output_kind="model_relative_threshold_evaluation",
            ),
            "input": {
                "tumor_reduction": tumor_reduction,
                "healthy_loss": healthy_loss,
                "time_to_recurrence": time_to_recurrence,
                "time_horizon": time_horizon,
            },
            "evaluation": evaluation,
            "applicable_rules": applicable_rules,
            "algorithm_version": algorithm.get("version"),
        })
    
    # Summary format - just thresholds and rule names
    if fmt == "summary":
        summary = {
            "governance": governance_metadata(
                epistemic_label=EpistemicLabel.HEURISTIC,
                output_kind="model_relative_threshold_catalog",
            ),
            "version": algorithm.get("version"),
            "last_updated": algorithm.get("last_updated"),
            "thresholds": algorithm.get("thresholds"),
            "rule_count": len(algorithm.get("decision_rules", [])),
            "rules_summary": [
                {
                    "id": r.get("id"),
                    "name_en": r.get("name_en"),
                    "name_it": r.get("name_it"),
                    "priority": r.get("priority"),
                }
                for r in algorithm.get("decision_rules", [])
            ],
            "risk_stratification_stages": list(algorithm.get("risk_stratification", {}).keys()),
            "disclaimer": algorithm.get("disclaimer", {}).get(lang, ""),
        }
        return JsonResponse(summary)
    
    # Full algorithm
    return JsonResponse(algorithm)

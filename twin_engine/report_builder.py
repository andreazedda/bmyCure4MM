from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from django.conf import settings

from .provenance import hash_json
from .validators import validate_no_direct_identifier_in_artifact_payload


def build_counterfactual_report_payload(
    counterfactual_run,
    *,
    baseline_summary,
    alternative_summary,
    comparison_metrics,
    metadata=None,
    toxicity_constraints=None,
    warnings=None,
) -> dict[str, Any]:
    patient_token = hash_json({"patient_id": counterfactual_run.patient_id})[:12]
    payload = {
        "label": "research simulation",
        "patient_reference": patient_token,
        "counterfactual_type": "mechanistic model counterfactual",
        "causal_interpretation": "unvalidated exploratory hypothesis",
        "counterfactual_run_id": counterfactual_run.id,
        "base_twin_state_id": counterfactual_run.base_twin_state_id,
        "status": counterfactual_run.status,
        "baseline_summary": baseline_summary,
        "alternative_summary": alternative_summary,
        "baseline_exposure_summary": (baseline_summary or {}).get("exposure_summary") or {},
        "alternative_exposure_summary": (alternative_summary or {}).get("exposure_summary") or {},
        "baseline_toxicity_dynamics": (baseline_summary or {}).get("toxicity_dynamics") or {},
        "alternative_toxicity_dynamics": (alternative_summary or {}).get("toxicity_dynamics") or {},
        "baseline_predicted_biomarkers": (baseline_summary or {}).get("predicted_biomarkers"),
        "predicted_biomarkers": (alternative_summary or {}).get("predicted_biomarkers"),
        "comparison_metrics": comparison_metrics,
        "schedule_comparison": (comparison_metrics or {}).get("schedule_comparison") or {},
        "toxicity_constraints": toxicity_constraints or {},
        "warnings": list(warnings or []),
        "provenance": metadata or {},
    }
    validate_no_direct_identifier_in_artifact_payload(payload)
    return payload


def write_json_artifact(prefix: str, payload: dict[str, Any], *, patient_id: int, run_id: int, folder_name: str) -> tuple[str, str]:
    validate_no_direct_identifier_in_artifact_payload(payload)

    media_root = Path(settings.MEDIA_ROOT)
    artifact_dir = media_root / folder_name
    artifact_dir.mkdir(parents=True, exist_ok=True)

    pseudonym = hash_json({"patient_id": patient_id, "run_id": run_id})[:12]
    filename = f"{prefix}_{pseudonym}_{run_id}.json"
    path = artifact_dir / filename
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    media_url = settings.MEDIA_URL.rstrip("/")
    return f"{media_url}/{folder_name}/{filename}", str(path)

from __future__ import annotations

from typing import Any

from django.core.exceptions import ValidationError

from simulator.permissions import is_editor


BANNED_IDENTIFIER_KEYS = {
    "first_name",
    "last_name",
    "full_name",
    "mrn",
    "medical_record_number",
    "patient_name",
    "notes",
    "note",
    "schedule_notes",
    "birth_date",
}

RESEARCH_INTERVENTION_KEYS = {
    "alternative_regimen_id",
    "drug_doses",
    "start_day",
    "duration_days",
    "schedule",
    "parameter_overrides",
    "random_seed",
    "causal_graph",
    "adjustment_set",
    "identification_status",
}


def validate_no_direct_identifier_in_artifact_payload(payload: Any) -> Any:
    offending_paths: list[str] = []

    def _walk(value: Any, path: list[str]) -> None:
        if isinstance(value, dict):
            for key, nested in value.items():
                key_text = str(key).lower()
                next_path = path + [str(key)]
                if key_text in BANNED_IDENTIFIER_KEYS:
                    offending_paths.append(".".join(next_path))
                _walk(nested, next_path)
            return
        if isinstance(value, list):
            for index, nested in enumerate(value):
                _walk(nested, path + [str(index)])

    _walk(payload, [])
    if offending_paths:
        raise ValidationError(
            "Artifact payload contains direct identifiers: " + ", ".join(sorted(offending_paths))
        )
    return payload


def validate_patient_access(user, patient) -> bool:
    if not user or not getattr(user, "is_authenticated", False):
        return False
    if getattr(user, "is_staff", False) or is_editor(user):
        return True
    if getattr(patient, "owner_id", None) == user.id:
        return True
    mrn = str(getattr(patient, "mrn", "") or "").upper()
    return getattr(patient, "owner_id", None) is None and mrn.startswith("DEMO")


def validate_research_run_inputs(patient, base_twin_state, intervention_definition: dict[str, Any], horizon_days: int) -> dict[str, Any]:
    if base_twin_state.patient_id != patient.id:
        raise ValidationError("Base twin state does not belong to the requested patient.")
    if horizon_days <= 0:
        raise ValidationError("horizon_days must be greater than zero.")
    if not isinstance(intervention_definition, dict):
        raise ValidationError("intervention_definition must be a JSON object.")

    unexpected_keys = sorted(set(intervention_definition) - RESEARCH_INTERVENTION_KEYS)
    if unexpected_keys:
        raise ValidationError(
            "Unexpected intervention_definition keys: " + ", ".join(unexpected_keys)
        )

    validate_no_direct_identifier_in_artifact_payload(intervention_definition)
    return intervention_definition


def validate_assessment_minimum_fields(assessment) -> dict[str, list[str]]:
    available: list[str] = []
    missing: list[str] = []

    if getattr(assessment, "date", None):
        available.append("date")
    else:
        missing.append("date")

    biomarker_fields = [
        "m_protein_g_dl",
        "flc_ratio",
        "hemoglobin_g_dl",
    ]
    for field_name in biomarker_fields:
        value = getattr(assessment, field_name, None)
        if value is None or value == "":
            missing.append(field_name)
        else:
            available.append(field_name)

    if "date" not in available:
        raise ValidationError("Assessment must include a date.")

    if not any(field in available for field in biomarker_fields):
        raise ValidationError(
            "Assessment must include at least one of m_protein_g_dl, flc_ratio, or hemoglobin_g_dl."
        )

    return {
        "available": available,
        "missing": missing,
    }

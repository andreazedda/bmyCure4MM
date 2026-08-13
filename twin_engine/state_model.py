from __future__ import annotations

from django.db import transaction

from simulator.twin import build_patient_twin

from .input_contract import build_twin_lineage
from .models import PatientTwinState
from .observation_model import default_observation_parameters
from .provenance import CURRENT_MODEL_VERSION, collect_twin_config_hash, record_simulation_metadata
from .validators import validate_assessment_minimum_fields


def initialize_from_assessment(assessment, user=None) -> PatientTwinState:
    validate_assessment_minimum_fields(assessment)

    preview_payload = preview_state_from_assessment(assessment)
    twin_payload = preview_payload["twin_payload"]
    state_vector = preview_payload["state_vector"]
    parameters = preview_payload["parameters"]
    lineage = build_twin_lineage(
        patient=assessment.patient,
        assessments=[assessment],
        therapies=[],
        purpose="initialization",
    )

    with transaction.atomic():
        state = PatientTwinState.objects.create(
            patient=assessment.patient,
            assessment=assessment,
            state_date=assessment.date,
            is_current=False,
            state_vector=state_vector,
            parameters=parameters,
            parameter_uncertainty={},
            risk_score=twin_payload.get("risk_score"),
            method=PatientTwinState.METHOD_INITIAL_RISK_MAPPING,
            model_version=CURRENT_MODEL_VERSION,
            config_hash=collect_twin_config_hash(),
            lineage=lineage,
            created_by=user,
        )
        state.source_assessments.add(assessment)
        set_current_state(state)
        record_simulation_metadata(
            twin_state=state,
            model_id="patient_twin_state_model",
            solver_name="initial_state_mapping",
            input_payload={
                "assessment_date": assessment.date.isoformat(),
                "computational_input_sha256": lineage["computational_input"]["sha256"],
            },
            solver_parameters={"mapping": "simulator.twin.build_patient_twin"},
            output_payload={
                "state_vector": state_vector,
                "parameters": parameters,
                "lineage_sha256": lineage["computational_input"]["sha256"],
            },
        )
    return state


def preview_state_from_assessment(assessment) -> dict[str, object]:
    twin_payload = build_patient_twin(assessment)
    state_vector = _infer_state_vector_from_assessment(assessment)
    parameters = {
        **twin_payload,
        "observation": default_observation_parameters(
            assessment,
            tumor_cells=state_vector["tumor_cells"],
            healthy_cells=state_vector["healthy_cells"],
        ),
        "assumptions": {
            "label": "research simulation",
            "counterfactual_type": "mechanistic model counterfactual",
            "validation_status": "unvalidated exploratory hypothesis",
        },
    }
    return {
        "twin_payload": twin_payload,
        "state_vector": state_vector,
        "parameters": parameters,
        "config_hash": collect_twin_config_hash(),
        "model_version": CURRENT_MODEL_VERSION,
    }


def get_current_twin_state(patient) -> PatientTwinState | None:
    return patient.twin_states.filter(is_current=True).first()


def set_current_state(state: PatientTwinState) -> PatientTwinState:
    with transaction.atomic():
        PatientTwinState.objects.filter(patient=state.patient, is_current=True).exclude(pk=state.pk).update(is_current=False)
        if not state.is_current:
            state.is_current = True
            state.save(update_fields=["is_current"])
    return state


def serialize_state(state: PatientTwinState) -> dict[str, object]:
    return {
        "id": state.id,
        "patient_id": state.patient_id,
        "assessment_id": state.assessment_id,
        "state_date": state.state_date.isoformat(),
        "is_current": state.is_current,
        "state_vector": state.state_vector,
        "parameters": state.parameters,
        "parameter_uncertainty": state.parameter_uncertainty,
        "risk_score": state.risk_score,
        "method": state.method,
        "model_version": state.model_version,
        "config_hash": state.config_hash,
        "lineage": state.lineage,
        "source_assessment_ids": list(state.source_assessments.values_list("id", flat=True)),
        "created_at": state.created_at.isoformat(),
        "created_by_id": state.created_by_id,
    }


def _infer_state_vector_from_assessment(assessment) -> dict[str, float]:
    m_protein = _float_or_default(getattr(assessment, "m_protein_g_dl", None), 0.0)
    flc_ratio = _float_or_default(getattr(assessment, "flc_ratio", None), 1.0)
    hemoglobin = _float_or_default(getattr(assessment, "hemoglobin_g_dl", None), 12.0)

    tumor_cells = max(1.0e8, 5.0e8 + (m_protein * 4.0e8) + max(flc_ratio - 1.0, 0.0) * 5.0e7)
    healthy_fraction = min(max(hemoglobin / 13.0, 0.2), 1.2)
    healthy_cells = max(1.0e10, 5.0e11 * healthy_fraction)

    return {
        "tumor_cells": float(tumor_cells),
        "healthy_cells": float(healthy_cells),
        "m_protein_g_dl": float(m_protein),
        "flc_ratio": float(flc_ratio),
        "hemoglobin_g_dl": float(hemoglobin),
    }


def _float_or_default(value, default: float) -> float:
    if value is None or value == "":
        return float(default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)

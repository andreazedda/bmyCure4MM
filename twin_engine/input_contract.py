from __future__ import annotations

import json
import math
from datetime import date, datetime
from decimal import Decimal
from typing import Any

from django.db.models import Q

from .models import (
    AdverseEvent,
    LongitudinalLabResult,
    TherapyInterruption,
)
from .provenance import (
    CURRENT_MODEL_VERSION,
    collect_drug_preset_hashes,
    collect_twin_config_hash,
    get_code_commit_hash,
    hash_json,
)


COMPUTATIONAL_INPUT_CONTRACT = "research-computational-input-v0.1"
TWIN_LINEAGE_CONTRACT = "research-twin-lineage-v0.1"
DATASET_BINDING_SCOPE = "structured_subset_only"


ASSESSMENT_COMPUTATIONAL_FIELDS = (
    "m_protein_g_dl",
    "flc_ratio",
    "hemoglobin_g_dl",
    "beta2m_mg_l",
    "ldH_u_l",
    "r_iss",
)


def _normalize(value: Any) -> Any:
    if isinstance(value, (date, datetime)):
        return value.isoformat()

    if isinstance(value, Decimal):
        if value == 0:
            return "0"
        return format(value.normalize(), "f")

    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Non-finite float is not allowed in computational input")
        return float(value)

    if isinstance(value, dict):
        return {
            str(key): _normalize(value[key])
            for key in sorted(value, key=str)
        }

    if isinstance(value, (list, tuple)):
        return [_normalize(item) for item in value]

    return value


def _canonical_sort_key(payload: dict[str, Any]) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def serialize_assessment_input(assessment) -> dict[str, Any]:
    payload = {
        "date": assessment.date,
    }

    for field_name in ASSESSMENT_COMPUTATIONAL_FIELDS:
        payload[field_name] = getattr(
            assessment,
            field_name,
            None,
        )

    return _normalize(payload)


def serialize_start_state_input(state) -> dict[str, Any] | None:
    if state is None:
        return None

    return _normalize(
        {
            "state_date": state.state_date,
            "state_vector": state.state_vector or {},
            "parameters": state.parameters or {},
            "risk_score": state.risk_score,
            "method": state.method,
            "model_version": state.model_version,
            "config_hash": state.config_hash,
        }
    )


def serialize_therapy_input(therapy) -> dict[str, Any]:
    regimen = therapy.regimen

    days_on = sorted(
        {
            int(day)
            for day in (therapy.days_on or [])
        }
    )

    return _normalize(
        {
            "start_date": therapy.start_date,
            "end_date": therapy.end_date,
            "regimen": {
                "name": regimen.name,
                "components": regimen.components,
            },
            "cycle_length_days": therapy.cycle_length_days,
            "days_on": days_on,
            "doses": therapy.doses or {},
        }
    )


def serialize_interruption_input(interruption) -> dict[str, Any]:
    return _normalize(
        {
            "start_date": interruption.start_date,
            "end_date": interruption.end_date,
            "drug": str(
                interruption.drug or ""
            ).strip().lower(),
        }
    )


def build_computational_input_manifest(
    *,
    assessments,
    therapies,
    interruptions,
    horizon_start_date,
    horizon_end_date,
    purpose: str,
    start_state=None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    assessment_payloads = [
        serialize_assessment_input(item)
        for item in assessments
    ]

    therapy_payloads = [
        serialize_therapy_input(item)
        for item in therapies
    ]

    interruption_payloads = [
        serialize_interruption_input(item)
        for item in interruptions
    ]

    assessment_payloads.sort(
        key=_canonical_sort_key
    )
    therapy_payloads.sort(
        key=_canonical_sort_key
    )
    interruption_payloads.sort(
        key=_canonical_sort_key
    )

    return _normalize(
        {
            "contract": COMPUTATIONAL_INPUT_CONTRACT,
            "purpose": purpose,
            "horizon": {
                "start_date": horizon_start_date,
                "end_date": horizon_end_date,
            },
            "start_state": serialize_start_state_input(
                start_state
            ),
            "assessments": assessment_payloads,
            "therapies": therapy_payloads,
            "interruptions": interruption_payloads,
            "extra": extra or {},
        }
    )


def computational_input_sha256(
    manifest: dict[str, Any],
) -> str:
    return hash_json(manifest)


def _dataset_provenance_payload(
    record_type: str,
    row,
) -> dict[str, Any]:
    if record_type == "therapy_interruption":
        return dict(
            ((row.evidence or {}).get(
                "_dataset_provenance",
                {},
            ))
        )

    return dict(row.provenance or {})


def collect_structured_dataset_binding(
    patient,
) -> dict[str, Any]:
    rows: list[tuple[str, Any]] = []

    rows.extend(
        ("lab_result", row)
        for row in LongitudinalLabResult.objects.filter(
            patient=patient
        )
    )

    rows.extend(
        ("adverse_event", row)
        for row in AdverseEvent.objects.filter(
            patient=patient
        )
    )

    rows.extend(
        ("therapy_interruption", row)
        for row in TherapyInterruption.objects.filter(
            patient=patient
        )
    )

    binding_counts: dict[
        tuple[str, str, str],
        int,
    ] = {}

    record_counts = {
        "lab_result": 0,
        "adverse_event": 0,
        "therapy_interruption": 0,
    }

    bound_record_count = 0
    unbound_record_count = 0

    for record_type, row in rows:
        record_counts[record_type] += 1

        provenance = _dataset_provenance_payload(
            record_type,
            row,
        )

        identity = (
            str(provenance.get("dataset_id") or ""),
            str(provenance.get("dataset_version") or ""),
            str(
                provenance.get(
                    "canonical_dataset_sha256"
                )
                or ""
            ),
        )

        if all(identity):
            binding_counts[identity] = (
                binding_counts.get(identity, 0) + 1
            )
            bound_record_count += 1
        else:
            unbound_record_count += 1

    base = {
        "coverage_scope": DATASET_BINDING_SCOPE,
        "assessment_records_covered": False,
        "therapy_course_records_covered": False,
        "record_counts": record_counts,
        "bound_record_count": bound_record_count,
        "unbound_record_count": unbound_record_count,
    }

    if not binding_counts:
        return {
            **base,
            "status": "absent",
        }

    if len(binding_counts) > 1:
        bindings = [
            {
                "dataset_id": identity[0],
                "dataset_version": identity[1],
                "canonical_dataset_sha256": identity[2],
                "record_count": count,
            }
            for identity, count in sorted(
                binding_counts.items()
            )
        ]

        return {
            **base,
            "status": "mixed",
            "bindings": bindings,
        }

    identity, count = next(
        iter(binding_counts.items())
    )

    status = (
        "bound"
        if unbound_record_count == 0
        else "partially_bound"
    )

    return {
        **base,
        "status": status,
        "dataset_id": identity[0],
        "dataset_version": identity[1],
        "canonical_dataset_sha256": identity[2],
        "dataset_bound_record_count": count,
    }


def _interruptions_for_horizon(
    patient,
    start_date,
    end_date,
):
    return list(
        patient.therapy_interruptions.filter(
            start_date__lte=end_date,
        )
        .filter(
            Q(end_date__isnull=True)
            | Q(end_date__gte=start_date)
        )
        .order_by(
            "start_date",
            "drug",
        )
    )


def build_twin_lineage(
    *,
    patient,
    assessments,
    therapies,
    purpose: str,
    parent_state=None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    assessment_list = sorted(
        list(assessments),
        key=lambda item: (
            item.date,
            item.pk or 0,
        ),
    )

    therapy_list = sorted(
        list(therapies),
        key=lambda item: (
            item.start_date,
            item.end_date or item.start_date,
            item.regimen.name,
            item.pk or 0,
        ),
    )

    if not assessment_list:
        raise ValueError(
            "At least one assessment is required "
            "to build twin lineage"
        )

    if parent_state is None:
        horizon_start_date = assessment_list[0].date
    else:
        horizon_start_date = parent_state.state_date

    horizon_end_date = assessment_list[-1].date

    interruptions = (
        _interruptions_for_horizon(
            patient,
            horizon_start_date,
            horizon_end_date,
        )
        if purpose == "calibration"
        else []
    )

    manifest = build_computational_input_manifest(
        assessments=assessment_list,
        therapies=therapy_list,
        interruptions=interruptions,
        horizon_start_date=horizon_start_date,
        horizon_end_date=horizon_end_date,
        purpose=purpose,
        start_state=parent_state,
        extra=extra,
    )

    input_sha = computational_input_sha256(
        manifest
    )

    parent_input_sha = ""

    if parent_state is not None:
        parent_input_sha = str(
            (
                (parent_state.lineage or {})
                .get(
                    "computational_input",
                    {},
                )
                .get(
                    "sha256",
                    "",
                )
            )
        )

    lineage = {
        "contract": TWIN_LINEAGE_CONTRACT,
        "purpose": purpose,
        "runtime": {
            "model_version": CURRENT_MODEL_VERSION,
            "twin_config_hash": collect_twin_config_hash(),
            "drug_preset_hashes": (
                collect_drug_preset_hashes()
            ),
            "code_commit_hash": get_code_commit_hash(),
        },
        "dataset_binding": (
            collect_structured_dataset_binding(
                patient
            )
        ),
        "computational_input": {
            "contract": COMPUTATIONAL_INPUT_CONTRACT,
            "sha256": input_sha,
            "assessment_count": len(
                assessment_list
            ),
            "therapy_count": len(
                therapy_list
            ),
            "interruption_count": len(
                interruptions
            ),
            "horizon_start_date": (
                horizon_start_date.isoformat()
            ),
            "horizon_end_date": (
                horizon_end_date.isoformat()
            ),
            "contains_local_database_ids": False,
            "raw_manifest_persisted": False,
        },
    }

    if parent_state is not None:
        lineage["parent"] = {
            # Local linkage for navigating this database only.
            # It is deliberately excluded from the
            # computational-input hash.
            "local_twin_state_id": parent_state.pk,
            "parent_computational_input_sha256": (
                parent_input_sha
            ),
            "model_version": parent_state.model_version,
            "config_hash": parent_state.config_hash,
        }

    return lineage

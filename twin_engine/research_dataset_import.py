from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any

from django.db import transaction
from jsonschema import Draft202012Validator, FormatChecker

from clinic.models import Assessment
from twin_engine.models import (
    AdverseEvent,
    LongitudinalLabResult,
    TherapyInterruption,
)


SCHEMA_PATH = (
    Path(__file__).resolve().parent
    / "schemas"
    / "research_dataset_v0_1.schema.json"
)

DATASET_PROVENANCE_KEY = "_dataset_provenance"

LEGACY_EVIDENCE_PROVENANCE_KEYS = {
    "source",
    "extraction_status",
    "case_label",
}

BANNED_IDENTIFIER_KEYS = {
    "mrn",
    "medical_record_number",
    "patient_name",
    "first_name",
    "last_name",
    "birth_date",
    "date_of_birth",
    "dob",
}

SUPPORTED_RECORD_TYPES = {
    "lab_result",
    "adverse_event",
    "therapy_interruption",
}


class DatasetImportError(ValueError):
    pass


def canonical_dataset_bytes(payload: dict[str, Any]) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def canonical_dataset_sha256(payload: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_dataset_bytes(payload)).hexdigest()


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)

    return digest.hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise DatasetImportError(f"{label} does not exist") from exc
    except json.JSONDecodeError as exc:
        raise DatasetImportError(
            f"{label} is not valid JSON at line {exc.lineno}"
        ) from exc

    if not isinstance(payload, dict):
        raise DatasetImportError(f"{label} must contain a JSON object")

    return payload


def _walk_identifier_keys(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = str(key).strip().lower()

            if normalized in BANNED_IDENTIFIER_KEYS:
                raise DatasetImportError(
                    f"direct identifier key is forbidden at {path}.{key}"
                )

            _walk_identifier_keys(child, f"{path}.{key}")

    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_identifier_keys(child, f"{path}[{index}]")


def _validate_schema(dataset: dict[str, Any]) -> None:
    schema = _read_json(SCHEMA_PATH, "canonical dataset schema")

    Draft202012Validator.check_schema(schema)

    validator = Draft202012Validator(
        schema,
        format_checker=FormatChecker(),
    )

    errors = sorted(
        validator.iter_errors(dataset),
        key=lambda error: list(error.absolute_path),
    )

    if not errors:
        return

    first = errors[0]
    path = ".".join(str(part) for part in first.absolute_path) or "$"

    # Deliberately do not print the private offending value.
    raise DatasetImportError(
        "dataset schema validation failed: "
        f"path={path} validator={first.validator} "
        f"errors={len(errors)}"
    )


def _semantic_key(record: dict[str, Any]) -> tuple[Any, ...]:
    record_type = record["record_type"]
    identity = record["identity"]

    if record_type == "lab_result":
        return (
            record_type,
            identity.get("case_ref"),
            identity.get("date"),
            identity.get("analyte"),
        )

    if record_type == "adverse_event":
        return (
            record_type,
            identity.get("case_ref"),
            identity.get("date"),
            identity.get("event_type"),
        )

    if record_type == "therapy_interruption":
        return (
            record_type,
            identity.get("case_ref"),
            identity.get("drug"),
            identity.get("start_date"),
            identity.get("end_date"),
        )

    raise DatasetImportError(
        f"unsupported record_type: {record_type}"
    )


def _validate_record_semantics(dataset: dict[str, Any]) -> None:
    source_ids = set()
    record_ids = set()
    assertion_ids = set()
    semantic_keys = set()

    case_ref = dataset["case_ref"]

    for source in dataset["sources"]:
        source_id = source["source_id"]

        if source_id in source_ids:
            raise DatasetImportError(
                f"duplicate source_id: {source_id}"
            )

        source_ids.add(source_id)

    for record in dataset["records"]:
        record_id = record["record_id"]
        record_type = record["record_type"]
        identity = record["identity"]
        payload = record["payload"]

        if record_type not in SUPPORTED_RECORD_TYPES:
            raise DatasetImportError(
                "Dataset v0.1 importer currently supports only "
                "lab_result, adverse_event and therapy_interruption"
            )

        if record_id in record_ids:
            raise DatasetImportError(
                f"duplicate record_id: {record_id}"
            )

        record_ids.add(record_id)

        if identity.get("case_ref") != case_ref:
            raise DatasetImportError(
                f"record case_ref mismatch: {record_id}"
            )

        semantic_key = _semantic_key(record)

        if semantic_key in semantic_keys:
            raise DatasetImportError(
                f"duplicate semantic identity: {record_id}"
            )

        semantic_keys.add(semantic_key)

        if record_type == "lab_result":
            if payload.get("date") != identity.get("date"):
                raise DatasetImportError(
                    f"date identity/payload mismatch: {record_id}"
                )

            if payload.get("analyte") != identity.get("analyte"):
                raise DatasetImportError(
                    f"analyte identity/payload mismatch: {record_id}"
                )

        elif record_type == "adverse_event":
            if payload.get("date") != identity.get("date"):
                raise DatasetImportError(
                    f"date identity/payload mismatch: {record_id}"
                )

            if payload.get("event_type") != identity.get("event_type"):
                raise DatasetImportError(
                    f"event_type identity/payload mismatch: {record_id}"
                )

        elif record_type == "therapy_interruption":
            for field in ("drug", "start_date", "end_date"):
                if payload.get(field) != identity.get(field):
                    raise DatasetImportError(
                        f"{field} identity/payload mismatch: {record_id}"
                    )

            if DATASET_PROVENANCE_KEY in (
                payload.get("evidence", {}) or {}
            ):
                raise DatasetImportError(
                    f"reserved evidence key used by payload: {record_id}"
                )

        for assertion in record["provenance"]:
            assertion_id = assertion["assertion_id"]

            if assertion_id in assertion_ids:
                raise DatasetImportError(
                    f"duplicate assertion_id: {assertion_id}"
                )

            assertion_ids.add(assertion_id)

            if assertion["source_id"] not in source_ids:
                raise DatasetImportError(
                    f"unknown source_id in provenance: {record_id}"
                )


def load_dataset_bundle(
    dataset_path: str | Path,
    manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    dataset_path = Path(dataset_path).expanduser().resolve()

    manifest_path = (
        Path(manifest_path).expanduser().resolve()
        if manifest_path
        else dataset_path.with_name("manifest.json")
    )

    dataset = _read_json(dataset_path, "dataset")
    manifest = _read_json(manifest_path, "manifest")

    _validate_schema(dataset)
    _walk_identifier_keys(dataset)
    _validate_record_semantics(dataset)

    canonical_sha = canonical_dataset_sha256(dataset)
    physical_sha = file_sha256(dataset_path)

    expected_canonical_sha = manifest.get(
        "canonical_dataset_sha256"
    )
    expected_file_sha = manifest.get("dataset_file_sha256")

    if expected_canonical_sha != canonical_sha:
        raise DatasetImportError(
            "canonical dataset SHA-256 does not match manifest"
        )

    if expected_file_sha != physical_sha:
        raise DatasetImportError(
            "dataset file SHA-256 does not match manifest"
        )

    for field in (
        "dataset_id",
        "dataset_version",
        "case_ref",
    ):
        if manifest.get(field) != dataset.get(field):
            raise DatasetImportError(
                f"dataset/manifest {field} mismatch"
            )

    record_counts = manifest.get("record_counts") or {}

    if (
        "total" in record_counts
        and int(record_counts["total"]) != len(dataset["records"])
    ):
        raise DatasetImportError(
            "manifest record total does not match dataset"
        )

    actual_by_type: dict[str, int] = {}

    for record in dataset["records"]:
        kind = record["record_type"]
        actual_by_type[kind] = actual_by_type.get(kind, 0) + 1

    for kind, actual_count in actual_by_type.items():
        if (
            kind in record_counts
            and int(record_counts[kind]) != actual_count
        ):
            raise DatasetImportError(
                f"manifest count mismatch for {kind}"
            )

    return {
        "dataset_path": str(dataset_path),
        "manifest_path": str(manifest_path),
        "dataset": dataset,
        "manifest": manifest,
        "canonical_dataset_sha256": canonical_sha,
        "dataset_file_sha256": physical_sha,
    }


def _choice_values(model, field_name: str) -> set[str]:
    field = model._meta.get_field(field_name)

    return {
        str(value)
        for value, _label in field.choices
    }


def _require_choice(
    model,
    field_name: str,
    value: str,
    record_id: str,
) -> None:
    allowed = _choice_values(model, field_name)

    if value not in allowed:
        raise DatasetImportError(
            f"invalid {field_name} for {record_id}"
        )


def _numbers_equal(left: Any, right: Any) -> bool:
    if left is None or right is None:
        return left is None and right is None

    try:
        return float(left) == float(right)
    except (TypeError, ValueError):
        return left == right


def _dataset_provenance(
    bundle: dict[str, Any],
    record: dict[str, Any],
) -> dict[str, Any]:
    dataset = bundle["dataset"]

    return {
        "dataset_id": dataset["dataset_id"],
        "dataset_version": dataset["dataset_version"],
        "canonical_dataset_sha256": bundle[
            "canonical_dataset_sha256"
        ],
        "record_id": record["record_id"],
        "assertions": record["provenance"],
    }


def _strip_interruption_metadata(
    evidence: dict[str, Any] | None,
) -> dict[str, Any]:
    cleaned = dict(evidence or {})

    cleaned.pop(DATASET_PROVENANCE_KEY, None)

    for key in LEGACY_EVIDENCE_PROVENANCE_KEYS:
        cleaned.pop(key, None)

    return cleaned


def _resolve_patient_therapy(
    patient,
    drug: str,
    start_date: str,
    end_date: str | None,
):
    start = date.fromisoformat(start_date)
    end = date.fromisoformat(end_date) if end_date else None

    therapies = patient.therapies.select_related("regimen").order_by(
        "start_date"
    )

    for therapy in therapies:
        therapy_end = therapy.end_date or end or start

        if therapy.start_date > (end or start):
            continue

        if therapy_end < start:
            continue

        components = (therapy.regimen.components or "").lower()

        dose_keys = {
            str(key).lower()
            for key in (therapy.doses or {}).keys()
        }

        if (
            drug.lower() in components
            or drug.lower() in dose_keys
        ):
            return therapy

    return None


def _plan_lab(patient, bundle, record):
    identity = record["identity"]
    payload = record["payload"]

    analyte = str(identity["analyte"])

    _require_choice(
        LongitudinalLabResult,
        "analyte",
        analyte,
        record["record_id"],
    )

    source_quality = str(
        payload.get(
            "source_quality",
            LongitudinalLabResult.SOURCE_QUALITY_EXTRACTED_DOCUMENT,
        )
    )

    _require_choice(
        LongitudinalLabResult,
        "source_quality",
        source_quality,
        record["record_id"],
    )

    result_date = date.fromisoformat(identity["date"])

    assessment = Assessment.objects.filter(
        patient=patient,
        date=result_date,
    ).first()

    provenance = _dataset_provenance(bundle, record)

    desired = {
        "assessment": assessment,
        "value": payload.get("value"),
        "unit": str(payload.get("unit", "")),
        "source_quality": source_quality,
        "provenance": provenance,
        "notes": (
            "Imported from structured research dataset "
            f"{bundle['dataset']['dataset_id']} "
            f"v{bundle['dataset']['dataset_version']}."
        ),
    }

    existing = LongitudinalLabResult.objects.filter(
        patient=patient,
        date=result_date,
        analyte=analyte,
    ).first()

    if existing is None:
        return {
            "action": "created",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "model": LongitudinalLabResult,
            "create": {
                "patient": patient,
                "date": result_date,
                "analyte": analyte,
                **desired,
            },
        }

    conflicts = []

    if not _numbers_equal(existing.value, desired["value"]):
        conflicts.append("value")

    if existing.unit != desired["unit"]:
        conflicts.append("unit")

    if conflicts:
        return {
            "action": "conflict",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "conflict_fields": conflicts,
        }

    updates = {}

    if existing.assessment_id != (
        assessment.id if assessment else None
    ):
        updates["assessment"] = assessment

    for field in (
        "source_quality",
        "provenance",
        "notes",
    ):
        if getattr(existing, field) != desired[field]:
            updates[field] = desired[field]

    if updates:
        return {
            "action": "changed",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "object": existing,
            "updates": updates,
        }

    return {
        "action": "unchanged",
        "record_id": record["record_id"],
        "record_type": record["record_type"],
    }


def _plan_adverse_event(patient, bundle, record):
    identity = record["identity"]
    payload = record["payload"]

    event_type = str(identity["event_type"])

    _require_choice(
        AdverseEvent,
        "event_type",
        event_type,
        record["record_id"],
    )

    event_date = date.fromisoformat(identity["date"])

    provenance = _dataset_provenance(bundle, record)

    desired = {
        "grade": str(payload.get("grade", "")),
        "suspected_drug": str(
            payload.get("suspected_drug", "")
        ),
        "observed_values": payload.get(
            "observed_values",
            {},
        )
        or {},
        "action_taken": str(
            payload.get("action_taken", "")
        ),
        "outcome": str(payload.get("outcome", "")),
        "provenance": provenance,
    }

    existing = AdverseEvent.objects.filter(
        patient=patient,
        date=event_date,
        event_type=event_type,
    ).first()

    if existing is None:
        return {
            "action": "created",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "model": AdverseEvent,
            "create": {
                "patient": patient,
                "date": event_date,
                "event_type": event_type,
                **desired,
            },
        }

    semantic_fields = (
        "grade",
        "suspected_drug",
        "observed_values",
        "action_taken",
        "outcome",
    )

    conflicts = [
        field
        for field in semantic_fields
        if getattr(existing, field) != desired[field]
    ]

    if conflicts:
        return {
            "action": "conflict",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "conflict_fields": conflicts,
        }

    if existing.provenance != provenance:
        return {
            "action": "changed",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "object": existing,
            "updates": {
                "provenance": provenance,
            },
        }

    return {
        "action": "unchanged",
        "record_id": record["record_id"],
        "record_type": record["record_type"],
    }


def _plan_interruption(patient, bundle, record):
    identity = record["identity"]
    payload = record["payload"]

    reason = str(payload["reason"])

    _require_choice(
        TherapyInterruption,
        "reason",
        reason,
        record["record_id"],
    )

    start_date = date.fromisoformat(identity["start_date"])

    end_date = (
        date.fromisoformat(identity["end_date"])
        if identity.get("end_date")
        else None
    )

    drug = str(identity["drug"])

    patient_therapy = _resolve_patient_therapy(
        patient,
        drug,
        identity["start_date"],
        identity.get("end_date"),
    )

    clinical_evidence = dict(
        payload.get("evidence", {}) or {}
    )

    desired_evidence = {
        **clinical_evidence,
        DATASET_PROVENANCE_KEY: _dataset_provenance(
            bundle,
            record,
        ),
    }

    source_quality = str(
        payload.get("source_quality", "clinical_record")
    )

    desired = {
        "patient_therapy": patient_therapy,
        "reason": reason,
        "evidence": desired_evidence,
        "action_taken": str(
            payload.get("action_taken", "")
        ),
        "source_quality": source_quality,
    }

    existing = TherapyInterruption.objects.filter(
        patient=patient,
        drug=drug,
        start_date=start_date,
        end_date=end_date,
    ).first()

    if existing is None:
        return {
            "action": "created",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "model": TherapyInterruption,
            "create": {
                "patient": patient,
                "drug": drug,
                "start_date": start_date,
                "end_date": end_date,
                **desired,
            },
        }

    current_clinical_evidence = _strip_interruption_metadata(
        existing.evidence
    )

    conflicts = []

    if existing.reason != reason:
        conflicts.append("reason")

    if current_clinical_evidence != clinical_evidence:
        conflicts.append("evidence")

    if existing.action_taken != desired["action_taken"]:
        conflicts.append("action_taken")

    if conflicts:
        return {
            "action": "conflict",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "conflict_fields": conflicts,
        }

    updates = {}

    if existing.patient_therapy_id != (
        patient_therapy.id if patient_therapy else None
    ):
        updates["patient_therapy"] = patient_therapy

    if existing.evidence != desired_evidence:
        updates["evidence"] = desired_evidence

    if existing.source_quality != source_quality:
        updates["source_quality"] = source_quality

    if updates:
        return {
            "action": "changed",
            "record_id": record["record_id"],
            "record_type": record["record_type"],
            "object": existing,
            "updates": updates,
        }

    return {
        "action": "unchanged",
        "record_id": record["record_id"],
        "record_type": record["record_type"],
    }


def _plan_record(patient, bundle, record):
    record_type = record["record_type"]

    if record_type == "lab_result":
        return _plan_lab(patient, bundle, record)

    if record_type == "adverse_event":
        return _plan_adverse_event(patient, bundle, record)

    if record_type == "therapy_interruption":
        return _plan_interruption(
            patient,
            bundle,
            record,
        )

    raise DatasetImportError(
        f"unsupported record_type: {record_type}"
    )


def import_dataset(
    patient,
    bundle: dict[str, Any],
    *,
    dry_run: bool,
) -> dict[str, Any]:
    dataset = bundle["dataset"]

    result: dict[str, Any] = {
        "patient_id": patient.id,
        "dry_run": dry_run,
        "dataset_id": dataset["dataset_id"],
        "dataset_version": dataset["dataset_version"],
        "case_ref": dataset["case_ref"],
        "canonical_dataset_sha256": bundle[
            "canonical_dataset_sha256"
        ],
        "records_total": len(dataset["records"]),
        "created": 0,
        "changed": 0,
        "unchanged": 0,
        "conflicts": 0,
        "by_type": {},
        "conflict_records": [],
        "applied": False,
    }

    plans = []

    for record in dataset["records"]:
        plan = _plan_record(patient, bundle, record)
        plans.append(plan)

        action = plan["action"]
        record_type = plan["record_type"]
        counter_key = (
            "conflicts"
            if action == "conflict"
            else action
        )

        result[counter_key] += 1

        by_type = result["by_type"].setdefault(
            record_type,
            {
                "created": 0,
                "changed": 0,
                "unchanged": 0,
                "conflicts": 0,
            },
        )

        by_type[counter_key] += 1

        if action == "conflict":
            result["conflict_records"].append(
                {
                    "record_id": plan["record_id"],
                    "fields": sorted(
                        plan["conflict_fields"]
                    ),
                }
            )

    if dry_run or result["conflicts"]:
        return result

    with transaction.atomic():
        for plan in plans:
            if plan["action"] == "created":
                plan["model"].objects.create(
                    **plan["create"]
                )

            elif plan["action"] == "changed":
                obj = plan["object"]
                updates = plan["updates"]

                for field, value in updates.items():
                    setattr(obj, field, value)

                obj.save(
                    update_fields=sorted(updates.keys())
                )

    result["applied"] = True

    return result

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from django.core.exceptions import ValidationError
from jsonschema import Draft202012Validator, FormatChecker

LEGACY_CONTRACT_VERSION = "legacy-unversioned-v0"
STATE_VECTOR_SCHEMA_VERSION = "twin-state-vector-v1"
PARAMETERS_SCHEMA_VERSION = "twin-parameters-v1"
PARAMETER_UNCERTAINTY_SCHEMA_VERSION = "parameter-uncertainty-v1"
LINEAGE_SCHEMA_VERSION = "twin-lineage-v1"
INTERVENTION_SCHEMA_VERSION = "intervention-v1"
SIMULATION_SUMMARY_SCHEMA_VERSION = "simulation-summary-v1"
COMPARISON_METRICS_SCHEMA_VERSION = "comparison-metrics-v1"
OBSERVATION_RESIDUAL_SCHEMA_VERSION = "observation-residual-v1"
REPORT_MANIFEST_SCHEMA_VERSION = "report-manifest-v1"

SCHEMA_FILES = {
    STATE_VECTOR_SCHEMA_VERSION: "state_vector_v1.schema.json",
    PARAMETERS_SCHEMA_VERSION: "parameters_v1.schema.json",
    PARAMETER_UNCERTAINTY_SCHEMA_VERSION: "parameter_uncertainty_v1.schema.json",
    LINEAGE_SCHEMA_VERSION: "lineage_v1.schema.json",
    INTERVENTION_SCHEMA_VERSION: "intervention_v1.schema.json",
    SIMULATION_SUMMARY_SCHEMA_VERSION: "simulation_summary_v1.schema.json",
    COMPARISON_METRICS_SCHEMA_VERSION: "comparison_metrics_v1.schema.json",
    OBSERVATION_RESIDUAL_SCHEMA_VERSION: "observation_residual_v1.schema.json",
    REPORT_MANIFEST_SCHEMA_VERSION: "report_manifest_v1.schema.json",
}


@lru_cache(maxsize=None)
def load_schema(version: str) -> dict[str, Any]:
    filename = SCHEMA_FILES.get(version)
    if filename is None:
        raise ValidationError(f"Unknown scientific schema version: {version}")
    path = Path(__file__).resolve().parent / "schemas" / filename
    schema = json.loads(path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    return schema


def validate_contract(payload: Any, version: str, *, field_name: str) -> None:
    if version == LEGACY_CONTRACT_VERSION:
        return
    schema = load_schema(version)
    errors = sorted(
        Draft202012Validator(schema, format_checker=FormatChecker()).iter_errors(payload),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if not errors:
        return
    details = []
    for error in errors[:5]:
        location = ".".join(str(part) for part in error.absolute_path) or "<root>"
        details.append(f"{field_name}.{location}: {error.message}")
    raise ValidationError(details)


def validate_twin_contracts(
    *,
    state_vector: Any,
    state_vector_schema_version: str,
    parameters: Any,
    parameters_schema_version: str,
    parameter_uncertainty: Any,
    parameter_uncertainty_schema_version: str,
    lineage: Any,
    lineage_schema_version: str,
) -> None:
    validate_contract(state_vector, state_vector_schema_version, field_name="state_vector")
    validate_contract(parameters, parameters_schema_version, field_name="parameters")
    validate_contract(
        parameter_uncertainty,
        parameter_uncertainty_schema_version,
        field_name="parameter_uncertainty",
    )
    validate_contract(lineage, lineage_schema_version, field_name="lineage")


def validate_residual_contract(instance: Any) -> None:
    validate_contract(
        {
            "predicted_values": instance.predicted_values,
            "observed_values": instance.observed_values,
            "residuals": instance.residuals,
            "normalized_residuals": instance.normalized_residuals,
            "biomarker_weights": instance.biomarker_weights,
        },
        instance.payload_schema_version,
        field_name="observation_residual",
    )


def validate_counterfactual_contracts(instance: Any) -> None:
    validate_contract(
        instance.intervention_definition,
        instance.intervention_schema_version,
        field_name="intervention_definition",
    )
    validate_contract(
        instance.simulation_summary,
        instance.simulation_summary_schema_version,
        field_name="simulation_summary",
    )
    validate_contract(
        instance.comparison_metrics,
        instance.comparison_metrics_schema_version,
        field_name="comparison_metrics",
    )

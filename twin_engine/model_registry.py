from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from django.core.exceptions import ValidationError
from django.utils.module_loading import import_string

from .contracts import SCHEMA_FILES

EVIDENCE_STATUSES = {
    "UNKNOWN",
    "NEEDS_EVIDENCE",
    "OBSERVED",
    "DERIVED",
    "USER_PROVIDED",
    "SIMULATED",
    "HEURISTIC",
    "LITERATURE_INFORMED",
    "HYPOTHETICAL",
    "VALIDATED_EXTERNAL",
}


@dataclass(frozen=True)
class ParameterDefinition:
    name: str
    unit: str
    evidence_status: str


@dataclass(frozen=True)
class ModelRegistration:
    model_id: str
    model_version: str
    input_schema_version: str
    output_schema_version: str
    endpoint_definition_version: str
    parameter_definitions: tuple[ParameterDefinition, ...]
    entry_point: str
    model_card: str
    compatible_major_versions: tuple[int, ...]
    invalidated_prior_versions: tuple[str, ...]

    @property
    def major_version(self) -> int:
        match = re.search(r"(?:^|[-v])v?(\d+)(?:\.|$)", self.model_version)
        if match is None:
            raise ValidationError(f"Model version has no parseable major version: {self.model_version}")
        return int(match.group(1))

    @property
    def units(self) -> dict[str, str]:
        return {item.name: item.unit for item in self.parameter_definitions}


def registry_path() -> Path:
    return Path(__file__).resolve().parent / "model_registry.json"


def registry_sha256() -> str:
    return hashlib.sha256(registry_path().read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def load_model_registry() -> dict[str, ModelRegistration]:
    payload: dict[str, Any] = json.loads(registry_path().read_text(encoding="utf-8"))
    if payload.get("registry_version") != "model-registry-v1":
        raise ValidationError("Unsupported model registry version")

    registrations: dict[str, ModelRegistration] = {}
    for item in payload.get("models", []):
        model_id = str(item.get("model_id") or "")
        if not model_id or model_id in registrations:
            raise ValidationError(f"Duplicate or empty model_id: {model_id!r}")
        for schema_field in ("input_schema_version", "output_schema_version"):
            schema_version = str(item.get(schema_field) or "")
            if schema_version not in SCHEMA_FILES:
                raise ValidationError(f"{model_id} references unknown {schema_field}: {schema_version}")
        parameters = tuple(
            ParameterDefinition(
                name=str(parameter["name"]),
                unit=str(parameter["unit"]),
                evidence_status=str(parameter["evidence_status"]),
            )
            for parameter in item.get("parameter_definitions", [])
        )
        parameter_names = [parameter.name for parameter in parameters]
        if len(parameter_names) != len(set(parameter_names)):
            raise ValidationError(f"{model_id} contains duplicate parameter definitions")
        for parameter in parameters:
            if not parameter.name or not parameter.unit:
                raise ValidationError(f"{model_id} has an incomplete parameter definition")
            if parameter.evidence_status not in EVIDENCE_STATUSES:
                raise ValidationError(
                    f"{model_id}.{parameter.name} has unknown evidence status: {parameter.evidence_status}"
                )
        registration = ModelRegistration(
            model_id=model_id,
            model_version=str(item["model_version"]),
            input_schema_version=str(item["input_schema_version"]),
            output_schema_version=str(item["output_schema_version"]),
            endpoint_definition_version=str(item["endpoint_definition_version"]),
            parameter_definitions=parameters,
            entry_point=str(item["entry_point"]),
            model_card=str(item["model_card"]),
            compatible_major_versions=tuple(int(value) for value in item["compatible_major_versions"]),
            invalidated_prior_versions=tuple(str(value) for value in item["invalidated_prior_versions"]),
        )
        if registration.major_version not in registration.compatible_major_versions:
            raise ValidationError(f"{model_id} does not declare compatibility with its own major version")
        model_card_path = registry_path().parent.parent / registration.model_card
        if not model_card_path.is_file():
            raise ValidationError(f"{model_id} model card is unavailable: {registration.model_card}")
        registrations[model_id] = registration
    if not registrations:
        raise ValidationError("Model registry contains no registrations")
    return registrations


def get_model_registration(model_id: str) -> ModelRegistration:
    try:
        return load_model_registry()[model_id]
    except KeyError as exc:
        raise ValidationError(f"Unregistered scientific model: {model_id}") from exc


def validate_registry_entry_points() -> None:
    for registration in load_model_registry().values():
        try:
            import_string(registration.entry_point)
        except ImportError as exc:
            raise ValidationError(
                f"{registration.model_id} entry point is unavailable: {registration.entry_point}"
            ) from exc

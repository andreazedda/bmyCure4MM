from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from django.core.exceptions import ValidationError

from .model_registry import get_model_registration
from .models import ResearchRunInvalidation, ResearchRunManifest
from .run_manifest import canonical_sha256

DIRECTLY_COMPARABLE = "DIRECTLY_COMPARABLE"
NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN = "NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN"


@dataclass(frozen=True)
class ComparabilityResult:
    status: str
    reasons: tuple[str, ...]
    compatibility: dict[str, bool]

    @property
    def directly_comparable(self) -> bool:
        return self.status == DIRECTLY_COMPARABLE


def compare_manifests(
    left: ResearchRunManifest,
    right: ResearchRunManifest,
) -> ComparabilityResult:
    reasons: list[str] = []
    compatibility: dict[str, bool] = {}

    compatibility["not_invalidated"] = not left.invalidations.exists() and not right.invalidations.exists()
    if not compatibility["not_invalidated"]:
        reasons.append("one or both run manifests have an immutable invalidation record")

    compatibility["model_identity"] = left.model_id == right.model_id
    if not compatibility["model_identity"]:
        reasons.append("model_id differs")

    compatibility["model_major_version"] = False
    compatibility["schema_versions"] = False
    compatibility["endpoint_definition"] = False
    compatibility["parameter_units"] = False
    if compatibility["model_identity"]:
        registration = get_model_registration(left.model_id)
        left_major = _major(left.model_version)
        right_major = _major(right.model_version)
        declared = set(registration.compatible_major_versions)
        compatibility["model_major_version"] = left_major == right_major and left_major in declared
        if not compatibility["model_major_version"]:
            reasons.append("model major versions are not declared compatible")
        same_registry = left.model_registry_sha256 == right.model_registry_sha256
        compatibility["schema_versions"] = same_registry
        compatibility["endpoint_definition"] = same_registry
        compatibility["parameter_units"] = same_registry
        if not same_registry:
            reasons.append(
                "model registry identity differs; schema, endpoint, and unit compatibility require rerun"
            )

    compatibility["dataset_identity"] = (
        left.dataset_id,
        left.dataset_version,
        left.dataset_sha256,
        left.record_subset_sha256,
    ) == (
        right.dataset_id,
        right.dataset_version,
        right.dataset_sha256,
        right.record_subset_sha256,
    )
    if not compatibility["dataset_identity"]:
        reasons.append("dataset or record-subset identity differs")
    compatibility["dataset_bound"] = (
        all(
            value not in {"", "UNBOUND", "UNAVAILABLE"}
            for value in (left.dataset_id, left.dataset_version, left.dataset_sha256)
        )
        and len(left.dataset_sha256) == 64
    )
    if not compatibility["dataset_bound"]:
        reasons.append("dataset identity is explicitly unbound or unavailable")

    compatibility["intended_use"] = left.intended_use_level == right.intended_use_level
    if not compatibility["intended_use"]:
        reasons.append("intended-use levels differ")

    status = DIRECTLY_COMPARABLE if all(compatibility.values()) else NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN
    return ComparabilityResult(status=status, reasons=tuple(reasons), compatibility=compatibility)


def invalidate_manifest(
    manifest: ResearchRunManifest,
    *,
    change_kind: str,
    replacement_identity: str,
    reason: str,
) -> ResearchRunInvalidation:
    field_by_change = {
        ResearchRunInvalidation.CHANGE_CLINICAL_VALUE: "record_subset_sha256",
        ResearchRunInvalidation.CHANGE_DATASET: "dataset_sha256",
        ResearchRunInvalidation.CHANGE_UNIT_MAPPING: "model_registry_sha256",
        ResearchRunInvalidation.CHANGE_MODEL: "model_version",
        ResearchRunInvalidation.CHANGE_PARAMETER_DEFAULT: "configuration_sha256",
        ResearchRunInvalidation.CHANGE_DEPENDENCY: "dependency_lock_sha256",
        ResearchRunInvalidation.CHANGE_CONFIGURATION: "configuration_sha256",
    }
    try:
        identity_field = field_by_change[change_kind]
    except KeyError as exc:
        raise ValidationError(f"Unsupported invalidation change kind: {change_kind}") from exc
    previous_identity = str(getattr(manifest, identity_field))
    if not replacement_identity.strip() or replacement_identity == previous_identity:
        raise ValidationError("Invalidation requires a distinct replacement identity")
    if not reason.strip():
        raise ValidationError("Invalidation requires a reason")
    payload: dict[str, Any] = {
        "run_id": manifest.run_id,
        "change_kind": change_kind,
        "identity_field": identity_field,
        "previous_identity": previous_identity,
        "replacement_identity": replacement_identity,
        "reason": reason,
    }
    return ResearchRunInvalidation.objects.create(
        manifest=manifest,
        change_kind=change_kind,
        previous_identity=previous_identity,
        replacement_identity=replacement_identity,
        change_sha256=canonical_sha256(payload),
        reason=reason,
    )


def _major(version: str) -> int:
    parts = version.rsplit("-v", 1)
    if len(parts) != 2:
        raise ValidationError(f"Model version has no parseable major version: {version}")
    try:
        return int(parts[1].split(".", 1)[0])
    except ValueError as exc:
        raise ValidationError(f"Model version has no parseable major version: {version}") from exc

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

from django.conf import settings
from django.core.exceptions import ValidationError
from django.utils import timezone

from mmportal.governance import (
    API_VERSION,
    APP_VERSION,
    EVIDENCE_GRAPH_VERSION,
    INTENDED_USE_LEVEL,
    MODEL_CARD_VERSION,
    REPORT_TEMPLATE_VERSION,
    VALIDATION_PROTOCOL_VERSION,
    EpistemicLabel,
)

from .contracts import REPORT_MANIFEST_SCHEMA_VERSION, validate_contract
from .model_registry import get_model_registration, registry_sha256
from .validators import validate_no_direct_identifier_in_artifact_payload

RUN_MANIFEST_CONTRACT_VERSION = "research-run-manifest-v1"
SIMULATION_METADATA_CONTRACT_VERSION = "simulation-metadata-v1"
UNBOUND_DATASET = "UNBOUND"
UNAVAILABLE = "UNAVAILABLE"
UNAVAILABLE_LOCAL_CONTAINER = "UNAVAILABLE_LOCAL"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_CONTAINER_DIGEST_RE = re.compile(r"^sha256:[0-9a-f]{64}$")


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)


def canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repository_root() -> Path:
    return Path(__file__).resolve().parent.parent


def dependency_lock_sha256() -> str:
    lock_path = repository_root() / "uv.lock"
    if not lock_path.is_file():
        raise ValidationError("Canonical dependency lock uv.lock is unavailable")
    return file_sha256(lock_path)


def current_git_sha() -> str:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repository_root(),
            capture_output=True,
            check=True,
            text=True,
        )
    except (OSError, subprocess.SubprocessError):
        return UNAVAILABLE
    value = completed.stdout.strip()
    return value if re.fullmatch(r"[0-9a-f]{40}", value) else UNAVAILABLE


def database_schema_version() -> str:
    from django.db.migrations.recorder import MigrationRecorder

    applied = sorted(
        f"{app}.{name}" for app, name in MigrationRecorder.Migration.objects.values_list("app", "name")
    )
    return f"db-schema-sha256:{canonical_sha256(applied)}"


def container_digest() -> str:
    value = os.environ.get("CONTAINER_IMAGE_DIGEST", "").strip()
    if not value:
        return UNAVAILABLE_LOCAL_CONTAINER
    if not _CONTAINER_DIGEST_RE.fullmatch(value):
        raise ValidationError("CONTAINER_IMAGE_DIGEST must be sha256:<64 lowercase hex>")
    return value


def _dataset_identity(twin_state: Any | None, input_payload: Any) -> tuple[str, str, str, str]:
    lineage = dict(getattr(twin_state, "lineage", None) or {})
    binding = dict(lineage.get("dataset_binding") or {})
    dataset_id = str(binding.get("dataset_id") or UNBOUND_DATASET)
    dataset_version = str(binding.get("dataset_version") or UNBOUND_DATASET)
    dataset_sha = str(binding.get("canonical_dataset_sha256") or UNBOUND_DATASET)
    computational_input = dict(lineage.get("computational_input") or {})
    subset_sha = str(computational_input.get("sha256") or canonical_sha256(input_payload))
    if dataset_sha != UNBOUND_DATASET and not _SHA256_RE.fullmatch(dataset_sha):
        raise ValidationError("Lineage dataset_sha256 is not a canonical SHA-256")
    if not _SHA256_RE.fullmatch(subset_sha):
        raise ValidationError("record_subset_sha256 is not a canonical SHA-256")
    return dataset_id, dataset_version, dataset_sha, subset_sha


@dataclass(frozen=True)
class PreparedRunIdentity:
    run_id: str
    status: str
    version_vector: dict[str, Any]


def prepare_run_identity(
    *,
    model_id: str,
    input_payload: Any,
    configuration_payload: Any,
    twin_state: Any | None = None,
    random_seed: int | None = None,
    epistemic_label: EpistemicLabel | str = EpistemicLabel.SIMULATED,
    status: str = "completed",
    created_at: datetime | None = None,
) -> PreparedRunIdentity:
    registration = get_model_registration(model_id)
    dataset_id, dataset_version, dataset_sha, subset_sha = _dataset_identity(twin_state, input_payload)
    timestamp = created_at or timezone.now()
    if timezone.is_naive(timestamp):
        raise ValidationError("Run manifest created_at must be timezone-aware")
    label = epistemic_label.value if isinstance(epistemic_label, EpistemicLabel) else str(epistemic_label)
    if label not in {item.value for item in EpistemicLabel}:
        raise ValidationError(f"Unknown epistemic label: {label}")
    if status not in {"completed", "failed"}:
        raise ValidationError(f"Unsupported run status: {status}")
    version_vector = {
        "app_version": APP_VERSION,
        "api_version": API_VERSION,
        "db_schema_version": database_schema_version(),
        "dataset_id": dataset_id,
        "dataset_version": dataset_version,
        "dataset_sha256": dataset_sha,
        "record_subset_sha256": subset_sha,
        "model_id": model_id,
        "model_version": registration.model_version,
        "model_card_version": MODEL_CARD_VERSION,
        "configuration_sha256": canonical_sha256(configuration_payload),
        "evidence_graph_version": EVIDENCE_GRAPH_VERSION,
        "validation_protocol_version": VALIDATION_PROTOCOL_VERSION,
        "report_template_version": REPORT_TEMPLATE_VERSION,
        "git_sha": current_git_sha(),
        "container_digest": container_digest(),
        "dependency_lock_sha256": dependency_lock_sha256(),
        "model_registry_sha256": registry_sha256(),
        "random_seed": int(random_seed or 0),
        "created_at": timestamp.isoformat(),
        "intended_use_level": INTENDED_USE_LEVEL,
        "epistemic_label": label,
    }
    identity_hash = canonical_sha256({"status": status, "version_vector": version_vector})
    compact_time = timestamp.strftime("%Y%m%dT%H%M%S%fZ")
    return PreparedRunIdentity(
        run_id=f"RUN_{compact_time}_{identity_hash[:12]}",
        status=status,
        version_vector=version_vector,
    )


def hash_only_artifact(name: str, role: str, payload: Any) -> dict[str, Any]:
    digest = canonical_sha256(payload)
    return {
        "name": name,
        "role": role,
        "media_type": "application/json",
        "uri": f"hash://sha256/{digest}",
        "sha256": digest,
        "hash_only": True,
    }


def file_artifact(name: str, role: str, path: str | Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    media_root = Path(settings.MEDIA_ROOT).resolve()
    try:
        relative = resolved.relative_to(media_root)
    except ValueError as exc:
        raise ValidationError("Research artifacts must be stored beneath MEDIA_ROOT") from exc
    if not resolved.is_file():
        raise ValidationError(f"Research artifact is unavailable: {resolved}")
    return {
        "name": name,
        "role": role,
        "media_type": "application/json",
        "uri": f"media://{relative.as_posix()}",
        "sha256": file_sha256(resolved),
        "hash_only": False,
    }


def _manifest_unsigned_payload(
    identity: PreparedRunIdentity,
    artifact_manifest: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "run_id": identity.run_id,
        "version_vector": identity.version_vector,
        "artifact_manifest": artifact_manifest,
    }


def build_manifest_payload(
    identity: PreparedRunIdentity,
    artifact_manifest: Iterable[dict[str, Any]],
) -> dict[str, Any]:
    artifacts = list(artifact_manifest)
    unsigned = _manifest_unsigned_payload(identity, artifacts)
    payload = {**unsigned, "manifest_sha256": canonical_sha256(unsigned)}
    validate_contract(payload, REPORT_MANIFEST_SCHEMA_VERSION, field_name="run_manifest")
    validate_no_direct_identifier_in_artifact_payload(payload)
    return payload


def manifest_storage_path(run_id: str) -> Path:
    if not re.fullmatch(r"RUN_[A-Za-z0-9_.-]+", run_id):
        raise ValidationError("Invalid research run_id")
    return Path(settings.MEDIA_ROOT) / "research_manifests" / f"{run_id}.json"


def persist_manifest_artifact(payload: dict[str, Any]) -> tuple[str, Path]:
    path = manifest_storage_path(str(payload["run_id"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)
    relative = path.resolve().relative_to(Path(settings.MEDIA_ROOT).resolve())
    return f"media://{relative.as_posix()}", path


def validate_manifest_instance(instance: Any) -> None:
    vector = instance.version_vector()
    identity = PreparedRunIdentity(
        run_id=instance.run_id,
        status=instance.status,
        version_vector=vector,
    )
    expected = build_manifest_payload(identity, instance.artifact_manifest)
    if instance.contract_version != RUN_MANIFEST_CONTRACT_VERSION:
        raise ValidationError("Unsupported research run manifest contract")
    if instance.manifest_sha256 != expected["manifest_sha256"]:
        raise ValidationError("Research run manifest hash does not match its contents")


def verify_manifest_artifacts(manifest: Any) -> None:
    validate_manifest_instance(manifest)
    manifest_path = _resolve_media_uri(manifest.manifest_artifact_uri)
    if not manifest_path.is_file():
        raise ValidationError("Run manifest artifact is unavailable")
    stored = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = build_manifest_payload(
        PreparedRunIdentity(manifest.run_id, manifest.status, manifest.version_vector()),
        manifest.artifact_manifest,
    )
    if stored != expected:
        raise ValidationError("Run manifest artifact content does not match the database record")
    for artifact in manifest.artifact_manifest:
        if artifact.get("hash_only"):
            expected_uri = f"hash://sha256/{artifact['sha256']}"
            if artifact.get("uri") != expected_uri or not _SHA256_RE.fullmatch(str(artifact["sha256"])):
                raise ValidationError(f"Invalid hash-only artifact reference: {artifact.get('name')}")
            continue
        path = _resolve_media_uri(str(artifact.get("uri") or ""))
        if not path.is_file() or file_sha256(path) != artifact.get("sha256"):
            raise ValidationError(f"Artifact integrity check failed: {artifact.get('name')}")


def _resolve_media_uri(uri: str) -> Path:
    if not uri.startswith("media://"):
        raise ValidationError("Artifact URI must use the media:// scheme")
    media_root = Path(settings.MEDIA_ROOT).resolve()
    candidate = (media_root / uri.removeprefix("media://")).resolve()
    try:
        candidate.relative_to(media_root)
    except ValueError as exc:
        raise ValidationError("Artifact URI escapes MEDIA_ROOT") from exc
    return candidate

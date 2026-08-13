from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any

from django.core.exceptions import ValidationError
from django.db import transaction

from mmportal.governance import CURRENT_RESEARCH_MODEL_VERSION

from .model_registry import get_model_registration
from .models import ResearchRunManifest, SimulationRunMetadata
from .run_manifest import (
    RUN_MANIFEST_CONTRACT_VERSION,
    SIMULATION_METADATA_CONTRACT_VERSION,
    PreparedRunIdentity,
    build_manifest_payload,
    file_artifact,
    hash_only_artifact,
    persist_manifest_artifact,
    prepare_run_identity,
)

CURRENT_MODEL_VERSION = CURRENT_RESEARCH_MODEL_VERSION


def hash_json(payload: Any) -> str:
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def hash_file(path: str | Path) -> str:
    file_path = Path(path)
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_code_commit_hash() -> str:
    repo_root = Path(__file__).resolve().parent.parent
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            capture_output=True,
            check=True,
            text=True,
        )
    except Exception:
        return ""
    return completed.stdout.strip()


def collect_drug_preset_hashes() -> dict[str, str]:
    presets_dir = Path(__file__).resolve().parent.parent / "simulator" / "presets" / "drugs"
    if not presets_dir.exists():
        return {}
    return {path.name: hash_file(path) for path in sorted(presets_dir.glob("*.yaml"))}


def collect_twin_config_hash() -> str:
    from .observation_model import DEFAULT_OBSERVATION_PARAMETERS

    twin_config_path = Path(__file__).resolve().parent.parent / "simulator" / "presets" / "twin_risk.yaml"
    payload = {
        "model_version": CURRENT_MODEL_VERSION,
        "twin_risk_hash": hash_file(twin_config_path) if twin_config_path.exists() else "",
        "observation_defaults": DEFAULT_OBSERVATION_PARAMETERS,
    }
    return hash_json(payload)


def prepare_simulation_identity(
    *,
    model_id: str,
    solver_name: str,
    input_payload: Any,
    solver_parameters: dict[str, Any] | None = None,
    twin_state=None,
    random_seed: int | None = None,
    status: str = ResearchRunManifest.STATUS_COMPLETED,
) -> PreparedRunIdentity:
    return prepare_run_identity(
        model_id=model_id,
        input_payload=input_payload,
        configuration_payload={
            "solver_name": solver_name,
            "solver_parameters": solver_parameters or {},
            "twin_config_sha256": collect_twin_config_hash(),
            "drug_preset_hashes": collect_drug_preset_hashes(),
        },
        twin_state=twin_state,
        random_seed=random_seed,
        status=status,
    )


def record_simulation_metadata(
    *,
    model_id: str = "patient_twin_state_model",
    model_version: str | None = None,
    solver_name: str,
    input_payload: Any,
    solver_parameters: dict[str, Any] | None = None,
    output_payload: Any | None = None,
    simulation_attempt=None,
    counterfactual_run=None,
    twin_state=None,
    random_seed: int | None = None,
    prepared_identity: PreparedRunIdentity | None = None,
    artifact_paths: list[tuple[str, str, str | Path]] | None = None,
    status: str = ResearchRunManifest.STATUS_COMPLETED,
) -> SimulationRunMetadata:
    registration = get_model_registration(model_id)
    if model_version is not None and model_version != registration.model_version:
        raise ValidationError(
            f"model_version {model_version!r} does not match registered {registration.model_version!r}"
        )
    identity = prepared_identity or prepare_simulation_identity(
        model_id=model_id,
        solver_name=solver_name,
        input_payload=input_payload,
        solver_parameters=solver_parameters,
        twin_state=twin_state,
        random_seed=random_seed,
        status=status,
    )
    if identity.version_vector["model_id"] != model_id:
        raise ValidationError("Prepared run identity model_id does not match the metadata request")
    if identity.status != status:
        raise ValidationError("Prepared run identity status does not match the metadata request")

    artifacts = [hash_only_artifact("computational_input", "input", input_payload)]
    if output_payload is not None:
        artifacts.append(hash_only_artifact("computational_output", "output", output_payload))
    artifacts.extend(file_artifact(name, role, path) for name, role, path in (artifact_paths or []))
    manifest_payload = build_manifest_payload(identity, artifacts)
    manifest_uri = ""
    manifest_path: Path | None = None
    try:
        with transaction.atomic():
            manifest_uri, manifest_path = persist_manifest_artifact(manifest_payload)
            vector = identity.version_vector
            manifest = ResearchRunManifest.objects.create(
                run_id=identity.run_id,
                contract_version=RUN_MANIFEST_CONTRACT_VERSION,
                status=identity.status,
                app_version=vector["app_version"],
                api_version=vector["api_version"],
                db_schema_version=vector["db_schema_version"],
                dataset_id=vector["dataset_id"],
                dataset_version=vector["dataset_version"],
                dataset_sha256=vector["dataset_sha256"],
                record_subset_sha256=vector["record_subset_sha256"],
                model_id=vector["model_id"],
                model_version=vector["model_version"],
                model_card_version=vector["model_card_version"],
                configuration_sha256=vector["configuration_sha256"],
                evidence_graph_version=vector["evidence_graph_version"],
                validation_protocol_version=vector["validation_protocol_version"],
                report_template_version=vector["report_template_version"],
                git_sha=vector["git_sha"],
                container_digest=vector["container_digest"],
                dependency_lock_sha256=vector["dependency_lock_sha256"],
                model_registry_sha256=vector["model_registry_sha256"],
                random_seed=vector["random_seed"],
                intended_use_level=vector["intended_use_level"],
                epistemic_label=vector["epistemic_label"],
                artifact_manifest=artifacts,
                manifest_artifact_uri=manifest_uri,
                manifest_sha256=manifest_payload["manifest_sha256"],
                created_at=datetime.fromisoformat(vector["created_at"]),
            )
            return SimulationRunMetadata.objects.create(
                contract_version=SIMULATION_METADATA_CONTRACT_VERSION,
                manifest=manifest,
                simulation_attempt=simulation_attempt,
                counterfactual_run=counterfactual_run,
                twin_state=twin_state,
                model_version=registration.model_version,
                code_commit_hash=vector["git_sha"],
                twin_config_hash=collect_twin_config_hash(),
                drug_preset_hashes=collect_drug_preset_hashes(),
                solver_name=solver_name,
                solver_parameters=solver_parameters or {},
                input_hash=hash_json(input_payload),
                output_hash=hash_json(output_payload) if output_payload is not None else "",
                random_seed=vector["random_seed"],
            )
    except Exception:
        if manifest_path is not None and manifest_path.is_file():
            manifest_path.unlink()
        raise

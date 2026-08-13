from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from mmportal.governance import CURRENT_RESEARCH_MODEL_VERSION

from .models import SimulationRunMetadata


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
    return {
        path.name: hash_file(path)
        for path in sorted(presets_dir.glob("*.yaml"))
    }


def collect_twin_config_hash() -> str:
    from .observation_model import DEFAULT_OBSERVATION_PARAMETERS

    twin_config_path = Path(__file__).resolve().parent.parent / "simulator" / "presets" / "twin_risk.yaml"
    payload = {
        "model_version": CURRENT_MODEL_VERSION,
        "twin_risk_hash": hash_file(twin_config_path) if twin_config_path.exists() else "",
        "observation_defaults": DEFAULT_OBSERVATION_PARAMETERS,
    }
    return hash_json(payload)


def record_simulation_metadata(
    *,
    model_version: str,
    solver_name: str,
    input_payload: Any,
    solver_parameters: dict[str, Any] | None = None,
    output_payload: Any | None = None,
    simulation_attempt=None,
    counterfactual_run=None,
    twin_state=None,
    random_seed: int | None = None,
) -> SimulationRunMetadata:
    return SimulationRunMetadata.objects.create(
        simulation_attempt=simulation_attempt,
        counterfactual_run=counterfactual_run,
        twin_state=twin_state,
        model_version=model_version,
        code_commit_hash=get_code_commit_hash(),
        twin_config_hash=collect_twin_config_hash(),
        drug_preset_hashes=collect_drug_preset_hashes(),
        solver_name=solver_name,
        solver_parameters=solver_parameters or {},
        input_hash=hash_json(input_payload),
        output_hash=hash_json(output_payload) if output_payload is not None else "",
        random_seed=random_seed,
    )

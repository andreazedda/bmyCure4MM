from __future__ import annotations

import contextlib
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import yaml
from django.conf import settings


class ToolRunError(RuntimeError):
    """Raised when an external tool fails to execute successfully."""

    def __init__(self, message: str, stdout: str = "", stderr: str = "") -> None:
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr


BASE_DIR = Path(__file__).resolve().parent.parent
PIPELINES_DIR = BASE_DIR / "pipelines"
MODULES_DIR = BASE_DIR / "modules"
LIGAND_DIR = BASE_DIR / "LigandSimilaritySearcher"


def _media_root() -> Path:
    root = Path(getattr(settings, "MEDIA_ROOT", BASE_DIR / "media"))
    root.mkdir(parents=True, exist_ok=True)
    return root


def _chem_media_root() -> Path:
    path = _media_root() / "chem"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _pdb_cache_dir() -> Path:
    path = _media_root() / "pdb_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _chem_cache_dir() -> Path:
    path = _media_root() / "chem_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _binding_cache_dir() -> Path:
    path = _chem_cache_dir() / "binding"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _drug_cache_dir() -> Path:
    path = _chem_cache_dir() / "drug"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _sim_cache_dir() -> Path:
    path = _chem_cache_dir() / "similarity"
    path.mkdir(parents=True, exist_ok=True)
    return path


@contextlib.contextmanager
def temporary_yaml_override(path: Path, mutator) -> Iterable[None]:
    original = path.read_text(encoding="utf-8")
    data = yaml.safe_load(original) if original.strip() else {}
    mutator(data)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    try:
        yield
    finally:
        path.write_text(original, encoding="utf-8")


def _compose_logs(stdout: str, stderr: str) -> str:
    logs = []
    if stdout:
        logs.append(stdout)
    if stderr:
        logs.append("--- STDERR ---")
        logs.append(stderr)
    return "\n".join(logs)


def _run_command(
    command: list[str],
    env: Optional[dict] = None,
    cwd: Optional[Path] = None,
    timeout: int = 900,
) -> tuple[int, str, str]:
    """Execute a command capturing stdout/stderr."""
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    try:
        result = subprocess.run(
            command,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            env=full_env,
            check=False,
            timeout=timeout,
        )
        stdout = result.stdout.strip()
        stderr = result.stderr.strip()
        return result.returncode, stdout, stderr
    except subprocess.TimeoutExpired as exc:  # pragma: no cover - defensive
        raise ToolRunError("Tool execution timed out", stdout=exc.stdout or "", stderr=exc.stderr or "") from exc


def _snapshot(paths: Iterable[Path]) -> Dict[Path, int]:
    snapshot: Dict[Path, int] = {}
    for path in paths:
        if path.exists():
            snapshot[path] = path.stat().st_mtime_ns
    return snapshot


def _new_files(before: Dict[Path, int], after: Iterable[Path]) -> list[Path]:
    updated: list[Path] = []
    for path in after:
        if path.exists():
            mtime = path.stat().st_mtime_ns
            if path not in before or before[path] != mtime:
                updated.append(path)
    return updated


def _resolve_outdir(outdir: Path | str | None, default_subdir: str) -> Path:
    if outdir:
        resolved = Path(outdir)
    else:
        env_override = os.environ.get("MM_OUTDIR")
    resolved = Path(env_override) if env_override else _chem_media_root() / default_subdir
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def run_settings_generator() -> str:
    script = PIPELINES_DIR / "processes" / "settings_generator.py"
    code, stdout, stderr = _run_command([sys.executable, str(script)], cwd=script.parent)
    if code != 0:
        raise ToolRunError("settings_generator failed", stdout, stderr)
    return _compose_logs(stdout, stderr)


def run_paths_generator() -> str:
    script = PIPELINES_DIR / "processes" / "paths_generator.py"
    code, stdout, stderr = _run_command([sys.executable, str(script)], cwd=script.parent)
    if code != 0:
        raise ToolRunError("paths_generator failed", stdout, stderr)
    return _compose_logs(stdout, stderr)


def run_drug_parameter_evaluator(
    smiles: str | None = None,
    cid: int | str | None = None,
    outdir: Path | str | None = None,
) -> Tuple[Optional[Path], str, str]:
    """Trigger the drug parameter evaluator script and collect the HTML report."""
    script = PIPELINES_DIR / "processes" / "drug_parameter_evaluator.py"
    config_path = PIPELINES_DIR / "configs" / "drug_parameter_evaluator.yaml"
    outdir_path = _resolve_outdir(outdir, "drug_tmp")

    outputs_dir = PIPELINES_DIR / "data" / "outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    before = _snapshot(outputs_dir.glob("*.html"))

    def mutate_config(data):
        if cid:
            data["drug_id"] = str(cid)

    env: dict[str, str] = {"MM_OUTDIR": str(outdir_path)}
    if smiles:
        env["MM_SMILES_OVERRIDE"] = smiles
    if cid:
        env["MM_CID_OVERRIDE"] = str(cid)

    with temporary_yaml_override(config_path, mutate_config):
        code, stdout, stderr = _run_command(
            [sys.executable, str(script)],
            env=env,
            cwd=script.parent,
        )
    if code != 0:
        raise ToolRunError(f"{script.name} failed", stdout, stderr)

    after = list(outputs_dir.glob("*.html"))
    new = _new_files(before, after)
    target_file: Optional[Path] = None
    if new:
        newest = max(new, key=lambda p: p.stat().st_mtime)
        dest_name = f"drug_{cid or 'custom'}_{newest.stem}.html"
        dest_path = outdir_path / dest_name
        shutil.copy2(newest, dest_path)
        target_file = dest_path
        cache_target = _drug_cache_dir() / dest_name
        shutil.copy2(newest, cache_target)
    return target_file, stdout, stderr


def _prepare_pdb_cache(pdb_id: str) -> None:
    cached = _pdb_cache_dir() / f"{pdb_id}.pdb"
    module_pdb = MODULES_DIR / "binding_visualizer" / f"{pdb_id}.pdb"
    if cached.exists() and not module_pdb.exists():
        shutil.copy2(cached, module_pdb)


def _store_pdb_cache(pdb_id: str) -> None:
    module_pdb = MODULES_DIR / "binding_visualizer" / f"{pdb_id}.pdb"
    if module_pdb.exists():
        cache_target = _pdb_cache_dir() / f"{pdb_id}.pdb"
        shutil.copy2(module_pdb, cache_target)


def run_binding_visualizer(
    pdb_id: str,
    ligand: str | None = "",
    outdir: Path | str | None = None,
) -> Tuple[Optional[Path], str, str]:
    """Run the binding visualizer workflow and collect the generated HTML viewer."""
    script = MODULES_DIR / "binding_visualizer" / "binding_visualizer.py"
    config_path = MODULES_DIR / "binding_visualizer" / "binding_visualizer.yaml"
    outdir_path = _resolve_outdir(outdir, "binding_tmp")

    should_generate_pdf = os.environ.get("MM_BINDING_PDF") == "1"

    _prepare_pdb_cache(pdb_id)

    module_outputs = MODULES_DIR / "binding_visualizer"
    before = _snapshot(module_outputs.glob("*_structure_viewer.html"))

    def mutate_config(data):
        data["pdb_id"] = pdb_id
        if ligand:
            data["ligand"] = ligand
        pdf_options = data.get("pdf_options", {})
        pdf_options["cleanup_aux_files"] = False
        data["pdf_options"] = pdf_options
        data["generate_pdf"] = should_generate_pdf

    env = {
        "MM_OUTDIR": str(outdir_path),
        "MM_PDB_ID": pdb_id,
        "MM_LIGAND": ligand or "",
        "BROWSER": "true",
    }

    with temporary_yaml_override(config_path, mutate_config):
        code, stdout, stderr = _run_command(
            [sys.executable, str(script)],
            env=env,
            cwd=script.parent,
        )
    if code != 0:
        raise ToolRunError(f"{script.name} failed", stdout, stderr)

    after = list(module_outputs.glob("*_structure_viewer.html"))
    new = _new_files(before, after)
    target_file: Optional[Path] = None
    if new:
        newest = max(new, key=lambda p: p.stat().st_mtime)
        dest_name = f"binding_{pdb_id}_{newest.stem}.html"
        dest_path = outdir_path / dest_name
        shutil.copy2(newest, dest_path)
        target_file = dest_path
        cache_target = _binding_cache_dir() / f"{pdb_id}.html"
        shutil.copy2(newest, cache_target)

    _store_pdb_cache(pdb_id)
    return target_file, stdout, stderr


def run_similarity_search(
    smiles: str,
    out_csv: Path | str | None = None,
) -> Tuple[Optional[Path], str, str]:
    """Execute ligand similarity search and copy resulting CSV."""
    script = LIGAND_DIR / "sources" / "ligand_similarity_searcher.py"
    config_path = LIGAND_DIR / "configs" / "ligand_similarity_searcher.yaml"
    lig_data_dir = LIGAND_DIR / "data"
    lig_data_dir.mkdir(parents=True, exist_ok=True)
    output_dir = lig_data_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    before = _snapshot(output_dir.glob("*.csv"))

    outdir_path = _resolve_outdir(out_csv, "similarity_tmp")

    def mutate_config(data):
        data["input_smiles"] = smiles

    env = {
        "MM_OUTDIR": str(outdir_path),
        "MM_SMILES_OVERRIDE": smiles,
    }

    with temporary_yaml_override(config_path, mutate_config):
        code, stdout, stderr = _run_command(
            [sys.executable, str(script)],
            env=env,
            cwd=script.parent,
        )
    if code != 0:
        raise ToolRunError(f"{script.name} failed", stdout, stderr)

    after = list(output_dir.glob("*.csv"))
    new = _new_files(before, after)
    target_file: Optional[Path] = None
    if new:
        newest = max(new, key=lambda p: p.stat().st_mtime)
        name_hash = hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:12]
        dest_name = f"sim_{name_hash}_{newest.stem}.csv"
        dest_path = outdir_path / dest_name
        shutil.copy2(newest, dest_path)
        target_file = dest_path
        cache_target = _sim_cache_dir() / dest_name
        shutil.copy2(newest, cache_target)
    return target_file, stdout, stderr

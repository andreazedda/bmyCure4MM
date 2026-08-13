#!/usr/bin/env python3
"""Emit the minimum dependency identity evidence available before issue #19."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = (
    "Django",
    "numpy",
    "scipy",
    "pandas",
    "matplotlib",
    "rdkit",
)


def _git_sha() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "NOT_INSTALLED"


def build_report() -> dict[str, Any]:
    lock_path = ROOT / "uv.lock"
    return {
        "schema_version": "dependency-environment-v1",
        "git_sha": _git_sha(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "container_digest": os.environ.get("CONTAINER_IMAGE_DIGEST", "UNAVAILABLE_LOCAL"),
        "uv_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "packages": {name: _package_version(name) for name in PACKAGES},
    }


def main() -> int:
    print(json.dumps(build_report(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

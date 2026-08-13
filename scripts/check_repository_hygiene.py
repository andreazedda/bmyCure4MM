#!/usr/bin/env python3
"""Fail on tracked generated/private artifacts or new malformed notebooks."""

from __future__ import annotations

import ast
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = ROOT / ".quality-baseline.json"
FORBIDDEN_PARTS = {"__pycache__", "local_private", "media", "outputs"}
FORBIDDEN_SUFFIXES = {".db", ".pyc", ".sqlite3"}


def _tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return [Path(value.decode()) for value in result.stdout.split(b"\0") if value]


def _notebook_errors(path: Path) -> list[str]:
    try:
        notebook: dict[str, Any] = json.loads((ROOT / path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"invalid JSON: {exc}"]

    errors: list[str] = []
    for index, cell in enumerate(notebook.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if any(line.lstrip().startswith(("%", "!")) for line in source.splitlines()):
            continue
        try:
            ast.parse(source)
        except SyntaxError as exc:
            errors.append(f"cell {index}: {exc.msg} at line {exc.lineno}")
    return errors


def main() -> int:
    baseline = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    known_invalid: dict[str, str] = baseline["known_invalid_notebooks"]
    failures: list[str] = []

    for path in _tracked_files():
        if FORBIDDEN_PARTS.intersection(path.parts) or path.suffix in FORBIDDEN_SUFFIXES:
            failures.append(f"forbidden tracked artifact: {path}")
        if path.suffix != ".ipynb":
            continue
        errors = _notebook_errors(path)
        if not errors:
            if str(path) in known_invalid:
                failures.append(f"resolved notebook remains in invalid baseline: {path}")
            continue
        digest = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        if known_invalid.get(str(path)) != digest:
            failures.append(f"untriaged notebook {path} ({digest}): {'; '.join(errors)}")

    print(
        json.dumps(
            {
                "status": "PASS" if not failures else "FAIL",
                "known_invalid_notebooks": sorted(known_invalid),
                "failures": failures,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())

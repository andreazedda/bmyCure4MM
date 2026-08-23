from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from scripts.audit_dependencies import TRIAGED_TEST_TOOL_ADVISORIES


def test_governed_numerical_baseline_has_no_drift() -> None:
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "-m", "scripts.check_numerical_baseline"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_supported_django_baseline_has_no_framework_advisory_exceptions() -> None:
    root = Path(__file__).resolve().parents[1]
    audit_source = (root / "scripts" / "audit_dependencies.py").read_text(encoding="utf-8")

    assert TRIAGED_TEST_TOOL_ADVISORIES == ("GHSA-6w46-j5rx-g56g",)
    assert "DJANGO_42" not in audit_source
    assert "PYSEC-2026-3717" not in audit_source

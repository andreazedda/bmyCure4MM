from __future__ import annotations

import subprocess
import sys
from datetime import date
from pathlib import Path

from scripts.audit_dependencies import (
    TEMPORARY_DJANGO_42_EXCEPTION_EXPIRES,
    _temporary_django_exception_is_current,
)


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


def test_temporary_django_42_triage_has_a_hard_expiry() -> None:
    assert _temporary_django_exception_is_current(TEMPORARY_DJANGO_42_EXCEPTION_EXPIRES)
    assert not _temporary_django_exception_is_current(date(2026, 10, 1))

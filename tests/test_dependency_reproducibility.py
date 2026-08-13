from __future__ import annotations

import subprocess
import sys
from pathlib import Path


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

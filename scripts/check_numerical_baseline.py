#!/usr/bin/env python3
"""Compare deterministic scientific outputs with the governed dependency baseline."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import scipy

from simulator.design.reporting import run_design_report
from simulator.models_simulation import MathematicalModel

ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = ROOT / "tests" / "fixtures" / "numerical_baseline_v1.json"


def characterize() -> dict[str, Any]:
    model = MathematicalModel(
        baseline_tumor_cells=1e9,
        baseline_healthy_cells=5e11,
        drug_doses={"lenalidomide": 25.0, "bortezomib": 1.3},
        pk_params={
            "lenalidomide": {"half_life": 6.0, "Vd": 75.0},
            "bortezomib": {"half_life": 40.0, "Vd": 500.0},
        },
        pd_params={
            "lenalidomide": {"Emax": 0.12, "EC50": 0.5},
            "bortezomib": {"Emax": 0.2, "EC50": 0.2},
        },
        growth_rates={"tumor": 0.023, "healthy": 0.015},
        interaction_matrix=np.array([[0.0, 0.05], [0.05, 0.0]]),
        time_span=(0.0, 60.0),
        evaluation_points=200,
    )
    frame = model.simulate()
    report = asdict(run_design_report(seed=7, steps=18, dt_days=7.0))
    report_json = json.dumps(report, sort_keys=True, separators=(",", ":"))

    return {
        "environment": {
            "python": f"{sys.version_info.major}.{sys.version_info.minor}",
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "ode": {
            "final_tumor_burden": float(frame["tumor_cells"].iloc[-1]),
            "minimum_tumor_burden": float(frame["tumor_cells"].min()),
            "final_healthy_cells": float(frame["healthy_cells"].iloc[-1]),
            "minimum_healthy_cells": float(frame["healthy_cells"].min()),
            "final_lenalidomide_concentration": float(frame["lenalidomide_concentration"].iloc[-1]),
            "final_bortezomib_concentration": float(frame["bortezomib_concentration"].iloc[-1]),
            "lenalidomide_auc": float(
                np.trapezoid(  # type: ignore[attr-defined]
                    frame["lenalidomide_concentration"], frame["time"]
                )
            ),
            "bortezomib_auc": float(
                np.trapezoid(  # type: ignore[attr-defined]
                    frame["bortezomib_concentration"], frame["time"]
                )
            ),
        },
        "seeded_design_report": {
            "canonical_sha256": hashlib.sha256(report_json.encode()).hexdigest(),
            "canonical_bytes": len(report_json),
            "selected_therapy": report["rationale"]["selected_therapy"],
            "end_state": report["rationale"]["end_state"],
            "relapse_mean_days": report["relapse_forecast"]["mean_days"],
            "relapse_p50_days": report["relapse_forecast"]["p50_days"],
            "relapse_p90_days": report["relapse_forecast"]["p90_days"],
        },
    }


def _compare(expected: Any, actual: Any, path: str, differences: list[dict[str, Any]]) -> None:
    if isinstance(expected, dict) and isinstance(actual, dict):
        for key in sorted(expected.keys() | actual.keys()):
            child = f"{path}.{key}" if path else key
            if key not in expected or key not in actual:
                differences.append({"path": child, "expected": expected.get(key), "actual": actual.get(key)})
            else:
                _compare(expected[key], actual[key], child, differences)
        return
    if isinstance(expected, float) and isinstance(actual, (float, int)):
        if not math.isclose(expected, float(actual), rel_tol=1e-7, abs_tol=1e-10):
            differences.append({"path": path, "expected": expected, "actual": actual})
        return
    if expected != actual:
        differences.append({"path": path, "expected": expected, "actual": actual})


def main() -> int:
    baseline = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    actual = characterize()
    differences: list[dict[str, Any]] = []
    _compare(baseline["expected"], actual, "", differences)
    result = {
        "baseline_id": baseline["baseline_id"],
        "status": "PASS" if not differences else "NUMERICAL_DRIFT",
        "differences": differences,
        "actual": actual,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not differences else 1


if __name__ == "__main__":
    raise SystemExit(main())

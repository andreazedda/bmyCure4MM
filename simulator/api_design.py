from __future__ import annotations

from dataclasses import asdict

from django.http import JsonResponse
from django.views.decorators.http import require_GET

from simulator.design.reporting import run_design_report


def _parse_int(value: str | None, default: int, *, min_value: int, max_value: int) -> int:
    try:
        parsed = int(value) if value is not None else default
    except (TypeError, ValueError):
        return default
    return max(min_value, min(max_value, parsed))


def _parse_float(value: str | None, default: float, *, min_value: float, max_value: float) -> float:
    try:
        parsed = float(value) if value is not None else default
    except (TypeError, ValueError):
        return default
    return max(min_value, min(max_value, parsed))


@require_GET
def design_report(request):
    """Generate a design/simulation report JSON.

    Query params:
      - seed: int (default 7)
      - steps: int (default 18, max 120)
      - dt_days: float (default 7.0, max 30)

    This endpoint is intentionally side-effect-free (no DB writes).
    """

    seed = _parse_int(request.GET.get("seed"), 7, min_value=0, max_value=1_000_000)
    steps = _parse_int(request.GET.get("steps"), 18, min_value=1, max_value=120)
    dt_days = _parse_float(request.GET.get("dt_days"), 7.0, min_value=0.5, max_value=30.0)

    report = run_design_report(seed=seed, steps=steps, dt_days=dt_days)

    resp = JsonResponse(asdict(report))
    # This is deterministic for a given seed, but can be expensive.
    resp["Cache-Control"] = "max-age=60, private"
    return resp

from __future__ import annotations

import dataclasses
import json
import random
import statistics
from dataclasses import asdict
from pathlib import Path

from .simulator import step
from .targeting import estimate_tp_fp
from .toxicity import estimate_toxicity
from .types import AntigenProfile, Clone, DesignReport, LogicGate, NormalCompartment, Therapy, TumorState


def _json_default(obj):
    if dataclasses.is_dataclass(obj):
        return asdict(obj)
    raise TypeError(f"Not JSON serializable: {type(obj)}")


def _pareto_front(points: list[dict]) -> list[dict]:
    """Maximize TP, minimize FP."""
    pts = sorted(points, key=lambda p: (-p["tp"], p["fp"]))
    front: list[dict] = []
    best_fp = float("inf")
    for p in pts:
        if p["fp"] <= best_fp:
            front.append(p)
            best_fp = p["fp"]
    return front


def build_demo_scenario(seed: int = 7):
    rng = random.Random(seed)

    clones = {
        "C0": Clone(
            name="C0",
            antigen_profile=AntigenProfile({"BCMA": 120.0, "GPRC5D": 30.0, "SLAMF7": 80.0}),
            sensitivities={"small_molecule": 0.65, "adc": 0.75, "car_t": 0.85},
        ),
        "C1": Clone(
            name="C1",
            antigen_profile=AntigenProfile({"BCMA": 60.0, "GPRC5D": 90.0, "SLAMF7": 40.0}),
            sensitivities={"small_molecule": 0.55, "adc": 0.65, "car_t": 0.70},
        ),
        "C2": Clone(
            name="C2",
            antigen_profile=AntigenProfile({"BCMA": 15.0, "GPRC5D": 55.0, "SLAMF7": 20.0}),
            sensitivities={"small_molecule": 0.35, "adc": 0.45, "car_t": 0.50},
        ),
    }

    # Start with a bulk-dominant disease but non-trivial reservoir
    state = TumorState(
        time_days=0.0,
        loads={"bulk": 1.0, "reservoir": 0.35},
        fractions={
            "bulk": {"C0": 0.55, "C1": 0.35, "C2": 0.10},
            "reservoir": {"C0": 0.25, "C1": 0.35, "C2": 0.40},
        },
    )

    normals = [
        NormalCompartment(
            name="Marrow progenitors",
            antigen_profile=AntigenProfile({"BCMA": 5.0, "GPRC5D": 0.0, "SLAMF7": 15.0}),
            vulnerability=0.9,
        ),
        NormalCompartment(
            name="Plasma cells (normal)",
            antigen_profile=AntigenProfile({"BCMA": 40.0, "GPRC5D": 5.0, "SLAMF7": 30.0}),
            vulnerability=0.7,
        ),
    ]

    # Candidate therapies: single antigen vs logic-gated dual targeting
    therapies = []

    therapies.append(
        Therapy(
            name="ADC-BCMA",
            modality="adc",
            target="BCMA",
            logic=None,
            params={
                "name": "ADC-BCMA",
                "target_antigen": "BCMA",
                "affinity": 1.0,
                "internalization_rate": 0.75,
                "linker_cleavage": 0.70,
                "payload_potency": 0.75,
                "bystander_diffusion": 0.40,
                "systemic_leakage": 0.60,
            },
        )
    )

    therapies.append(
        Therapy(
            name="ADC-(BCMA AND GPRC5D)",
            modality="adc",
            target=None,
            logic=LogicGate(
                kind="AND",
                positive_antigens=["BCMA", "GPRC5D"],
                negative_antigens=[],
                thresholds={"BCMA": 50.0, "GPRC5D": 40.0},
                noise=0.0,
            ),
            params={
                "name": "ADC-(BCMA AND GPRC5D)",
                "target_antigen": "BCMA",
                "affinity": 1.0,
                "internalization_rate": 0.72,
                "linker_cleavage": 0.75,
                "payload_potency": 0.72,
                "bystander_diffusion": 0.35,
                "systemic_leakage": 0.55,
            },
        )
    )

    therapies.append(
        Therapy(
            name="CAR-T BCMA",
            modality="car_t",
            target="BCMA",
            logic=None,
            params={
                "name": "CAR-T BCMA",
                "target_antigen": "BCMA",
                "activation_threshold": 45.0,
                "expansion_rate": 0.65,
                "exhaustion_rate": 0.25,
                "cytokine_per_kill": 0.60,
            },
        )
    )

    return rng, clones, state, normals, therapies


def run_design_report(*, seed: int = 7, steps: int = 18, dt_days: float = 7.0) -> DesignReport:
    rng, clones, state, normals, therapies = build_demo_scenario(seed=seed)

    tumor_profiles = [c.antigen_profile for c in clones.values()]
    normal_profiles = [c.antigen_profile for c in normals]

    # Tradeoffs: TP/FP + toxicity per therapy
    tradeoff_points: list[dict] = []
    for therapy in therapies:
        target_antigen = str(therapy.params.get("target_antigen") or therapy.target or "")
        affinity = float(therapy.params.get("affinity", 1.0))
        tp_fp = estimate_tp_fp(
            tumor_profiles,
            normal_profiles,
            logic=therapy.logic,
            target_antigen=target_antigen if target_antigen else None,
            affinity=affinity,
        )
        tox_scores, _ = estimate_toxicity(therapy, normal_compartments=normals)
        tradeoff_points.append(
            {
                "therapy": therapy.name,
                "tp": float(tp_fp["tp"]),
                "fp": float(tp_fp["fp"]),
                "toxicity": {
                    "marrow_reserve_depletion": float(tox_scores["marrow_reserve_depletion"]),
                    "on_target_off_tumor": float(tox_scores["on_target_off_tumor"]),
                    "immune_mediated": float(tox_scores["immune_mediated"]),
                    "payload_mediated": float(tox_scores["payload_mediated"]),
                    "off_target_polypharmacology": float(tox_scores["off_target_polypharmacology"]),
                },
            }
        )

    pareto = _pareto_front([{k: v for k, v in p.items() if k in ("therapy", "tp", "fp")} for p in tradeoff_points])

    # Pick a therapy to simulate: best on Pareto (TP desc, FP asc)
    pareto_sorted = sorted(pareto, key=lambda p: (-p["tp"], p["fp"]))
    best_name = pareto_sorted[0]["therapy"] if pareto_sorted else therapies[0].name
    best_therapy = next(t for t in therapies if t.name == best_name)

    # Dynamics simulation
    outcomes = []
    for _ in range(steps):
        state, outcome, clones = step(
            state,
            clones=clones,
            therapy=best_therapy,
            normal_compartments=normals,
            dt_days=dt_days,
            rng=rng,
        )
        outcomes.append(asdict(outcome))

    # Relapse forecast: sample time-to-relapse from per-step relapse_risk
    relapse_times_days = []
    for _ in range(200):
        t = 0
        for o in outcomes:
            hazard = float(min(0.95, max(0.0, o.get("relapse_risk", 0.0))))
            if rng.random() < hazard * 0.25:
                break
            t += 1
        relapse_times_days.append(float(t * dt_days))

    relapse_forecast = {
        "samples_days": relapse_times_days,
        "mean_days": float(statistics.fmean(relapse_times_days)) if relapse_times_days else 0.0,
        "p50_days": float(statistics.median(relapse_times_days)) if relapse_times_days else 0.0,
        "p90_days": float(sorted(relapse_times_days)[int(0.9 * (len(relapse_times_days) - 1))]) if relapse_times_days else 0.0,
    }

    tox_scores, tox_why = estimate_toxicity(best_therapy, normal_compartments=normals)
    toxicity_decomposition = {
        "therapy": best_therapy.name,
        "scores": {k: float(v) for k, v in tox_scores.items()},
        "explanations": tox_why,
    }

    last = outcomes[-1] if outcomes else {}
    rationale = {
        "summary": "Pareto-optimal TP/FP selection with mechanistic toxicity + reservoir persistence considerations.",
        "selected_therapy": best_therapy.name,
        "end_state": {
            "time_days": float(last.get("time_days", 0.0)),
            "tumor_load_bulk": float(last.get("tumor_load_bulk", 0.0)),
            "tumor_load_reservoir": float(last.get("tumor_load_reservoir", 0.0)),
            "relapse_risk": float(last.get("relapse_risk", 0.0)),
        },
    }

    return DesignReport(
        inputs={
            "seed": seed,
            "steps": steps,
            "dt_days": dt_days,
            "scenario": "demo_mm_bulk_reservoir",
        },
        tradeoff_points=tradeoff_points,
        pareto_front=pareto,
        dynamics=outcomes,
        relapse_forecast=relapse_forecast,
        toxicity_decomposition=toxicity_decomposition,
        rationale=rationale,
    )


def write_report_json(report: DesignReport, output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(report), indent=2, default=_json_default) + "\n", encoding="utf-8")
    return path

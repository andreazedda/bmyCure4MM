from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple
import math
import optuna
from django.contrib.auth import get_user_model
from django.utils import timezone
from . import models
from .objectives import objectives, constraints
from .search_space import default_space

@dataclass
class TrialResult:
    params: Dict
    summary: Dict
    values: Tuple[float, float, float]  # efficacy, safety, exposure
    feasible: bool

def _suggest(trial: optuna.trial.Trial, space: Dict) -> Dict:
    params = {}
    for name, cfg in space.items():
        t = cfg["type"]
        if t == "float":
            params[name] = trial.suggest_float(name, cfg["low"], cfg["high"], step=cfg.get("step", None))
        elif t == "int":
            params[name] = trial.suggest_int(name, cfg["low"], cfg["high"], step=cfg.get("step", 1))
        else:
            raise ValueError(f"Unsupported type {t} for {name}")
    return params

def _apply_schedule(base: Dict, p: Dict) -> Dict:
    # Translate knobs into the model inputs; keep backward compatibility
    base = dict(base)
    base.update(
        lenalidomide_dose=p["lenalidomide_dose"],
        bortezomib_dose=p["bortezomib_dose"],
        daratumumab_dose=p["daratumumab_dose"],
        time_horizon=p["time_horizon"],
        interaction_strength=p["interaction_strength"],
    )
    # Optionally annotate schedule knobs for report
    base["_schedule"] = {
        "len_on_days": p["len_on_days"],
        "bor_weekly": p["bor_weekly"],
        "dara_interval": p["dara_interval"],
    }
    return base

def run_study(*, scenario: models.Scenario, user_id: int, n_trials: int = 100, seed: int = 42) -> Dict:
    space = default_space()
    pruner = optuna.pruners.MedianPruner(n_startup_trials=10, n_warmup_steps=0)
    # MOTPESampler is the multi-objective version of TPESampler
    sampler = optuna.samplers.TPESampler(seed=seed, multivariate=True)
    study = optuna.create_study(directions=["maximize", "maximize", "maximize"], sampler=sampler, pruner=pruner)

    base_params = {
        "baseline_tumor_cells": 1.0e9,
        "baseline_healthy_cells": 5.0e11,
        "tumor_growth_rate": 0.023,
        "healthy_growth_rate": 0.015,
    }

    results: List[TrialResult] = []

    def objective(trial: optuna.trial.Trial):
        p = _suggest(trial, space)
        model_params = _apply_schedule(base_params, p)
        attempt = models.SimulationAttempt.objects.create(
            scenario=scenario,
            user_id=user_id,
            parameters=model_params,
            submitted=timezone.now(),
        )
        summary = attempt.run_model()
        feas = all(constraints(summary).values())
        obj = objectives(summary)
        vals = (obj["efficacy"], obj["safety"], obj["exposure"])
        results.append(TrialResult(params=model_params, summary=summary, values=vals, feasible=feas))
        # Penalize infeasible by large negative
        if not feas:
            vals = tuple(v - 1e6 for v in vals)
        return vals

    study.optimize(objective, n_trials=n_trials, show_progress_bar=False)

    # Pareto front (non-dominated)
    trials = [r for r in results if r.feasible]
    def dominated(a, b):  # a dominated by b?
        return all(bv >= av for av, bv in zip(a.values, b.values)) and any(bv > av for av, bv in zip(a.values, b.values))
    pareto = []
    for r in trials:
        if not any(dominated(r, o) for o in trials):
            pareto.append(r)

    # Serialize compact payload for UI
    return {
        "pareto": [
            {
                "params": r.params,
                "efficacy": r.values[0],
                "safety": 1 + r.values[1],  # Convert back from -healthy_loss to safety score
                "exposure": -r.values[2],   # Convert back from negative exposure
            }
            for r in sorted(pareto, key=lambda x: (-x.values[0], x.values[1]))
        ],
        "study": study,
        "n_trials": len(results),
    }

from __future__ import annotations
import numpy as np
from typing import Dict, List
from . import models

def simulate_replicates(attempt: models.SimulationAttempt, n: int = 50, seed: int = 123) -> Dict:
    rng = np.random.default_rng(seed)
    base = dict(attempt.parameters)
    samples: List[Dict] = []
    for _ in range(n):
        perturbed = dict(base)
        # ±20% growth variability, ±10% pk half-life noise
        perturbed["tumor_growth_rate"] *= float(rng.normal(1.0, 0.1))
        perturbed["healthy_growth_rate"] *= float(rng.normal(1.0, 0.1))
        perturbed["interaction_strength"] *= max(0.0, float(rng.normal(1.0, 0.05)))
        a = models.SimulationAttempt.objects.create(
            scenario=attempt.scenario, user=attempt.user, parameters=perturbed
        )
        samples.append(a.run_model())
    # Summaries
    tr = np.array([s["tumor_reduction"] for s in samples])
    hl = np.array([s["healthy_loss"] for s in samples])
    return {
        "n": n,
        "efficacy_mean": float(tr.mean()), "efficacy_p05": float(np.quantile(tr, 0.05)), "efficacy_p95": float(np.quantile(tr, 0.95)),
        "safety_mean": float(hl.mean()), "safety_p05": float(np.quantile(hl, 0.05)), "safety_p95": float(np.quantile(hl, 0.95)),
    }

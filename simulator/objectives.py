from __future__ import annotations

from typing import Dict

# Map simulation output -> scalar objectives (higher is better unless negated)
def objectives(summary: Dict) -> Dict[str, float]:
    return {
        "efficacy": float(summary.get("tumor_reduction", 0.0)),           # maximize
        "safety": -float(summary.get("healthy_loss", 1.0)),               # minimize loss -> maximize negative
        "exposure": -float(sum((summary.get("auc") or {}).values())),     # minimize AUC sum
    }

# Hard constraints considered "feasible"
def constraints(summary: Dict) -> Dict[str, bool]:
    return {
        "toxicity_ok": float(summary.get("healthy_loss", 1.0)) <= 0.25,
    }

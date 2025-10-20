from __future__ import annotations
from typing import Dict, Any

# Discrete/continuous bounds for optimizer
def default_space() -> Dict[str, Any]:
    return {
        "lenalidomide_dose": {"type": "float", "low": 0.0, "high": 40.0, "step": 2.5},
        "bortezomib_dose": {"type": "float", "low": 0.0, "high": 1.6, "step": 0.1},
        "daratumumab_dose": {"type": "float", "low": 0.0, "high": 16.0, "step": 1.0},
        "time_horizon": {"type": "int", "low": 56, "high": 224, "step": 28},   # 2–8 cycles
        # Simple scheduling knobs
        "len_on_days": {"type": "int", "low": 14, "high": 28, "step": 7},      # days on per 28
        "bor_weekly": {"type": "int", "low": 0, "high": 2, "step": 1},         # 0/1/2 per week
        "dara_interval": {"type": "int", "low": 7, "high": 28, "step": 7},     # q1–4 weeks
        "interaction_strength": {"type": "float", "low": 0.0, "high": 0.15, "step": 0.01},
    }

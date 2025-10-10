"""Clinical regimen presets for simulator workflows."""

from __future__ import annotations

PRESETS: dict[str, dict[str, object]] = {
    "VRd": {
        "label": "VRd (Bortezomib + Lenalidomide + Dexamethasone)",
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 1.2e9,
            "baseline_healthy_cells": 4.5e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 0.0,
            "time_horizon": 168.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.08,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": [1, 4, 8, 11]},
            "daratumumab_dose": {"type": "days_list", "days": []},
        },
    },
    "Dara-Rd": {
        "label": "Daratumumab + Lenalidomide",
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 9.5e8,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 0.0,
            "daratumumab_dose": 16.0,
            "time_horizon": 196.0,
            "tumor_growth_rate": 0.02,
            "healthy_growth_rate": 0.016,
            "interaction_strength": 0.12,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": []},
            "daratumumab_dose": {"type": "qweekly", "weeks": [1, 2, 3, 4], "day": 1},
        },
    },
    "KRd": {
        "label": "Carfilzomib + Lenalidomide",
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 4.8e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.1,
            "daratumumab_dose": 0.0,
            "time_horizon": 168.0,
            "tumor_growth_rate": 0.024,
            "healthy_growth_rate": 0.014,
            "interaction_strength": 0.1,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": [1, 2, 8, 9, 15, 16]},
            "daratumumab_dose": {"type": "days_list", "days": []},
        },
    },
}

"""Hypothetical model-input presets for exploratory simulator workflows.

The numerical values and schedules below are model inputs, not patient-specific
instructions, clinical safety limits, or evidence of treatment benefit.
"""

from __future__ import annotations


def _story_en(name: str) -> dict[str, object]:
    return {
        "title": f"Hypothetical model preset: {name}",
        "background": (
            "This preset groups named drug-exposure inputs for an educational "
            "mechanistic scenario. It does not establish a clinical indication."
        ),
        "challenge": (
            "Explore how the current equations respond to the fixed assumptions, "
            "then inspect uncertainty and model-relative diagnostics."
        ),
        "why_these_drugs": (
            "Drug names map to configured PK/PD profiles. Their inclusion is a "
            "hypothetical model choice, not a treatment-selection claim."
        ),
        "expected_outcome": (
            "No clinical outcome is expected or predicted. The simulator emits "
            "model-relative compartment trajectories under stated assumptions."
        ),
        "learning_points": [
            "Schedule values are hypothetical inputs, not dosing instructions.",
            "Dose ranges are simulator bounds, not clinical safety limits.",
            "Outputs are SIMULATED and do not identify patient benefit or causality.",
            "Interpret comparisons only relative to the current model and version.",
        ],
    }


def _story_it(name: str) -> dict[str, object]:
    return {
        "title": f"Preset ipotetico del modello: {name}",
        "background": (
            "Questo preset raggruppa input di esposizione farmacologica per uno "
            "scenario meccanicistico educativo. Non stabilisce un'indicazione clinica."
        ),
        "challenge": (
            "Esplora la risposta delle equazioni alle ipotesi fissate, quindi "
            "esamina incertezza e diagnostica relativa al modello."
        ),
        "why_these_drugs": (
            "I nomi dei farmaci mappano profili PK/PD configurati. La loro presenza "
            "è una scelta ipotetica del modello, non una selezione terapeutica."
        ),
        "expected_outcome": (
            "Non viene previsto alcun esito clinico. Il simulatore produce traiettorie "
            "dei compartimenti relative al modello e alle ipotesi dichiarate."
        ),
        "learning_points": [
            "I calendari sono input ipotetici, non istruzioni di dosaggio.",
            "Gli intervalli di dose sono limiti del simulatore, non limiti di sicurezza clinica.",
            "Gli output sono SIMULATED e non identificano beneficio o causalità.",
            "Le comparazioni valgono solo rispetto al modello e alla versione correnti.",
        ],
    }


PRESETS: dict[str, dict[str, object]] = {
    "VRd": {
        "label": "VRd (Bortezomib + Lenalidomide + Dexamethasone)",
        "description_en": "Hypothetical named model-input preset; not a clinical selection.",
        "description_it": "Preset ipotetico di input del modello; non è una selezione clinica.",
        "story_en": _story_en("VRd"),
        "story_it": _story_it("VRd"),
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 1.2e9,
            "baseline_healthy_cells": 4.5e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 0.0,
            "carfilzomib_dose": 0.0,
            "time_horizon": 168.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.08,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": [1, 4, 8, 11]},
            "daratumumab_dose": {"type": "days_list", "days": []},
            "carfilzomib_dose": {"type": "days_list", "days": []},
        },
    },
    "Dara-Rd": {
        "label": "Daratumumab + Lenalidomide",
        "description_en": "Hypothetical named model-input preset; not a clinical selection.",
        "description_it": "Preset ipotetico di input del modello; non è una selezione clinica.",
        "story_en": _story_en("Dara-Rd"),
        "story_it": _story_it("Dara-Rd"),
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 9.5e8,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 0.0,
            "daratumumab_dose": 16.0,
            "carfilzomib_dose": 0.0,
            "time_horizon": 196.0,
            "tumor_growth_rate": 0.02,
            "healthy_growth_rate": 0.016,
            "interaction_strength": 0.12,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": []},
            "daratumumab_dose": {"type": "qweekly", "weeks": [1, 2, 3, 4], "day": 1},
            "carfilzomib_dose": {"type": "days_list", "days": []},
        },
    },
    "KRd": {
        "label": "Carfilzomib + Lenalidomide",
        "description_en": "Hypothetical named model-input preset; not a clinical selection.",
        "description_it": "Preset ipotetico di input del modello; non è una selezione clinica.",
        "story_en": _story_en("KRd"),
        "story_it": _story_it("KRd"),
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 4.8e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 0.0,
            "daratumumab_dose": 0.0,
            "carfilzomib_dose": 27.0,
            "time_horizon": 168.0,
            "tumor_growth_rate": 0.024,
            "healthy_growth_rate": 0.014,
            "interaction_strength": 0.1,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": []},
            "daratumumab_dose": {"type": "days_list", "days": []},
            "carfilzomib_dose": {"type": "days_list", "days": [1, 2, 8, 9, 15, 16]},
        },
    },
}

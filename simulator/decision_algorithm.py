"""Transparent model-relative diagnostic thresholds for research simulations.

The rules in this module classify model outputs for navigation and falsification.
They are not a prescribing algorithm and never emit a dose, duration, or regimen
instruction.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from mmportal.governance import EpistemicLabel, governance_metadata


ALGORITHM_VERSION = "2.0.0"
ALGORITHM_UPDATED = "2026-08-12"


DECISION_ALGORITHM: dict[str, Any] = {
    "version": ALGORITHM_VERSION,
    "last_updated": ALGORITHM_UPDATED,
    "governance": governance_metadata(
        epistemic_label=EpistemicLabel.HEURISTIC,
        output_kind="model_relative_threshold_catalog",
    ),
    "description": {
        "en": (
            "Transparent heuristic thresholds for classifying mechanistic "
            "simulation outputs. The catalog supports research navigation, "
            "not treatment selection."
        ),
        "it": (
            "Soglie euristiche trasparenti per classificare gli output di "
            "simulazioni meccanicistiche. Il catalogo supporta la navigazione "
            "della ricerca, non la scelta del trattamento."
        ),
    },
    "thresholds": {
        "efficacy": {
            "simulated_high_impact_zone": {
                "condition": "tumor_reduction >= 0.50",
                "description_en": "Model-relative reduction at or above 0.50",
                "description_it": "Riduzione relativa al modello pari o superiore a 0,50",
                "badge": "info",
            },
            "simulated_intermediate_zone": {
                "condition": "0.30 <= tumor_reduction < 0.50",
                "description_en": "Model-relative reduction between 0.30 and 0.50",
                "description_it": "Riduzione relativa al modello tra 0,30 e 0,50",
                "badge": "secondary",
            },
            "simulated_low_impact_zone": {
                "condition": "0 <= tumor_reduction < 0.30",
                "description_en": "Model-relative reduction below 0.30",
                "description_it": "Riduzione relativa al modello inferiore a 0,30",
                "badge": "warning",
            },
            "model_regrowth_signal": {
                "condition": "tumor_reduction < 0",
                "description_en": "Modeled endpoint above starting tumor burden",
                "description_it": "Endpoint modellato sopra il carico tumorale iniziale",
                "badge": "warning",
            },
        },
        "toxicity": {
            "below_constraint_watch_zone": {
                "condition": "healthy_loss < 0.20",
                "description_en": "Modeled healthy-cell loss below 0.20",
                "description_it": "Perdita modellata di cellule sane inferiore a 0,20",
                "badge": "secondary",
            },
            "constraint_watch_zone": {
                "condition": "0.20 <= healthy_loss < 0.30",
                "description_en": "Modeled healthy-cell loss between 0.20 and 0.30",
                "description_it": "Perdita modellata di cellule sane tra 0,20 e 0,30",
                "badge": "warning",
            },
            "constraint_flag": {
                "condition": "healthy_loss >= 0.30",
                "description_en": "Modeled healthy-cell loss at or above 0.30",
                "description_it": "Perdita modellata di cellule sane pari o superiore a 0,30",
                "badge": "warning",
            },
        },
        "durability": {
            "long_model_durability_signal": {
                "condition": "time_to_recurrence >= 180",
                "unit": "days",
                "description_en": "Modeled recurrence time at or above 180 days",
                "description_it": "Tempo di recidiva modellato pari o superiore a 180 giorni",
                "badge": "secondary",
            },
            "intermediate_model_durability_signal": {
                "condition": "90 <= time_to_recurrence < 180",
                "unit": "days",
                "description_en": "Modeled recurrence time between 90 and 180 days",
                "description_it": "Tempo di recidiva modellato tra 90 e 180 giorni",
                "badge": "secondary",
            },
            "short_model_durability_signal": {
                "condition": "time_to_recurrence < 90",
                "unit": "days",
                "description_en": "Modeled recurrence time below 90 days",
                "description_it": "Tempo di recidiva modellato inferiore a 90 giorni",
                "badge": "warning",
            },
        },
    },
    "decision_rules": [
        {
            "id": "RULE_001",
            "name_en": "Healthy-cell constraint flag",
            "name_it": "Segnale di vincolo sulle cellule sane",
            "trigger_condition": "healthy_loss >= 0.30",
            "priority": "high",
            "icon": "⚠️",
            "action_en": "Record constraint_flag and inspect model assumptions and uncertainty.",
            "action_it": "Registra constraint_flag e verifica assunzioni e incertezza del modello.",
            "rationale_en": "The simulated metric crosses a configured heuristic threshold.",
            "rationale_it": "La metrica simulata supera una soglia euristica configurata.",
            "evidence_status": "HEURISTIC_NEEDS_EVIDENCE",
        },
        {
            "id": "RULE_002",
            "name_en": "Healthy-cell constraint watch zone",
            "name_it": "Zona di attenzione del vincolo sulle cellule sane",
            "trigger_condition": "0.20 <= healthy_loss < 0.30",
            "priority": "medium",
            "icon": "⚠️",
            "action_en": "Record constraint_watch_zone and inspect sensitivity to assumptions.",
            "action_it": "Registra constraint_watch_zone e verifica la sensibilità alle assunzioni.",
            "rationale_en": "The simulated metric lies in a configured heuristic interval.",
            "rationale_it": "La metrica simulata si trova in un intervallo euristico configurato.",
            "evidence_status": "HEURISTIC_NEEDS_EVIDENCE",
        },
        {
            "id": "RULE_003",
            "name_en": "Simulated low-impact zone",
            "name_it": "Zona simulata a basso impatto",
            "trigger_condition": "tumor_reduction < 0.30 AND healthy_loss < 0.30",
            "priority": "high",
            "icon": "📉",
            "action_en": "Record simulated_low_impact_zone and compare assumptions or model variants.",
            "action_it": "Registra simulated_low_impact_zone e confronta assunzioni o varianti del modello.",
            "rationale_en": "The model-relative reduction is below a configured heuristic threshold.",
            "rationale_it": "La riduzione relativa al modello è sotto una soglia euristica configurata.",
            "evidence_status": "HEURISTIC_NEEDS_EVIDENCE",
        },
        {
            "id": "RULE_004",
            "name_en": "Model regrowth signal",
            "name_it": "Segnale di ricrescita del modello",
            "trigger_condition": "tumor_reduction < 0",
            "priority": "critical",
            "icon": "📈",
            "action_en": "Record model_regrowth_signal and test model/input falsification hypotheses.",
            "action_it": "Registra model_regrowth_signal e testa ipotesi di falsificazione del modello/input.",
            "rationale_en": "The modeled endpoint exceeds its starting tumor burden.",
            "rationale_it": "L'endpoint modellato supera il carico tumorale iniziale.",
            "evidence_status": "HEURISTIC_NEEDS_EVIDENCE",
        },
        {
            "id": "RULE_005",
            "name_en": "Short model-durability signal",
            "name_it": "Segnale di breve durabilità del modello",
            "trigger_condition": "time_to_recurrence < 90 AND time_horizon < 200",
            "priority": "medium",
            "icon": "⏱️",
            "action_en": "Record short_model_durability_signal and report horizon dependence.",
            "action_it": "Registra short_model_durability_signal e riporta la dipendenza dall'orizzonte.",
            "rationale_en": "The simulated recurrence time crosses a configured heuristic threshold.",
            "rationale_it": "Il tempo di recidiva simulato supera una soglia euristica configurata.",
            "evidence_status": "HEURISTIC_NEEDS_EVIDENCE",
        },
        {
            "id": "RULE_006",
            "name_en": "Simulated high-impact zone without constraint flag",
            "name_it": "Zona simulata ad alto impatto senza segnale di vincolo",
            "trigger_condition": "tumor_reduction >= 0.50 AND healthy_loss < 0.20",
            "priority": "low",
            "icon": "ℹ️",
            "action_en": "Record simulated_high_impact_zone and test robustness across uncertainty.",
            "action_it": "Registra simulated_high_impact_zone e testa la robustezza rispetto all'incertezza.",
            "rationale_en": "Two model-relative metrics fall within configured heuristic zones.",
            "rationale_it": "Due metriche relative al modello rientrano in zone euristiche configurate.",
            "evidence_status": "HEURISTIC_NEEDS_EVIDENCE",
        },
    ],
    "risk_stratification": {
        "R-ISS_I": {
            "criteria": "β2M < 3.5 mg/L AND albumin ≥ 3.5 g/dL AND standard-risk cytogenetics AND normal LDH",
            "description_en": "Population-level R-ISS I context; not an individual outcome prediction.",
            "description_it": "Contesto R-ISS I di popolazione; non una previsione individuale.",
            "epistemic_label": EpistemicLabel.LITERATURE_INFORMED.value,
        },
        "R-ISS_II": {
            "criteria": "Not R-ISS I or III",
            "description_en": "Population-level R-ISS II context; not an individual outcome prediction.",
            "description_it": "Contesto R-ISS II di popolazione; non una previsione individuale.",
            "epistemic_label": EpistemicLabel.LITERATURE_INFORMED.value,
        },
        "R-ISS_III": {
            "criteria": "β2M ≥ 5.5 mg/L AND (high-risk cytogenetics OR elevated LDH)",
            "description_en": "Population-level R-ISS III context; not an individual outcome prediction.",
            "description_it": "Contesto R-ISS III di popolazione; non una previsione individuale.",
            "epistemic_label": EpistemicLabel.LITERATURE_INFORMED.value,
        },
    },
    "high_risk_cytogenetics": {
        "status": "literature_context_only",
        "epistemic_label": EpistemicLabel.LITERATURE_INFORMED.value,
        "limitation": "No treatment-selection inference is authorized by this catalog.",
    },
    "data_sources": [],
    "disclaimer": {
        "en": (
            "Research prototype only. Thresholds classify mechanistic model outputs; "
            "they do not identify a causal effect, validate a patient-specific prediction, "
            "or authorize a treatment change."
        ),
        "it": (
            "Solo prototipo di ricerca. Le soglie classificano output di un modello "
            "meccanicistico; non identificano un effetto causale, non validano una "
            "previsione individuale e non autorizzano modifiche terapeutiche."
        ),
    },
}


def get_algorithm() -> dict[str, Any]:
    return deepcopy(DECISION_ALGORITHM)


def get_threshold(category: str, level: str) -> dict[str, Any] | None:
    return deepcopy(DECISION_ALGORITHM.get("thresholds", {}).get(category, {}).get(level))


def get_rule(rule_id: str) -> dict[str, Any] | None:
    for rule in DECISION_ALGORITHM.get("decision_rules", []):
        if rule.get("id") == rule_id:
            return deepcopy(rule)
    return None


def get_risk_stratification(r_iss: str) -> dict[str, Any] | None:
    key = f"R-ISS_{r_iss.upper()}" if r_iss else None
    if key:
        return deepcopy(DECISION_ALGORITHM.get("risk_stratification", {}).get(key))
    return None


def evaluate_metrics(
    tumor_reduction: float | None,
    healthy_loss: float | None,
    time_to_recurrence: float | None,
) -> dict[str, str]:
    result: dict[str, str] = {}
    if tumor_reduction is None:
        result["efficacy"] = "unknown"
    elif tumor_reduction >= 0.50:
        result["efficacy"] = "simulated_high_impact_zone"
    elif tumor_reduction >= 0.30:
        result["efficacy"] = "simulated_intermediate_zone"
    elif tumor_reduction >= 0:
        result["efficacy"] = "simulated_low_impact_zone"
    else:
        result["efficacy"] = "model_regrowth_signal"

    if healthy_loss is None:
        result["toxicity"] = "unknown"
    elif healthy_loss < 0.20:
        result["toxicity"] = "below_constraint_watch_zone"
    elif healthy_loss < 0.30:
        result["toxicity"] = "constraint_watch_zone"
    else:
        result["toxicity"] = "constraint_flag"

    if time_to_recurrence is None:
        result["durability"] = "unknown"
    elif time_to_recurrence >= 180:
        result["durability"] = "long_model_durability_signal"
    elif time_to_recurrence >= 90:
        result["durability"] = "intermediate_model_durability_signal"
    else:
        result["durability"] = "short_model_durability_signal"
    return result


def get_applicable_rules(
    tumor_reduction: float | None,
    healthy_loss: float | None,
    time_to_recurrence: float | None,
    time_horizon: int = 168,
) -> list[dict[str, Any]]:
    applicable = []
    if healthy_loss is not None and healthy_loss >= 0.30:
        applicable.append(get_rule("RULE_001"))
    elif healthy_loss is not None and 0.20 <= healthy_loss < 0.30:
        applicable.append(get_rule("RULE_002"))
    if tumor_reduction is not None and tumor_reduction < 0.30:
        if healthy_loss is None or healthy_loss < 0.30:
            applicable.append(get_rule("RULE_003"))
    if tumor_reduction is not None and tumor_reduction < 0:
        applicable.append(get_rule("RULE_004"))
    if time_to_recurrence is not None and time_to_recurrence < 90 and time_horizon < 200:
        applicable.append(get_rule("RULE_005"))
    if (
        tumor_reduction is not None
        and tumor_reduction >= 0.50
        and healthy_loss is not None
        and healthy_loss < 0.20
    ):
        applicable.append(get_rule("RULE_006"))
    priority_order = {"critical": 0, "high": 1, "medium": 2, "low": 3}
    rules = [rule for rule in applicable if rule is not None]
    rules.sort(key=lambda rule: priority_order.get(rule.get("priority", "low"), 99))
    return rules

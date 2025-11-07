"""Bilingual KPI tooltips for simulator interfaces."""

from __future__ import annotations

KPI: dict[str, dict[str, str]] = {
    "tumor_reduction": {
        "en": "Tumor Reduction = 1 − (Tumor_end / Tumor_start). Higher is better.",
        "it": "Riduzione Tumorale = 1 − (Tumore_fine / Tumore_inizio). Più alto è meglio.",
    },
    "healthy_loss": {
        "en": "Healthy Loss = (Healthy_start − Healthy_end) / Healthy_start. Keep ≤ 0.25.",
        "it": "Perdita Sani = (Sani_inizio − Sani_fine) / Sani_inizio. Mantieni ≤ 0.25.",
    },
    "auc": {
        "en": "AUC captures total drug exposure over time. Lower exposure reduces toxicity risk.",
        "it": "L’AUC rappresenta l’esposizione totale al farmaco nel tempo. Ridurla limita la tossicità.",
    },
    "durability_index": {
        "en": "Durability = time below baseline ÷ time horizon. Measures response persistence.",
        "it": "Durabilità = tempo sotto la baseline ÷ orizzonte. Misura la persistenza della risposta.",
    },
    "time_to_recurrence": {
        "en": "Time to recurrence marks when tumor burden rises above 50% of baseline.",
        "it": "Il Tempo alla recidiva indica quando la massa tumorale torna oltre il 50% della baseline.",
    },
}


def tip(key: str) -> dict[str, str]:
    """Return KPI copy for tooltips or fall back to empty strings."""
    return KPI.get(key, {"en": "", "it": ""})

from __future__ import annotations

from django import template

from simulator.permissions import is_editor

register = template.Library()


@register.simple_tag
def simulator_is_editor(user) -> bool:
    """Return True when the given user can manage simulator content."""
    return is_editor(user)


@register.filter
def dict_get(mapping, key):
    """Safe dictionary lookup in templates."""
    if isinstance(mapping, dict):
        return mapping.get(key, "")
    return ""


@register.filter
def italian_field_label(field):
    """Translate visible learning labels; preserve scientific field names and units."""
    labels = {
        "baseline_tumor_cells": "Cellule tumorali iniziali",
        "baseline_healthy_cells": "Plasmacellule sane iniziali",
        "lenalidomide_dose": "Dose giornaliera di lenalidomide (mg)",
        "bortezomib_dose": "Dose settimanale di bortezomib (mg/m²)",
        "daratumumab_dose": "Dose di daratumumab (mg/kg)",
        "carfilzomib_dose": "Dose di carfilzomib (mg/m²)",
        "time_horizon": "Orizzonte della simulazione (giorni)",
        "cohort_size": "Dimensione della coorte virtuale",
        "tumor_growth_rate": "Tasso di crescita tumorale (giorno⁻¹)",
        "healthy_growth_rate": "Tasso di crescita delle cellule sane (giorno⁻¹)",
        "interaction_strength": "Intensità dell’interazione tra farmaci",
        "creatinine_clearance": "Clearance della creatinina (ml/min)",
        "neuropathy_grade": "Grado di neuropatia",
        "anc": "Conta assoluta dei neutrofili (×10⁹/L)",
        "platelets": "Piastrine (×10⁹/L)",
        "seed": "Seme casuale",
        "carrying_capacity_tumor": "Capacità portante tumorale (cellule)",
        "carrying_capacity_healthy": "Capacità portante delle cellule sane (cellule)",
        "immune_compromise_index": "Indice di compromissione immunitaria",
    }
    return labels.get(field.name, field.label)

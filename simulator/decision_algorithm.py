"""
Transparent Decision Algorithm for MM Treatment Optimization.

This module exposes the decision logic in a machine-readable format,
making the algorithm fully transparent and auditable.
"""

from __future__ import annotations

from typing import Any

# Version and metadata
ALGORITHM_VERSION = "1.0.0"
ALGORITHM_UPDATED = "2026-01-17"

# Complete algorithm definition
DECISION_ALGORITHM: dict[str, Any] = {
    "version": ALGORITHM_VERSION,
    "last_updated": ALGORITHM_UPDATED,
    "description": {
        "en": "Transparent decision algorithm for Multiple Myeloma treatment optimization. "
              "All thresholds and rules are evidence-based and fully disclosed.",
        "it": "Algoritmo decisionale trasparente per ottimizzazione trattamento Mieloma Multiplo. "
              "Tutte le soglie e regole sono basate sull'evidenza e completamente divulgate.",
    },
    
    "thresholds": {
        "efficacy": {
            "good": {
                "condition": "tumor_reduction >= 0.50",
                "description_en": "≥50% tumor reduction - Good response",
                "description_it": "≥50% riduzione tumorale - Buona risposta",
                "badge": "success",
            },
            "moderate": {
                "condition": "0.30 <= tumor_reduction < 0.50",
                "description_en": "30-50% tumor reduction - Moderate response",
                "description_it": "30-50% riduzione tumorale - Risposta moderata",
                "badge": "warning",
            },
            "poor": {
                "condition": "0 <= tumor_reduction < 0.30",
                "description_en": "<30% tumor reduction - Poor response",
                "description_it": "<30% riduzione tumorale - Risposta scarsa",
                "badge": "danger",
            },
            "failure": {
                "condition": "tumor_reduction < 0",
                "description_en": "Tumor growth - Treatment failure",
                "description_it": "Crescita tumorale - Fallimento terapeutico",
                "badge": "dark",
            },
        },
        "toxicity": {
            "acceptable": {
                "condition": "healthy_loss < 0.20",
                "description_en": "<20% healthy cell loss - Acceptable toxicity",
                "description_it": "<20% perdita cellule sane - Tossicità accettabile",
                "badge": "success",
            },
            "borderline": {
                "condition": "0.20 <= healthy_loss < 0.30",
                "description_en": "20-30% healthy cell loss - Borderline toxicity",
                "description_it": "20-30% perdita cellule sane - Tossicità borderline",
                "badge": "warning",
            },
            "high": {
                "condition": "healthy_loss >= 0.30",
                "description_en": "≥30% healthy cell loss - High toxicity",
                "description_it": "≥30% perdita cellule sane - Alta tossicità",
                "badge": "danger",
            },
        },
        "durability": {
            "good": {
                "condition": "time_to_recurrence >= 180",
                "unit": "days",
                "description_en": "≥6 months to recurrence - Good durability",
                "description_it": "≥6 mesi a recidiva - Buona durabilità",
                "badge": "success",
            },
            "moderate": {
                "condition": "90 <= time_to_recurrence < 180",
                "unit": "days",
                "description_en": "3-6 months to recurrence - Moderate durability",
                "description_it": "3-6 mesi a recidiva - Durabilità moderata",
                "badge": "warning",
            },
            "poor": {
                "condition": "time_to_recurrence < 90",
                "unit": "days",
                "description_en": "<3 months to recurrence - Poor durability",
                "description_it": "<3 mesi a recidiva - Scarsa durabilità",
                "badge": "danger",
            },
        },
    },
    
    "decision_rules": [
        {
            "id": "RULE_001",
            "name_en": "High Toxicity Management",
            "name_it": "Gestione Alta Tossicità",
            "trigger_condition": "healthy_loss >= 0.30",
            "priority": "high",
            "icon": "⚠️",
            "action_en": "Reduce drug doses by 20-30% or shorten treatment duration",
            "action_it": "Riduci dosi farmaci del 20-30% o accorcia durata trattamento",
            "rationale_en": "Excessive damage to healthy plasma cells compromises immune function "
                           "and increases infection risk. Dose reduction preserves efficacy while "
                           "improving tolerability.",
            "rationale_it": "Danno eccessivo alle plasmacellule sane compromette funzione immunitaria "
                           "e aumenta rischio infezioni. Riduzione dose preserva efficacia migliorando "
                           "tollerabilità.",
            "evidence": [
                "NCCN Multiple Myeloma Guidelines v2.2024",
                "Palumbo A et al. Blood 2014 (dose adjustments)",
                "IFM recommendations for frail patients",
            ],
        },
        {
            "id": "RULE_002",
            "name_en": "Moderate Toxicity Monitoring",
            "name_it": "Monitoraggio Tossicità Moderata",
            "trigger_condition": "0.20 <= healthy_loss < 0.30",
            "priority": "medium",
            "icon": "👁️",
            "action_en": "Consider reducing doses by 10-15% if patient shows clinical toxicity signs",
            "action_it": "Considera riduzione dosi del 10-15% se paziente mostra segni clinici di tossicità",
            "rationale_en": "Borderline toxicity requires close monitoring. Preemptive dose reduction "
                           "may prevent treatment discontinuation.",
            "rationale_it": "Tossicità borderline richiede monitoraggio stretto. Riduzione preventiva "
                           "può prevenire sospensione trattamento.",
            "evidence": [
                "IMWG guidelines for toxicity management",
                "Clinical experience from MAIA/CASSIOPEIA trials",
            ],
        },
        {
            "id": "RULE_003",
            "name_en": "Poor Efficacy Escalation",
            "name_it": "Escalation per Scarsa Efficacia",
            "trigger_condition": "tumor_reduction < 0.30 AND healthy_loss < 0.30",
            "priority": "high",
            "icon": "📈",
            "action_en": "Increase doses by 15-25% or extend treatment to 224-280 days",
            "action_it": "Aumenta dosi del 15-25% o estendi trattamento a 224-280 giorni",
            "rationale_en": "Suboptimal response with tolerable toxicity allows room for intensification. "
                           "Deeper responses correlate with longer progression-free survival.",
            "rationale_it": "Risposta subottimale con tossicità tollerabile permette intensificazione. "
                           "Risposte più profonde correlano con sopravvivenza libera da progressione più lunga.",
            "evidence": [
                "IMWG response criteria 2016",
                "MAIA trial (deeper responses = better PFS)",
                "CASSIOPEIA trial (MRD negativity benefit)",
            ],
        },
        {
            "id": "RULE_004",
            "name_en": "Treatment Failure - Regimen Switch",
            "name_it": "Fallimento Terapeutico - Cambio Regime",
            "trigger_condition": "tumor_reduction < 0",
            "priority": "critical",
            "icon": "🚨",
            "action_en": "Switch to alternative regimen immediately. Consider: DPd, KRd, Isa-Pd, or BCMA-targeted therapy",
            "action_it": "Cambia regime terapeutico immediatamente. Considera: DPd, KRd, Isa-Pd, o terapia anti-BCMA",
            "rationale_en": "Tumor growth indicates primary resistance. Continuing ineffective therapy "
                           "wastes time and exposes patient to toxicity without benefit.",
            "rationale_it": "Crescita tumorale indica resistenza primaria. Continuare terapia inefficace "
                           "spreca tempo ed espone paziente a tossicità senza beneficio.",
            "evidence": [
                "IMWG progressive disease criteria",
                "NCCN salvage therapy recommendations",
                "ICARIA-MM, IKEMA trials (isatuximab)",
                "KarMMa, CARTITUDE trials (BCMA-targeted)",
            ],
            "suggested_alternatives": [
                {"name": "DPd", "indication": "daratumumab-refractory excluded"},
                {"name": "KRd", "indication": "bortezomib-exposed"},
                {"name": "Isa-Pd", "indication": "daratumumab-refractory"},
                {"name": "Ide-cel / Cilta-cel", "indication": "triple-class exposed"},
            ],
        },
        {
            "id": "RULE_005",
            "name_en": "Early Recurrence Prevention",
            "name_it": "Prevenzione Recidiva Precoce",
            "trigger_condition": "time_to_recurrence < 90 AND time_horizon < 200",
            "priority": "medium",
            "icon": "⏱️",
            "action_en": "Extend treatment horizon to 224-280 days. Consider maintenance therapy.",
            "action_it": "Estendi orizzonte a 224-280 giorni. Considera terapia di mantenimento.",
            "rationale_en": "Longer treatment duration delays recurrence. Maintenance therapy with "
                           "lenalidomide significantly improves progression-free survival.",
            "rationale_it": "Durata maggiore ritarda recidiva. Mantenimento con lenalidomide "
                           "migliora significativamente sopravvivenza libera da progressione.",
            "evidence": [
                "FIRST trial (lenalidomide continuous vs fixed)",
                "TOURMALINE-MM3 (ixazomib maintenance)",
                "Myeloma XI (lenalidomide maintenance)",
            ],
        },
        {
            "id": "RULE_006",
            "name_en": "Optimal Balance - Fine Tuning",
            "name_it": "Equilibrio Ottimale - Ottimizzazione Fine",
            "trigger_condition": "tumor_reduction >= 0.50 AND healthy_loss < 0.20",
            "priority": "low",
            "icon": "✅",
            "action_en": "Current regimen is optimal. Consider ±10% dose adjustments or compare alternative regimens.",
            "action_it": "Regime attuale ottimale. Considera variazioni ±10% dosi o confronta regimi alternativi.",
            "rationale_en": "Good efficacy with acceptable toxicity indicates optimal balance. "
                           "Minor adjustments may further personalize therapy.",
            "rationale_it": "Buona efficacia con tossicità accettabile indica equilibrio ottimale. "
                           "Piccole modifiche possono personalizzare ulteriormente la terapia.",
            "evidence": [
                "Treat-to-target principles in MM",
                "Personalized medicine approaches",
            ],
        },
    ],
    
    "risk_stratification": {
        "R-ISS_I": {
            "criteria": "β2M < 3.5 mg/L AND albumin ≥ 3.5 g/dL AND standard-risk cytogenetics AND normal LDH",
            "median_PFS_months": 66,
            "median_OS_months": None,  # Not reached in most studies
            "five_year_survival_percent": 82,
            "description_en": "Low risk - Excellent prognosis. Standard triplet therapy recommended.",
            "description_it": "Basso rischio - Prognosi eccellente. Raccomandata terapia triplet standard.",
        },
        "R-ISS_II": {
            "criteria": "Not R-ISS I or III",
            "median_PFS_months": 42,
            "median_OS_months": 83,
            "five_year_survival_percent": 62,
            "description_en": "Intermediate risk - Standard prognosis. Consider quadruplet for fit patients.",
            "description_it": "Rischio intermedio - Prognosi standard. Considera quadruplet per pazienti fit.",
        },
        "R-ISS_III": {
            "criteria": "β2M ≥ 5.5 mg/L AND (high-risk cytogenetics OR elevated LDH)",
            "median_PFS_months": 29,
            "median_OS_months": 43,
            "five_year_survival_percent": 40,
            "description_en": "High risk - Aggressive disease. Intensive therapy + early transplant if eligible.",
            "description_it": "Alto rischio - Malattia aggressiva. Terapia intensiva + trapianto precoce se eleggibile.",
        },
    },
    
    "high_risk_cytogenetics": {
        "del17p": {
            "name": "del(17p)",
            "gene": "TP53",
            "impact_en": "Tumor suppressor loss. 40% reduction in PFS.",
            "impact_it": "Perdita oncosoppressore. 40% riduzione PFS.",
            "management_en": "Consider VRd + daratumumab quadruplet. Early transplant.",
            "management_it": "Considera quadruplet VRd + daratumumab. Trapianto precoce.",
        },
        "t_4_14": {
            "name": "t(4;14)",
            "gene": "FGFR3/MMSET",
            "impact_en": "25% reduction in PFS. Sensitive to proteasome inhibitors.",
            "impact_it": "25% riduzione PFS. Sensibile a inibitori proteasoma.",
            "management_en": "Bortezomib-based regimen preferred.",
            "management_it": "Regime basato su bortezomib preferito.",
        },
        "t_14_16": {
            "name": "t(14;16)",
            "gene": "MAF",
            "impact_en": "50% reduction in PFS. Very high risk.",
            "impact_it": "50% riduzione PFS. Rischio molto alto.",
            "management_en": "Aggressive quadruplet therapy. Clinical trial if available.",
            "management_it": "Terapia quadruplet aggressiva. Trial clinico se disponibile.",
        },
        "gain_1q21": {
            "name": "1q21 gain/amplification",
            "gene": "CKS1B",
            "impact_en": "15% reduction in PFS. Proliferation advantage.",
            "impact_it": "15% riduzione PFS. Vantaggio proliferativo.",
            "management_en": "May benefit from carfilzomib-based regimens.",
            "management_it": "Può beneficiare da regimi basati su carfilzomib.",
        },
    },
    
    "data_sources": [
        {
            "name": "R-ISS Staging",
            "reference": "Palumbo A, et al. J Clin Oncol. 2015;33(26):2863-9",
            "doi": "10.1200/JCO.2015.61.2267",
        },
        {
            "name": "MAIA Trial",
            "reference": "Facon T, et al. N Engl J Med. 2019;380(22):2104-2115",
            "doi": "10.1056/NEJMoa1817249",
        },
        {
            "name": "CASSIOPEIA Trial",
            "reference": "Moreau P, et al. Lancet. 2019;394(10192):29-38",
            "doi": "10.1016/S0140-6736(19)31240-1",
        },
        {
            "name": "NCCN Guidelines",
            "reference": "NCCN Clinical Practice Guidelines in Oncology: Multiple Myeloma v2.2024",
            "url": "https://www.nccn.org/guidelines/guidelines-detail?category=1&id=1445",
        },
    ],
    
    "disclaimer": {
        "en": "This algorithm provides decision support based on clinical trial data and guidelines. "
             "It is NOT a substitute for clinical judgment. Individual patient factors, "
             "comorbidities, preferences, and local protocols must be considered. "
             "Always validate recommendations with treating physician.",
        "it": "Questo algoritmo fornisce supporto decisionale basato su dati di trial clinici e linee guida. "
             "NON sostituisce il giudizio clinico. Fattori individuali del paziente, "
             "comorbidità, preferenze e protocolli locali devono essere considerati. "
             "Validare sempre le raccomandazioni con il medico curante.",
    },
}


def get_algorithm() -> dict[str, Any]:
    """Return the complete decision algorithm."""
    return DECISION_ALGORITHM


def get_threshold(category: str, level: str) -> dict[str, Any] | None:
    """Get a specific threshold definition."""
    return DECISION_ALGORITHM.get("thresholds", {}).get(category, {}).get(level)


def get_rule(rule_id: str) -> dict[str, Any] | None:
    """Get a specific decision rule by ID."""
    for rule in DECISION_ALGORITHM.get("decision_rules", []):
        if rule.get("id") == rule_id:
            return rule
    return None


def get_risk_stratification(r_iss: str) -> dict[str, Any] | None:
    """Get risk stratification data for an R-ISS stage."""
    key = f"R-ISS_{r_iss.upper()}" if r_iss else None
    if key:
        return DECISION_ALGORITHM.get("risk_stratification", {}).get(key)
    return None


def evaluate_metrics(
    tumor_reduction: float | None,
    healthy_loss: float | None,
    time_to_recurrence: float | None,
) -> dict[str, str]:
    """
    Evaluate simulation metrics against algorithm thresholds.
    
    Returns dict with 'efficacy', 'toxicity', 'durability' labels.
    """
    result = {}
    
    # Efficacy
    if tumor_reduction is None:
        result["efficacy"] = "unknown"
    elif tumor_reduction >= 0.50:
        result["efficacy"] = "good"
    elif tumor_reduction >= 0.30:
        result["efficacy"] = "moderate"
    elif tumor_reduction >= 0:
        result["efficacy"] = "poor"
    else:
        result["efficacy"] = "failure"
    
    # Toxicity
    if healthy_loss is None:
        result["toxicity"] = "unknown"
    elif healthy_loss < 0.20:
        result["toxicity"] = "acceptable"
    elif healthy_loss < 0.30:
        result["toxicity"] = "borderline"
    else:
        result["toxicity"] = "high"
    
    # Durability
    if time_to_recurrence is None:
        result["durability"] = "unknown"
    elif time_to_recurrence >= 180:
        result["durability"] = "good"
    elif time_to_recurrence >= 90:
        result["durability"] = "moderate"
    else:
        result["durability"] = "poor"
    
    return result


def get_applicable_rules(
    tumor_reduction: float | None,
    healthy_loss: float | None,
    time_to_recurrence: float | None,
    time_horizon: int = 168,
) -> list[dict[str, Any]]:
    """
    Get decision rules applicable to given simulation results.
    
    Returns list of applicable rules sorted by priority.
    """
    applicable = []
    
    # Rule 001: High toxicity
    if healthy_loss is not None and healthy_loss >= 0.30:
        applicable.append(get_rule("RULE_001"))
    
    # Rule 002: Moderate toxicity
    elif healthy_loss is not None and 0.20 <= healthy_loss < 0.30:
        applicable.append(get_rule("RULE_002"))
    
    # Rule 003: Poor efficacy (with tolerable toxicity)
    if tumor_reduction is not None and tumor_reduction < 0.30:
        if healthy_loss is None or healthy_loss < 0.30:
            applicable.append(get_rule("RULE_003"))
    
    # Rule 004: Treatment failure
    if tumor_reduction is not None and tumor_reduction < 0:
        applicable.append(get_rule("RULE_004"))
    
    # Rule 005: Early recurrence
    if time_to_recurrence is not None and time_to_recurrence < 90:
        if time_horizon < 200:
            applicable.append(get_rule("RULE_005"))
    
    # Rule 006: Optimal balance
    if (tumor_reduction is not None and tumor_reduction >= 0.50 and
        healthy_loss is not None and healthy_loss < 0.20):
        applicable.append(get_rule("RULE_006"))
    
    # Sort by priority
    priority_order = {"critical": 0, "high": 1, "medium": 2, "low": 3}
    applicable = [r for r in applicable if r is not None]
    applicable.sort(key=lambda r: priority_order.get(r.get("priority", "low"), 99))
    
    return applicable

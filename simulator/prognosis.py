"""
Prognosis estimation module for Multiple Myeloma patients.

Provides survival estimates based on:
- R-ISS staging (I, II, III)
- High-risk cytogenetics [t(4;14), t(14;16), del(17p), 1q gain, etc.]
- Patient characteristics (age, ECOG, comorbidities)
- Treatment response

All estimates are based on published literature and clinical trials.
References are provided for transparency.

DISCLAIMER: These are statistical estimates for educational purposes only.
Individual patient outcomes may vary significantly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal
import math


@dataclass
class PrognosisEstimate:
    """Container for prognosis estimates with confidence intervals."""
    
    median_pfs_months: float  # Median Progression-Free Survival
    median_os_months: float   # Median Overall Survival
    
    # Survival probabilities at key timepoints
    pfs_12m: float  # Probability of PFS at 12 months
    pfs_24m: float  # Probability of PFS at 24 months
    pfs_36m: float  # Probability of PFS at 36 months
    
    os_12m: float   # Probability of OS at 12 months
    os_24m: float   # Probability of OS at 24 months
    os_36m: float   # Probability of OS at 36 months
    os_60m: float   # Probability of OS at 5 years
    
    # Risk category
    risk_category: Literal["standard", "intermediate", "high", "very_high"]
    
    # Confidence/reliability score (0-1)
    confidence: float
    
    # Modifiers applied
    modifiers_applied: list[str]
    
    # Primary reference
    reference: str
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "median_pfs_months": self.median_pfs_months,
            "median_os_months": self.median_os_months,
            "survival_probabilities": {
                "pfs": {
                    "12_months": round(self.pfs_12m, 2),
                    "24_months": round(self.pfs_24m, 2),
                    "36_months": round(self.pfs_36m, 2),
                },
                "os": {
                    "12_months": round(self.os_12m, 2),
                    "24_months": round(self.os_24m, 2),
                    "36_months": round(self.os_36m, 2),
                    "60_months": round(self.os_60m, 2),
                },
            },
            "risk_category": self.risk_category,
            "confidence": round(self.confidence, 2),
            "modifiers_applied": self.modifiers_applied,
            "reference": self.reference,
        }


# ══════════════════════════════════════════════════════════════════════════════
# BASELINE SURVIVAL DATA (R-ISS based)
# Source: Palumbo et al. JCO 2015; IMWG consensus
# ══════════════════════════════════════════════════════════════════════════════

BASELINE_BY_RISS = {
    "I": {
        "median_pfs": 66,      # months
        "median_os": 120,      # months (not reached in many studies, estimated)
        "pfs_12m": 0.90,
        "pfs_24m": 0.80,
        "pfs_36m": 0.70,
        "os_12m": 0.97,
        "os_24m": 0.93,
        "os_36m": 0.88,
        "os_60m": 0.75,
        "risk": "standard",
        "reference": "Palumbo et al. JCO 2015; 5-year OS ~82%",
    },
    "II": {
        "median_pfs": 42,
        "median_os": 83,
        "pfs_12m": 0.82,
        "pfs_24m": 0.65,
        "pfs_36m": 0.52,
        "os_12m": 0.93,
        "os_24m": 0.83,
        "os_36m": 0.72,
        "os_60m": 0.55,
        "risk": "intermediate",
        "reference": "Palumbo et al. JCO 2015; 5-year OS ~62%",
    },
    "III": {
        "median_pfs": 29,
        "median_os": 43,
        "pfs_12m": 0.70,
        "pfs_24m": 0.48,
        "pfs_36m": 0.35,
        "os_12m": 0.82,
        "os_24m": 0.62,
        "os_36m": 0.47,
        "os_60m": 0.28,
        "risk": "high",
        "reference": "Palumbo et al. JCO 2015; 5-year OS ~40%",
    },
}

# ══════════════════════════════════════════════════════════════════════════════
# CYTOGENETIC RISK MODIFIERS
# Source: IMWG 2016 consensus, MAYO mSMART 3.0
# ══════════════════════════════════════════════════════════════════════════════

CYTOGENETIC_HAZARD_RATIOS = {
    # High-risk cytogenetics (HR for PFS and OS)
    "del(17p)": {"pfs_hr": 2.0, "os_hr": 2.5, "category": "high"},
    "t(4;14)": {"pfs_hr": 1.8, "os_hr": 1.8, "category": "high"},
    "t(14;16)": {"pfs_hr": 2.0, "os_hr": 2.2, "category": "high"},
    "t(14;20)": {"pfs_hr": 1.8, "os_hr": 2.0, "category": "high"},
    "1q21_gain": {"pfs_hr": 1.5, "os_hr": 1.6, "category": "intermediate"},
    "1p_deletion": {"pfs_hr": 1.4, "os_hr": 1.5, "category": "intermediate"},
    
    # Standard risk (no significant impact)
    "t(11;14)": {"pfs_hr": 1.0, "os_hr": 0.9, "category": "standard"},  # May be favorable
    "hyperdiploidy": {"pfs_hr": 0.9, "os_hr": 0.85, "category": "standard"},
    
    # Double-hit / Ultra-high risk combinations
    "double_hit": {"pfs_hr": 3.0, "os_hr": 4.0, "category": "very_high"},
}


# ══════════════════════════════════════════════════════════════════════════════
# AGE AND ECOG MODIFIERS
# ══════════════════════════════════════════════════════════════════════════════

def get_age_modifier(age: int) -> tuple[float, float]:
    """
    Returns (PFS multiplier, OS multiplier) based on age.
    Older age generally associated with worse outcomes, but also
    often excludes intensive therapy.
    """
    if age < 65:
        return (1.0, 1.0)  # Reference group
    elif age < 75:
        return (0.95, 0.90)  # Slight decrease
    elif age < 80:
        return (0.85, 0.75)  # Moderate decrease
    else:
        return (0.75, 0.60)  # Significant decrease


def get_ecog_modifier(ecog: int) -> tuple[float, float]:
    """
    Returns (PFS multiplier, OS multiplier) based on ECOG status.
    Poor performance status strongly associated with worse outcomes.
    """
    modifiers = {
        0: (1.0, 1.0),      # Fully active
        1: (0.95, 0.95),    # Restricted strenuous activity
        2: (0.85, 0.80),    # Ambulatory, unable to work
        3: (0.70, 0.55),    # Limited self-care
        4: (0.50, 0.30),    # Completely disabled
    }
    return modifiers.get(ecog, (1.0, 1.0))


# ══════════════════════════════════════════════════════════════════════════════
# TREATMENT RESPONSE MODIFIERS
# ══════════════════════════════════════════════════════════════════════════════

RESPONSE_MODIFIERS = {
    "sCR": {"pfs_mult": 1.3, "os_mult": 1.4, "description": "Stringent Complete Response"},
    "CR": {"pfs_mult": 1.2, "os_mult": 1.3, "description": "Complete Response"},
    "VGPR": {"pfs_mult": 1.1, "os_mult": 1.15, "description": "Very Good Partial Response"},
    "PR": {"pfs_mult": 1.0, "os_mult": 1.0, "description": "Partial Response (reference)"},
    "SD": {"pfs_mult": 0.7, "os_mult": 0.8, "description": "Stable Disease"},
    "PD": {"pfs_mult": 0.4, "os_mult": 0.5, "description": "Progressive Disease"},
}

# MRD status (Minimal Residual Disease)
MRD_MODIFIERS = {
    "negative": {"pfs_mult": 1.5, "os_mult": 1.4, "description": "MRD Negative (<10⁻⁵)"},
    "positive": {"pfs_mult": 1.0, "os_mult": 1.0, "description": "MRD Positive"},
    "unknown": {"pfs_mult": 1.0, "os_mult": 1.0, "description": "MRD Not Tested"},
}


# ══════════════════════════════════════════════════════════════════════════════
# MAIN PROGNOSIS FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def estimate_prognosis(
    r_iss: Literal["I", "II", "III"] = "II",
    cytogenetics: list[str] | None = None,
    age: int | None = None,
    ecog: int | None = None,
    response: str | None = None,
    mrd_status: str | None = None,
    line_of_therapy: int = 1,
) -> PrognosisEstimate:
    """
    Estimate prognosis for a MM patient based on multiple factors.
    
    Args:
        r_iss: R-ISS stage (I, II, or III)
        cytogenetics: List of cytogenetic abnormalities (codes)
        age: Patient age in years
        ecog: ECOG performance status (0-4)
        response: Treatment response (sCR, CR, VGPR, PR, SD, PD)
        mrd_status: Minimal residual disease status ("negative", "positive", "unknown")
        line_of_therapy: Current line of therapy (1=frontline, 2+=relapsed)
    
    Returns:
        PrognosisEstimate with survival predictions and confidence
    """
    cytogenetics = cytogenetics or []
    modifiers_applied = []
    
    # Start with baseline R-ISS values
    baseline = BASELINE_BY_RISS.get(r_iss, BASELINE_BY_RISS["II"])
    
    median_pfs = baseline["median_pfs"]
    median_os = baseline["median_os"]
    
    pfs_mult = 1.0
    os_mult = 1.0
    
    risk_level = baseline["risk"]
    
    # ── Apply cytogenetic modifiers ─────────────────────────────────────────
    for cyto in cytogenetics:
        cyto_lower = cyto.lower().replace("-", "_").replace("(", "").replace(")", "")
        
        # Try to match known cytogenetics
        matched = False
        for known, data in CYTOGENETIC_HAZARD_RATIOS.items():
            if known.lower().replace("(", "").replace(")", "") in cyto_lower or cyto_lower in known.lower():
                pfs_mult /= data["pfs_hr"]  # HR > 1 means worse prognosis
                os_mult /= data["os_hr"]
                modifiers_applied.append(f"Cytogenetics: {cyto} (HR {data['pfs_hr']})")
                
                # Update risk level
                if data["category"] == "very_high":
                    risk_level = "very_high"
                elif data["category"] == "high" and risk_level not in ["very_high"]:
                    risk_level = "high"
                elif data["category"] == "intermediate" and risk_level == "standard":
                    risk_level = "intermediate"
                    
                matched = True
                break
        
        if not matched and cyto:
            modifiers_applied.append(f"Cytogenetics: {cyto} (unknown impact)")
    
    # Check for double-hit (≥2 high-risk cytogenetics)
    high_risk_cyto_count = sum(
        1 for c in cytogenetics 
        for known, data in CYTOGENETIC_HAZARD_RATIOS.items()
        if data["category"] in ["high", "very_high"] and known.lower() in c.lower()
    )
    if high_risk_cyto_count >= 2:
        pfs_mult *= 0.6  # Additional penalty for double-hit
        os_mult *= 0.5
        risk_level = "very_high"
        modifiers_applied.append("Double-hit myeloma (≥2 high-risk cytogenetics)")
    
    # ── Apply age modifier ──────────────────────────────────────────────────
    if age is not None:
        age_pfs, age_os = get_age_modifier(age)
        pfs_mult *= age_pfs
        os_mult *= age_os
        if age_pfs < 1.0:
            modifiers_applied.append(f"Age {age} years (adjustment)")
    
    # ── Apply ECOG modifier ─────────────────────────────────────────────────
    if ecog is not None:
        ecog_pfs, ecog_os = get_ecog_modifier(ecog)
        pfs_mult *= ecog_pfs
        os_mult *= ecog_os
        if ecog >= 2:
            modifiers_applied.append(f"ECOG {ecog} (significant impact)")
        elif ecog == 1:
            modifiers_applied.append(f"ECOG {ecog} (minor adjustment)")
    
    # ── Apply response modifier ─────────────────────────────────────────────
    if response and response in RESPONSE_MODIFIERS:
        resp_data = RESPONSE_MODIFIERS[response]
        pfs_mult *= resp_data["pfs_mult"]
        os_mult *= resp_data["os_mult"]
        modifiers_applied.append(f"Response: {resp_data['description']}")
    
    # ── Apply MRD modifier ──────────────────────────────────────────────────
    if mrd_status and mrd_status in MRD_MODIFIERS:
        mrd_data = MRD_MODIFIERS[mrd_status]
        pfs_mult *= mrd_data["pfs_mult"]
        os_mult *= mrd_data["os_mult"]
        if mrd_status == "negative":
            modifiers_applied.append("MRD Negative (favorable)")
    
    # ── Apply line of therapy adjustment ────────────────────────────────────
    if line_of_therapy > 1:
        # Each relapse shortens survival expectations
        lot_factor = 0.7 ** (line_of_therapy - 1)
        pfs_mult *= lot_factor
        os_mult *= max(lot_factor, 0.5)  # Floor at 50% for OS
        modifiers_applied.append(f"Line {line_of_therapy} therapy (relapsed)")
    
    # ── Calculate adjusted values ───────────────────────────────────────────
    adjusted_pfs = median_pfs * pfs_mult
    adjusted_os = median_os * os_mult
    
    # Calculate survival probabilities using exponential decay approximation
    # S(t) = exp(-λt) where λ = ln(2) / median
    def survival_prob(median: float, months: float) -> float:
        if median <= 0:
            return 0.0
        lambda_rate = math.log(2) / median
        return math.exp(-lambda_rate * months)
    
    # Calculate adjusted probabilities
    pfs_12m = survival_prob(adjusted_pfs, 12)
    pfs_24m = survival_prob(adjusted_pfs, 24)
    pfs_36m = survival_prob(adjusted_pfs, 36)
    
    os_12m = survival_prob(adjusted_os, 12)
    os_24m = survival_prob(adjusted_os, 24)
    os_36m = survival_prob(adjusted_os, 36)
    os_60m = survival_prob(adjusted_os, 60)
    
    # ── Calculate confidence score ──────────────────────────────────────────
    # Higher confidence when more data points are provided
    confidence_factors = [
        0.5,  # Base confidence
        0.15 if r_iss else 0,
        0.10 if cytogenetics else 0,
        0.08 if age is not None else 0,
        0.07 if ecog is not None else 0,
        0.10 if response else 0,
    ]
    confidence = min(sum(confidence_factors), 0.95)
    
    return PrognosisEstimate(
        median_pfs_months=round(adjusted_pfs, 1),
        median_os_months=round(adjusted_os, 1),
        pfs_12m=pfs_12m,
        pfs_24m=pfs_24m,
        pfs_36m=pfs_36m,
        os_12m=os_12m,
        os_24m=os_24m,
        os_36m=os_36m,
        os_60m=os_60m,
        risk_category=risk_level,
        confidence=confidence,
        modifiers_applied=modifiers_applied,
        reference=baseline["reference"],
    )


def get_prognosis_explanation(estimate: PrognosisEstimate, lang: str = "en") -> dict:
    """
    Generate human-readable explanation of prognosis estimate.
    
    Args:
        estimate: PrognosisEstimate object
        lang: Language code ("en" or "it")
    
    Returns:
        Dictionary with explanation text
    """
    if lang == "it":
        risk_labels = {
            "standard": "Standard",
            "intermediate": "Intermedio",
            "high": "Alto",
            "very_high": "Molto Alto",
        }
        
        return {
            "summary": f"Rischio {risk_labels[estimate.risk_category]}",
            "pfs_text": f"Sopravvivenza libera da progressione mediana stimata: {estimate.median_pfs_months:.0f} mesi",
            "os_text": f"Sopravvivenza globale mediana stimata: {estimate.median_os_months:.0f} mesi",
            "timeline": {
                "3_months": {
                    "pfs": f"A 3 mesi: {survival_prob_text(estimate.pfs_12m * 1.3, 'it')} senza progressione",
                    "os": f"A 3 mesi: {survival_prob_text(estimate.os_12m * 1.2, 'it')} sopravvivenza",
                },
                "6_months": {
                    "pfs": f"A 6 mesi: {survival_prob_text(estimate.pfs_12m * 1.15, 'it')} senza progressione",
                    "os": f"A 6 mesi: {survival_prob_text(estimate.os_12m * 1.1, 'it')} sopravvivenza",
                },
                "12_months": {
                    "pfs": f"A 12 mesi: {survival_prob_text(estimate.pfs_12m, 'it')} senza progressione",
                    "os": f"A 12 mesi: {survival_prob_text(estimate.os_12m, 'it')} sopravvivenza",
                },
                "24_months": {
                    "pfs": f"A 24 mesi: {survival_prob_text(estimate.pfs_24m, 'it')} senza progressione",
                    "os": f"A 24 mesi: {survival_prob_text(estimate.os_24m, 'it')} sopravvivenza",
                },
            },
            "confidence": f"Affidabilità stima: {estimate.confidence * 100:.0f}%",
            "disclaimer": "Queste sono stime statistiche basate su dati di popolazione. L'esito individuale può variare significativamente.",
        }
    else:
        risk_labels = {
            "standard": "Standard",
            "intermediate": "Intermediate",
            "high": "High",
            "very_high": "Very High",
        }
        
        return {
            "summary": f"{risk_labels[estimate.risk_category]} Risk",
            "pfs_text": f"Estimated median progression-free survival: {estimate.median_pfs_months:.0f} months",
            "os_text": f"Estimated median overall survival: {estimate.median_os_months:.0f} months",
            "timeline": {
                "3_months": {
                    "pfs": f"At 3 months: {survival_prob_text(estimate.pfs_12m * 1.3, 'en')} progression-free",
                    "os": f"At 3 months: {survival_prob_text(estimate.os_12m * 1.2, 'en')} survival",
                },
                "6_months": {
                    "pfs": f"At 6 months: {survival_prob_text(estimate.pfs_12m * 1.15, 'en')} progression-free",
                    "os": f"At 6 months: {survival_prob_text(estimate.os_12m * 1.1, 'en')} survival",
                },
                "12_months": {
                    "pfs": f"At 12 months: {survival_prob_text(estimate.pfs_12m, 'en')} progression-free",
                    "os": f"At 12 months: {survival_prob_text(estimate.os_12m, 'en')} survival",
                },
                "24_months": {
                    "pfs": f"At 24 months: {survival_prob_text(estimate.pfs_24m, 'en')} progression-free",
                    "os": f"At 24 months: {survival_prob_text(estimate.os_24m, 'en')} survival",
                },
            },
            "confidence": f"Estimate confidence: {estimate.confidence * 100:.0f}%",
            "disclaimer": "These are statistical estimates based on population data. Individual outcomes may vary significantly.",
        }


def survival_prob_text(prob: float, lang: str) -> str:
    """Convert probability to human-readable text."""
    prob = min(prob, 0.99)  # Cap at 99%
    pct = int(prob * 100)
    
    if lang == "it":
        if pct >= 90:
            return f"~{pct}% probabilità"
        elif pct >= 70:
            return f"~{pct}% probabilità"
        elif pct >= 50:
            return f"~{pct}% probabilità"
        else:
            return f"~{pct}% probabilità"
    else:
        if pct >= 90:
            return f"~{pct}% probability"
        elif pct >= 70:
            return f"~{pct}% probability"
        elif pct >= 50:
            return f"~{pct}% probability"
        else:
            return f"~{pct}% probability"


# ══════════════════════════════════════════════════════════════════════════════
# PROGNOSIS COMPARISON FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def compare_scenarios(
    base_params: dict,
    scenarios: list[dict],
) -> list[dict]:
    """
    Compare prognosis across different treatment scenarios.
    
    Args:
        base_params: Base patient parameters (r_iss, cytogenetics, age, ecog)
        scenarios: List of scenario dicts with response, mrd_status, line_of_therapy
    
    Returns:
        List of comparison results with differences
    """
    base_estimate = estimate_prognosis(**base_params)
    
    results = []
    for i, scenario in enumerate(scenarios):
        scenario_params = {
            key: value for key, value in scenario.items() if key != "name"
        }
        combined_params = {**base_params, **scenario_params}
        scenario_estimate = estimate_prognosis(**combined_params)
        
        pfs_diff = scenario_estimate.median_pfs_months - base_estimate.median_pfs_months
        os_diff = scenario_estimate.median_os_months - base_estimate.median_os_months
        
        results.append({
            "scenario_id": i + 1,
            "scenario_name": scenario.get("name", f"Scenario {i + 1}"),
            "params": scenario_params,
            "estimate": scenario_estimate.to_dict(),
            "vs_baseline": {
                "pfs_diff_months": round(pfs_diff, 1),
                "os_diff_months": round(os_diff, 1),
                "pfs_improvement": f"+{pfs_diff:.0f} months" if pfs_diff > 0 else f"{pfs_diff:.0f} months",
                "os_improvement": f"+{os_diff:.0f} months" if os_diff > 0 else f"{os_diff:.0f} months",
            },
        })
    
    return results

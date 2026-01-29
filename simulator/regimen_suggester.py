"""
Regimen Suggester Module for Multiple Myeloma.

Suggests appropriate treatment regimens based on:
- Patient characteristics (age, frailty, comorbidities)
- Disease characteristics (R-ISS, cytogenetics, tumor burden)
- Treatment history (line of therapy, prior exposures)
- Current response status

Based on IMWG guidelines, NCCN recommendations, and mSMART criteria.

DISCLAIMER: These suggestions are for educational purposes only.
Treatment decisions must be made by qualified healthcare professionals.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass
class RegimenSuggestion:
    """A suggested treatment regimen with supporting information."""
    
    name: str                           # Regimen abbreviation (e.g., VRd, DRd)
    full_name: str                      # Full regimen name
    components: list[str]               # Individual drugs
    indication: str                     # When to use
    strength: Literal["strong", "moderate", "conditional"]  # Recommendation strength
    evidence_level: str                 # Evidence quality
    key_trials: list[str]              # Supporting clinical trials
    considerations: list[str]          # Important notes
    contraindications: list[str]       # When to avoid
    expected_response_rate: str        # ORR/CR rates if known
    
    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "full_name": self.full_name,
            "components": self.components,
            "indication": self.indication,
            "strength": self.strength,
            "evidence_level": self.evidence_level,
            "key_trials": self.key_trials,
            "considerations": self.considerations,
            "contraindications": self.contraindications,
            "expected_response_rate": self.expected_response_rate,
        }


# ══════════════════════════════════════════════════════════════════════════════
# REGIMEN DATABASE
# ══════════════════════════════════════════════════════════════════════════════

FRONTLINE_REGIMENS = {
    # Transplant-eligible
    "VRd": RegimenSuggestion(
        name="VRd",
        full_name="Bortezomib, Lenalidomide, Dexamethasone",
        components=["Bortezomib", "Lenalidomide", "Dexamethasone"],
        indication="Transplant-eligible NDMM",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["SWOG S0777", "IFM 2009"],
        considerations=[
            "Standard of care for transplant-eligible patients",
            "Weekly bortezomib reduces neuropathy risk",
            "Continue lenalidomide maintenance post-transplant"
        ],
        contraindications=["Severe neuropathy", "Renal failure requiring dialysis (adjust)"],
        expected_response_rate="ORR ~90%, ≥VGPR ~70%, CR ~30%"
    ),
    
    "Dara-VRd": RegimenSuggestion(
        name="Dara-VRd",
        full_name="Daratumumab + Bortezomib, Lenalidomide, Dexamethasone",
        components=["Daratumumab", "Bortezomib", "Lenalidomide", "Dexamethasone"],
        indication="Transplant-eligible NDMM",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["GRIFFIN", "PERSEUS"],
        considerations=[
            "Deepest responses in transplant-eligible patients",
            "Higher MRD negativity rates",
            "Infusion reactions with first dose (use SC formulation if available)"
        ],
        contraindications=["Prior daratumumab exposure"],
        expected_response_rate="ORR ~99%, ≥VGPR ~90%, sCR ~65%"
    ),
    
    "KRd": RegimenSuggestion(
        name="KRd",
        full_name="Carfilzomib, Lenalidomide, Dexamethasone",
        components=["Carfilzomib", "Lenalidomide", "Dexamethasone"],
        indication="Transplant-eligible NDMM (alternative to VRd)",
        strength="moderate",
        evidence_level="Category 2A",
        key_trials=["ENDURANCE"],
        considerations=[
            "Alternative for patients with pre-existing neuropathy",
            "Monitor for cardiac toxicity",
            "Requires IV administration"
        ],
        contraindications=["Significant cardiac disease", "Uncontrolled hypertension"],
        expected_response_rate="ORR ~87%, ≥VGPR ~70%"
    ),
    
    # Transplant-ineligible
    "DRd": RegimenSuggestion(
        name="DRd",
        full_name="Daratumumab, Lenalidomide, Dexamethasone",
        components=["Daratumumab", "Lenalidomide", "Dexamethasone"],
        indication="Transplant-ineligible NDMM",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["MAIA"],
        considerations=[
            "Preferred regimen for transplant-ineligible patients",
            "Continue until progression",
            "Excellent long-term outcomes"
        ],
        contraindications=["Prior daratumumab exposure", "Severe renal impairment (adjust Rd)"],
        expected_response_rate="ORR ~93%, ≥VGPR ~79%, CR ~48%"
    ),
    
    "VMP": RegimenSuggestion(
        name="VMP",
        full_name="Bortezomib, Melphalan, Prednisone",
        components=["Bortezomib", "Melphalan", "Prednisone"],
        indication="Transplant-ineligible NDMM, limited access settings",
        strength="moderate",
        evidence_level="Category 1",
        key_trials=["VISTA"],
        considerations=[
            "Fixed-duration treatment (9 cycles)",
            "Good option when lenalidomide not available",
            "Weekly bortezomib preferred"
        ],
        contraindications=["Severe cytopenias", "Severe neuropathy"],
        expected_response_rate="ORR ~71%, CR ~30%"
    ),
    
    "Dara-VMP": RegimenSuggestion(
        name="Dara-VMP",
        full_name="Daratumumab + Bortezomib, Melphalan, Prednisone",
        components=["Daratumumab", "Bortezomib", "Melphalan", "Prednisone"],
        indication="Transplant-ineligible NDMM",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["ALCYONE"],
        considerations=[
            "Fixed duration induction followed by daratumumab maintenance",
            "Good for patients who cannot tolerate continuous lenalidomide"
        ],
        contraindications=["Prior daratumumab exposure", "Severe cytopenias"],
        expected_response_rate="ORR ~91%, ≥VGPR ~72%, CR ~43%"
    ),
    
    "Rd": RegimenSuggestion(
        name="Rd",
        full_name="Lenalidomide, Dexamethasone",
        components=["Lenalidomide", "Dexamethasone"],
        indication="Frail elderly or patients with contraindications to other agents",
        strength="moderate",
        evidence_level="Category 1",
        key_trials=["FIRST"],
        considerations=[
            "Well-tolerated in frail patients",
            "Continuous until progression",
            "Reduce dexamethasone dose in elderly"
        ],
        contraindications=["Renal impairment requiring dialysis (dose adjust)"],
        expected_response_rate="ORR ~75%, ≥VGPR ~44%"
    ),
}

RELAPSED_REGIMENS = {
    "DPd": RegimenSuggestion(
        name="DPd",
        full_name="Daratumumab, Pomalidomide, Dexamethasone",
        components=["Daratumumab", "Pomalidomide", "Dexamethasone"],
        indication="Relapsed MM after lenalidomide-based therapy",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["APOLLO"],
        considerations=[
            "Preferred for lenalidomide-refractory patients",
            "Can be used for daratumumab-naïve or re-treatment"
        ],
        contraindications=["Pomalidomide allergy"],
        expected_response_rate="ORR ~69%, ≥VGPR ~51%"
    ),
    
    "DVd": RegimenSuggestion(
        name="DVd",
        full_name="Daratumumab, Bortezomib, Dexamethasone",
        components=["Daratumumab", "Bortezomib", "Dexamethasone"],
        indication="First relapse, especially if lenalidomide-refractory",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["CASTOR"],
        considerations=[
            "Good for patients relapsing on lenalidomide maintenance",
            "Weekly bortezomib reduces neuropathy"
        ],
        contraindications=["Severe neuropathy", "Prior bortezomib intolerance"],
        expected_response_rate="ORR ~83%, ≥VGPR ~63%"
    ),
    
    "KPd": RegimenSuggestion(
        name="KPd",
        full_name="Carfilzomib, Pomalidomide, Dexamethasone",
        components=["Carfilzomib", "Pomalidomide", "Dexamethasone"],
        indication="Relapsed MM, particularly bortezomib-exposed",
        strength="moderate",
        evidence_level="Category 2A",
        key_trials=["EMN011"],
        considerations=[
            "Good option for bortezomib-refractory patients",
            "Monitor cardiac function"
        ],
        contraindications=["Significant cardiac disease", "Severe pulmonary disease"],
        expected_response_rate="ORR ~87%, ≥VGPR ~69%"
    ),
    
    "IsaPd": RegimenSuggestion(
        name="Isa-Pd",
        full_name="Isatuximab, Pomalidomide, Dexamethasone",
        components=["Isatuximab", "Pomalidomide", "Dexamethasone"],
        indication="Relapsed/refractory MM (2+ prior lines)",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["ICARIA-MM"],
        considerations=[
            "Anti-CD38 alternative to daratumumab",
            "May have activity in some daratumumab-refractory patients"
        ],
        contraindications=["Isatuximab allergy"],
        expected_response_rate="ORR ~60%, ≥VGPR ~32%"
    ),
    
    "Selinexor-Vd": RegimenSuggestion(
        name="SVd",
        full_name="Selinexor, Bortezomib, Dexamethasone",
        components=["Selinexor", "Bortezomib", "Dexamethasone"],
        indication="Relapsed/refractory MM (1-3 prior lines)",
        strength="moderate",
        evidence_level="Category 1",
        key_trials=["BOSTON"],
        considerations=[
            "Once-weekly selinexor and bortezomib",
            "Manage GI toxicity with supportive care"
        ],
        contraindications=["Severe thrombocytopenia", "Active uncontrolled infection"],
        expected_response_rate="ORR ~76%, ≥VGPR ~44%"
    ),
    
    "Elotuzumab-Pd": RegimenSuggestion(
        name="Elo-Pd",
        full_name="Elotuzumab, Pomalidomide, Dexamethasone",
        components=["Elotuzumab", "Pomalidomide", "Dexamethasone"],
        indication="Relapsed MM (2+ prior lines)",
        strength="moderate",
        evidence_level="Category 1",
        key_trials=["ELOQUENT-3"],
        considerations=[
            "SLAMF7-targeted therapy",
            "Generally well-tolerated"
        ],
        contraindications=["Elotuzumab allergy"],
        expected_response_rate="ORR ~53%, ≥VGPR ~20%"
    ),
}

NOVEL_THERAPIES = {
    "Ide-cel": RegimenSuggestion(
        name="Ide-cel",
        full_name="Idecabtagene vicleucel (CAR-T)",
        components=["Anti-BCMA CAR-T cells"],
        indication="Relapsed/refractory MM (4+ prior lines, including anti-CD38)",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["KarMMa", "KarMMa-3"],
        considerations=[
            "One-time infusion",
            "Requires specialized center",
            "CRS and neurotoxicity risk",
            "Bridging therapy often needed"
        ],
        contraindications=["Active infection", "Poor performance status", "Organ dysfunction"],
        expected_response_rate="ORR ~73%, CR ~33%"
    ),
    
    "Cilta-cel": RegimenSuggestion(
        name="Cilta-cel",
        full_name="Ciltacabtagene autoleucel (CAR-T)",
        components=["Anti-BCMA CAR-T cells"],
        indication="Relapsed/refractory MM (1+ prior lines) or triple-class exposed",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["CARTITUDE-1", "CARTITUDE-4"],
        considerations=[
            "High response rates even in heavily pretreated",
            "Movement and neurocognitive toxicity risk",
            "Manufacturing time required"
        ],
        contraindications=["Active infection", "Poor performance status"],
        expected_response_rate="ORR ~98%, sCR ~83%"
    ),
    
    "Teclistamab": RegimenSuggestion(
        name="Teclistamab",
        full_name="Teclistamab (BCMA bispecific)",
        components=["Teclistamab"],
        indication="Relapsed/refractory MM (4+ prior lines)",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["MajesTEC-1"],
        considerations=[
            "Off-the-shelf bispecific antibody",
            "Step-up dosing required",
            "CRS risk, monitor closely"
        ],
        contraindications=["Active infection", "Severe cytopenias"],
        expected_response_rate="ORR ~63%, ≥VGPR ~59%"
    ),
    
    "Elranatamab": RegimenSuggestion(
        name="Elranatamab",
        full_name="Elranatamab (BCMA bispecific)",
        components=["Elranatamab"],
        indication="Relapsed/refractory MM (4+ prior lines)",
        strength="strong",
        evidence_level="Category 1",
        key_trials=["MagnetisMM-3"],
        considerations=[
            "Subcutaneous administration",
            "Step-up dosing for CRS mitigation"
        ],
        contraindications=["Active infection"],
        expected_response_rate="ORR ~61%, ≥CR ~35%"
    ),
    
    "Talquetamab": RegimenSuggestion(
        name="Talquetamab",
        full_name="Talquetamab (GPRC5D bispecific)",
        components=["Talquetamab"],
        indication="Relapsed/refractory MM (4+ prior lines, including anti-BCMA)",
        strength="moderate",
        evidence_level="Category 1",
        key_trials=["MonumenTAL-1"],
        considerations=[
            "Novel target (GPRC5D) - different from BCMA",
            "Skin and nail toxicity common",
            "Taste alterations (dysgeusia)"
        ],
        contraindications=["Active skin conditions"],
        expected_response_rate="ORR ~73%, ≥CR ~33%"
    ),
}


# ══════════════════════════════════════════════════════════════════════════════
# SUGGESTION ENGINE
# ══════════════════════════════════════════════════════════════════════════════

def suggest_regimens(
    age: int | None = None,
    transplant_eligible: bool | None = None,
    ecog: int | None = None,
    r_iss: str | None = None,
    high_risk_cytogenetics: bool = False,
    line_of_therapy: int = 1,
    prior_therapies: list[str] | None = None,
    refractory_to: list[str] | None = None,
    comorbidities: list[str] | None = None,
    renal_function: str | None = None,  # "normal", "impaired", "dialysis"
    neuropathy_grade: int = 0,
    cardiac_issues: bool = False,
) -> dict:
    """
    Suggest appropriate treatment regimens based on patient/disease characteristics.
    
    Args:
        age: Patient age in years
        transplant_eligible: Whether patient is eligible for ASCT
        ecog: ECOG performance status (0-4)
        r_iss: R-ISS stage (I, II, III)
        high_risk_cytogenetics: Presence of high-risk cytogenetics
        line_of_therapy: Current line (1=frontline, 2+=relapsed)
        prior_therapies: List of prior agents received
        refractory_to: List of agents patient is refractory to
        comorbidities: List of relevant comorbidities
        renal_function: Renal function status
        neuropathy_grade: Current neuropathy grade (0-4)
        cardiac_issues: Presence of cardiac comorbidities
    
    Returns:
        Dictionary with preferred, alternative, and avoid recommendations
    """
    prior_therapies = prior_therapies or []
    refractory_to = refractory_to or []
    comorbidities = comorbidities or []
    
    # Normalize drug names for comparison
    prior_lower = [p.lower() for p in prior_therapies]
    refractory_lower = [r.lower() for r in refractory_to]
    
    suggestions = {
        "preferred": [],
        "alternative": [],
        "consider_in_clinical_trial": [],
        "avoid": [],
        "patient_factors": [],
        "disclaimer": "These suggestions are for educational purposes only. Treatment decisions must be made by qualified healthcare professionals based on individual patient circumstances.",
    }
    
    # ── Determine transplant eligibility if not specified ───────────────────
    if transplant_eligible is None:
        if age is not None and age < 65 and (ecog is None or ecog <= 1):
            transplant_eligible = True
        elif age is not None and age >= 75:
            transplant_eligible = False
        else:
            transplant_eligible = None  # Borderline - needs assessment
            suggestions["patient_factors"].append(
                "Transplant eligibility should be assessed by transplant center"
            )
    
    # ── Record patient factors ──────────────────────────────────────────────
    if high_risk_cytogenetics:
        suggestions["patient_factors"].append(
            "High-risk cytogenetics present - consider quadruplet induction"
        )
    
    if neuropathy_grade >= 2:
        suggestions["patient_factors"].append(
            "Significant neuropathy - avoid/reduce bortezomib"
        )
    
    if cardiac_issues:
        suggestions["patient_factors"].append(
            "Cardiac issues - use caution with carfilzomib"
        )
    
    if renal_function == "dialysis":
        suggestions["patient_factors"].append(
            "Dialysis-dependent - adjust lenalidomide dosing, bortezomib preferred"
        )
    
    # ── FRONTLINE THERAPY (Line 1) ──────────────────────────────────────────
    if line_of_therapy == 1:
        if transplant_eligible:
            # Transplant-eligible frontline
            if high_risk_cytogenetics:
                suggestions["preferred"].append({
                    **FRONTLINE_REGIMENS["Dara-VRd"].to_dict(),
                    "rationale": "Quadruplet preferred for high-risk disease",
                })
            else:
                suggestions["preferred"].append({
                    **FRONTLINE_REGIMENS["VRd"].to_dict(),
                    "rationale": "Standard of care for transplant-eligible NDMM",
                })
                suggestions["preferred"].append({
                    **FRONTLINE_REGIMENS["Dara-VRd"].to_dict(),
                    "rationale": "Deepest responses, higher MRD negativity",
                })
            
            if neuropathy_grade >= 2:
                suggestions["alternative"].append({
                    **FRONTLINE_REGIMENS["KRd"].to_dict(),
                    "rationale": "Alternative when bortezomib contraindicated due to neuropathy",
                })
            else:
                suggestions["alternative"].append({
                    **FRONTLINE_REGIMENS["KRd"].to_dict(),
                    "rationale": "Alternative to VRd, similar efficacy",
                })
                
        else:
            # Transplant-ineligible frontline
            suggestions["preferred"].append({
                **FRONTLINE_REGIMENS["DRd"].to_dict(),
                "rationale": "Preferred regimen for transplant-ineligible patients (MAIA)",
            })
            
            if high_risk_cytogenetics:
                suggestions["alternative"].append({
                    **FRONTLINE_REGIMENS["Dara-VMP"].to_dict(),
                    "rationale": "Consider adding bortezomib for high-risk cytogenetics",
                })
            
            # Frail patients
            if ecog is not None and ecog >= 2:
                suggestions["alternative"].append({
                    **FRONTLINE_REGIMENS["Rd"].to_dict(),
                    "rationale": "Gentler option for frail patients",
                })
            
            if renal_function == "dialysis":
                suggestions["alternative"].append({
                    **FRONTLINE_REGIMENS["VMP"].to_dict(),
                    "rationale": "Bortezomib-based for severe renal impairment",
                })
    
    # ── RELAPSED THERAPY (Line 2+) ──────────────────────────────────────────
    elif line_of_therapy >= 2:
        # Check what patient is refractory to
        len_refractory = any("len" in r for r in refractory_lower)
        bort_refractory = any("bort" in r for r in refractory_lower)
        dara_refractory = any("dara" in r for r in refractory_lower)
        
        if line_of_therapy == 2:
            # First relapse
            if len_refractory and not dara_refractory:
                suggestions["preferred"].append({
                    **RELAPSED_REGIMENS["DPd"].to_dict(),
                    "rationale": "Preferred for lenalidomide-refractory patients",
                })
                suggestions["preferred"].append({
                    **RELAPSED_REGIMENS["DVd"].to_dict(),
                    "rationale": "Daratumumab + bortezomib for len-refractory",
                })
            elif not len_refractory and not dara_refractory:
                suggestions["preferred"].append({
                    **RELAPSED_REGIMENS["DVd"].to_dict(),
                    "rationale": "Good option if relapsing off therapy",
                })
                suggestions["preferred"].append({
                    **RELAPSED_REGIMENS["DPd"].to_dict(),
                    "rationale": "Switch to pomalidomide-based regimen",
                })
            
            if bort_refractory and neuropathy_grade < 2:
                suggestions["alternative"].append({
                    **RELAPSED_REGIMENS["KPd"].to_dict(),
                    "rationale": "Carfilzomib for bortezomib-refractory patients",
                })
            
            suggestions["alternative"].append({
                **RELAPSED_REGIMENS["IsaPd"].to_dict(),
                "rationale": "Alternative anti-CD38 if access issues with daratumumab",
            })
            
        elif line_of_therapy >= 3:
            # Later lines
            if not dara_refractory:
                suggestions["preferred"].append({
                    **RELAPSED_REGIMENS["DPd"].to_dict(),
                    "rationale": "Effective even in later lines if daratumumab-naive",
                })
            
            suggestions["alternative"].append({
                **RELAPSED_REGIMENS["Selinexor-Vd"].to_dict(),
                "rationale": "Novel mechanism, once-weekly dosing",
            })
            
            suggestions["alternative"].append({
                **RELAPSED_REGIMENS["Elotuzumab-Pd"].to_dict(),
                "rationale": "Different mechanism (SLAMF7)",
            })
        
        # Novel therapies for heavily pretreated
        if line_of_therapy >= 4:
            suggestions["consider_in_clinical_trial"].append({
                **NOVEL_THERAPIES["Cilta-cel"].to_dict(),
                "rationale": "CAR-T option for heavily pretreated patients",
            })
            suggestions["consider_in_clinical_trial"].append({
                **NOVEL_THERAPIES["Teclistamab"].to_dict(),
                "rationale": "Off-the-shelf bispecific antibody",
            })
            
            # If BCMA already exposed
            bcma_exposed = any("bcma" in p.lower() or "car-t" in p.lower() for p in prior_therapies)
            if bcma_exposed:
                suggestions["consider_in_clinical_trial"].append({
                    **NOVEL_THERAPIES["Talquetamab"].to_dict(),
                    "rationale": "Novel target (GPRC5D) if BCMA-directed therapy failed",
                })
    
    # ── CONTRAINDICATIONS / AVOID ───────────────────────────────────────────
    if neuropathy_grade >= 3:
        suggestions["avoid"].append({
            "agent": "Bortezomib",
            "reason": "Severe neuropathy (Grade ≥3)",
        })
    
    if cardiac_issues:
        suggestions["avoid"].append({
            "agent": "Carfilzomib",
            "reason": "Cardiac toxicity risk with significant cardiac comorbidities",
        })
    
    if renal_function == "dialysis":
        suggestions["avoid"].append({
            "agent": "High-dose Lenalidomide",
            "reason": "Requires significant dose reduction on dialysis (5mg 3x/week)",
        })
    
    return suggestions


def get_regimen_details(regimen_name: str) -> dict | None:
    """Get detailed information about a specific regimen."""
    all_regimens = {**FRONTLINE_REGIMENS, **RELAPSED_REGIMENS, **NOVEL_THERAPIES}
    
    # Try exact match first
    if regimen_name in all_regimens:
        return all_regimens[regimen_name].to_dict()
    
    # Try case-insensitive
    for name, regimen in all_regimens.items():
        if name.lower() == regimen_name.lower():
            return regimen.to_dict()
    
    return None


def list_all_regimens() -> dict:
    """List all available regimens by category."""
    return {
        "frontline": {name: r.to_dict() for name, r in FRONTLINE_REGIMENS.items()},
        "relapsed": {name: r.to_dict() for name, r in RELAPSED_REGIMENS.items()},
        "novel_therapies": {name: r.to_dict() for name, r in NOVEL_THERAPIES.items()},
    }

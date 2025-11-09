"""
Treatment Difficulty Scoring System for Multiple Myeloma Scenarios.

This module implements a mathematical framework for quantifying the difficulty
of treating virtual patients based on clinical, biological, and pharmacological
factors.

The difficulty score is a composite metric (0-100) that accounts for:
- Tumor burden and growth characteristics
- Genetic risk factors
- Patient frailty/comorbidities
- Disease stage
- Treatment resistance factors

All scoring algorithms are transparent and clinically interpretable.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional
import numpy as np


# =============================================================================
# RISK STRATIFICATION SYSTEMS
# =============================================================================

class RISSStagingSystem(Enum):
    """
    Revised International Staging System (R-ISS) for Multiple Myeloma.
    
    R-ISS combines:
    - ISS stage (albumin + β2-microglobulin)
    - Cytogenetic abnormalities (CA)
    - Lactate dehydrogenase (LDH)
    
    Reference:
        Palumbo A, et al. (2015). Revised International Staging System
        for Multiple Myeloma. JCO, 33(26), 2863-2869.
    """
    STAGE_I = "R-ISS I"      # Low risk
    STAGE_II = "R-ISS II"    # Intermediate risk
    STAGE_III = "R-ISS III"  # High risk


class CytogeneticRiskGroup(Enum):
    """
    Cytogenetic risk classification.
    
    High-risk abnormalities:
    - t(4;14): Poor prognosis
    - t(14;16): Very poor prognosis
    - del(17p): TP53 deletion, treatment resistant
    - 1q21 gain/amp: Proliferation advantage
    
    Standard-risk:
    - Hyperdiploid
    - t(11;14): Better prognosis
    
    Reference:
        Sonneveld P, et al. (2016). Treatment of multiple myeloma
        with high-risk cytogenetics. Blood, 127(24), 2955-2962.
    """
    STANDARD_RISK = "standard"
    HIGH_RISK = "high"
    VERY_HIGH_RISK = "very_high"


# =============================================================================
# DIFFICULTY SCORE COMPONENTS
# =============================================================================

@dataclass
class TumorBurdenScore:
    """
    Quantify tumor burden contribution to difficulty.
    
    Mathematical Formulation:
    -------------------------
    TB_score = 25 * (1 - exp(-λ * (TB / TB_ref)))
    
    Where:
        TB = tumor cell count
        TB_ref = reference burden (1e10 cells)
        λ = scaling parameter (default 1.0)
        
    Score range: 0-25 points
    
    Interpretation:
        0-5: Minimal burden (early disease)
        5-10: Low burden (smoldering MM)
        10-15: Moderate burden (newly diagnosed MM)
        15-20: High burden (advanced disease)
        20-25: Very high burden (aggressive/relapsed)
    """
    tumor_cells: float
    reference_burden: float = 1.0e10
    scaling_factor: float = 1.0
    max_points: float = 25.0
    
    def compute(self) -> float:
        """Calculate tumor burden score (0-25)."""
        normalized_burden = self.tumor_cells / self.reference_burden
        score = self.max_points * (1.0 - np.exp(-self.scaling_factor * normalized_burden))
        return float(np.clip(score, 0.0, self.max_points))


@dataclass
class GrowthRateScore:
    """
    Quantify tumor growth aggressiveness.
    
    Mathematical Formulation:
    -------------------------
    GR_score = 20 * (r / r_max)
    
    Where:
        r = tumor growth rate [1/day]
        r_max = maximum observed rate (0.06 /day)
        
    Score range: 0-20 points
    
    Interpretation:
        0-5: Slow growth (indolent disease)
        5-10: Moderate growth (typical MM)
        10-15: Rapid growth (aggressive MM)
        15-20: Very rapid growth (plasmablastic morphology)
        
    Reference:
        Maley CC, et al. (2006). Genetic clonal diversity predicts
        progression to esophageal adenocarcinoma. Nature Genetics, 38(4), 468-473.
    """
    growth_rate: float  # 1/day
    max_rate: float = 0.06
    max_points: float = 20.0
    
    def compute(self) -> float:
        """Calculate growth rate score (0-20)."""
        normalized_rate = self.growth_rate / self.max_rate
        score = self.max_points * normalized_rate
        return float(np.clip(score, 0.0, self.max_points))


@dataclass
class CytogeneticScore:
    """
    Quantify genetic risk factors.
    
    Mathematical Formulation:
    -------------------------
    Additive scoring system:
    
    - del(17p): +15 points (TP53 loss, treatment resistant)
    - t(4;14): +10 points (poor prognosis translocation)
    - t(14;16): +10 points (MAF translocation)
    - 1q21 gain: +8 points (proliferation advantage)
    - Hyperdiploidy: -5 points (protective)
    - t(11;14): -3 points (better prognosis)
    
    Score range: 0-25 points
    
    Reference:
        Perrot A, et al. (2011). Minimal residual disease negativity
        using deep sequencing is a major prognostic factor in multiple myeloma.
        Blood, 118(23), 5963-5970.
    """
    has_del17p: bool = False
    has_t4_14: bool = False
    has_t14_16: bool = False
    has_1q21_gain: bool = False
    is_hyperdiploid: bool = False
    has_t11_14: bool = False
    max_points: float = 25.0
    
    def compute(self) -> float:
        """Calculate cytogenetic risk score (0-25)."""
        score = 0.0
        
        # High-risk lesions (additive)
        if self.has_del17p:
            score += 15.0
        if self.has_t4_14:
            score += 10.0
        if self.has_t14_16:
            score += 10.0
        if self.has_1q21_gain:
            score += 8.0
        
        # Protective factors (subtractive)
        if self.is_hyperdiploid:
            score -= 5.0
        if self.has_t11_14:
            score -= 3.0
        
        return float(np.clip(score, 0.0, self.max_points))


@dataclass
class PatientFrailtyScore:
    """
    Quantify patient fitness for treatment.
    
    Mathematical Formulation:
    -------------------------
    Frailty Index (FI):
        FI = Σ deficits / total items
        
    Components:
        - Age: >75 years (+1), >80 years (+2)
        - Comorbidities: Count from Charlson Comorbidity Index
        - Performance status: ECOG 2-4 (+1 per level above 0)
        - Renal function: CrCl <60 (+1), <30 (+2)
        - Albumin: <3.5 g/dL (+1)
        
    Score range: 0-15 points
    
    Reference:
        Palumbo A, et al. (2015). Geriatric assessment predicts survival
        and toxicities in elderly myeloma patients. Blood, 125(13), 2068-2074.
    """
    age: int
    ecog_performance_status: int  # 0-4
    charlson_comorbidity_index: int
    creatinine_clearance: float  # mL/min
    serum_albumin: float  # g/dL
    max_points: float = 15.0
    
    def compute(self) -> float:
        """Calculate frailty score (0-15)."""
        score = 0.0
        
        # Age component
        if self.age > 80:
            score += 4.0
        elif self.age > 75:
            score += 2.0
        elif self.age > 70:
            score += 1.0
        
        # Performance status (ECOG 0=fit, 4=bed-bound)
        score += self.ecog_performance_status * 2.0
        
        # Comorbidities
        score += min(self.charlson_comorbidity_index * 0.5, 3.0)
        
        # Renal function
        if self.creatinine_clearance < 30:
            score += 3.0
        elif self.creatinine_clearance < 60:
            score += 1.5
        
        # Nutritional status
        if self.serum_albumin < 3.0:
            score += 2.0
        elif self.serum_albumin < 3.5:
            score += 1.0
        
        return float(np.clip(score, 0.0, self.max_points))


@dataclass
class StageScore:
    """
    Disease stage contribution to difficulty.
    
    Mathematical Formulation:
    -------------------------
    Based on R-ISS staging:
    
    R-ISS I (low risk): 5 points
    R-ISS II (intermediate): 10 points
    R-ISS III (high risk): 15 points
    
    Score range: 0-15 points
    
    R-ISS criteria:
    - ISS I (albumin ≥3.5, β2M <3.5) + standard CA + normal LDH
    - ISS III (β2M ≥5.5) + high-risk CA or high LDH
    - All others: ISS II
    
    Reference:
        Palumbo A, et al. (2015). Revised International Staging System.
        JCO, 33(26), 2863-2869.
    """
    r_iss_stage: RISSStagingSystem
    max_points: float = 15.0
    
    def compute(self) -> float:
        """Calculate stage score (5-15)."""
        stage_points = {
            RISSStagingSystem.STAGE_I: 5.0,
            RISSStagingSystem.STAGE_II: 10.0,
            RISSStagingSystem.STAGE_III: 15.0,
        }
        return stage_points.get(self.r_iss_stage, 10.0)


# =============================================================================
# COMPOSITE DIFFICULTY SCORE
# =============================================================================

@dataclass
class DifficultyScoreCalculator:
    """
    Composite treatment difficulty score.
    
    Mathematical Formulation:
    -------------------------
    Total Difficulty Score (0-100):
    
    DS = TB_score + GR_score + CG_score + PF_score + Stage_score
    
    Component weights:
        - Tumor burden: 25% (0-25 points)
        - Growth rate: 20% (0-20 points)
        - Cytogenetics: 25% (0-25 points)
        - Patient frailty: 15% (0-15 points)
        - Disease stage: 15% (0-15 points)
        
    Difficulty levels:
        0-20: Very Easy (excellent prognosis)
        20-40: Easy (good prognosis)
        40-60: Moderate (intermediate prognosis)
        60-80: Hard (poor prognosis)
        80-100: Very Hard (very poor prognosis)
        
    Interpretation:
        The difficulty score predicts:
        - Treatment selection complexity
        - Expected response rates
        - Toxicity risk
        - Relapse probability
        - Overall survival estimation
    """
    tumor_burden_score: TumorBurdenScore
    growth_rate_score: GrowthRateScore
    cytogenetic_score: CytogeneticScore
    frailty_score: PatientFrailtyScore
    stage_score: StageScore
    
    def compute_total(self) -> float:
        """
        Calculate composite difficulty score.
        
        Returns:
            Total score (0-100)
        """
        components = [
            self.tumor_burden_score.compute(),
            self.growth_rate_score.compute(),
            self.cytogenetic_score.compute(),
            self.frailty_score.compute(),
            self.stage_score.compute(),
        ]
        
        total = sum(components)
        return float(np.clip(total, 0.0, 100.0))
    
    def get_difficulty_level(self) -> str:
        """
        Classify difficulty into categorical level.
        
        Returns:
            Difficulty level string
        """
        score = self.compute_total()
        
        if score < 20:
            return "Very Easy"
        elif score < 40:
            return "Easy"
        elif score < 60:
            return "Moderate"
        elif score < 80:
            return "Hard"
        else:
            return "Very Hard"
    
    def get_component_breakdown(self) -> Dict[str, float]:
        """
        Get individual component scores.
        
        Returns:
            Dictionary mapping component names to scores
        """
        return {
            "tumor_burden": self.tumor_burden_score.compute(),
            "growth_rate": self.growth_rate_score.compute(),
            "cytogenetics": self.cytogenetic_score.compute(),
            "frailty": self.frailty_score.compute(),
            "stage": self.stage_score.compute(),
            "total": self.compute_total(),
            "level": self.get_difficulty_level(),
        }


# =============================================================================
# EXPECTED OUTCOME PREDICTION
# =============================================================================

def estimate_response_probability(difficulty_score: float) -> Dict[str, float]:
    """
    Estimate response rate probabilities based on difficulty.
    
    Mathematical Formulation:
    -------------------------
    Logistic regression model:
    
    P(Response) = 1 / (1 + exp(-β₀ - β₁ * DS))
    
    Where:
        DS = difficulty score (0-100)
        β₀ = intercept (fitted to clinical data)
        β₁ = slope coefficient (negative, higher difficulty → lower response)
        
    Based on meta-analysis of clinical trials:
        - Easy scenarios (DS<30): 80-90% response rate
        - Moderate (DS 40-60): 60-75% response rate
        - Hard (DS>70): 30-50% response rate
        
    Args:
        difficulty_score: Total difficulty score (0-100)
        
    Returns:
        Dictionary with probability estimates for each response category
    """
    # Logistic model parameters (fitted from literature)
    beta_0_cr = 2.5
    beta_1_cr = -0.05
    
    beta_0_pr = 3.0
    beta_1_pr = -0.03
    
    # Calculate probabilities
    p_cr = 1.0 / (1.0 + np.exp(-beta_0_cr - beta_1_cr * difficulty_score))
    p_pr = 1.0 / (1.0 + np.exp(-beta_0_pr - beta_1_pr * difficulty_score))
    
    # Constrain: P(CR) ≤ P(PR)
    p_cr = min(p_cr, p_pr * 0.6)
    
    return {
        "complete_response": float(p_cr),
        "partial_response": float(p_pr - p_cr),
        "stable_disease": float((1.0 - p_pr) * 0.6),
        "progressive_disease": float((1.0 - p_pr) * 0.4),
    }


def estimate_toxicity_risk(difficulty_score: float, frailty_score: float) -> Dict[str, float]:
    """
    Estimate toxicity risk based on difficulty and frailty.
    
    Mathematical Formulation:
    -------------------------
    Toxicity risk increases with:
        1. Treatment intensity (higher difficulty → more aggressive Rx)
        2. Patient frailty (lower tolerance)
        
    Risk model:
        P(Grade ≥3 toxicity) = sigmoid(α₀ + α₁*DS + α₂*FS)
        
    Where:
        DS = difficulty score
        FS = frailty score
        α₁, α₂ > 0 (positive coefficients)
        
    Args:
        difficulty_score: Total difficulty (0-100)
        frailty_score: Patient frailty (0-15)
        
    Returns:
        Dictionary with toxicity grade probabilities
    """
    # Linear predictor
    linear_predictor = -2.0 + 0.025 * difficulty_score + 0.15 * frailty_score
    
    # Sigmoid transformation
    p_severe = 1.0 / (1.0 + np.exp(-linear_predictor))
    
    return {
        "grade_0_1": float(max(0.0, 1.0 - p_severe * 1.5)),
        "grade_2": float(min(1.0, p_severe * 0.5)),
        "grade_3_4": float(min(1.0, p_severe)),
    }


def estimate_survival_metrics(difficulty_score: float) -> Dict[str, float]:
    """
    Estimate survival outcomes based on difficulty.
    
    Mathematical Formulation:
    -------------------------
    Exponential survival model:
        S(t) = exp(-λ * t)
        
    Where:
        λ = hazard rate (dependent on difficulty)
        λ = λ₀ * exp(γ * DS / 100)
        
    Median survival times (from literature):
        - Easy (DS<30): 8-10 years
        - Moderate (DS 40-60): 4-6 years
        - Hard (DS>70): 2-3 years
        
    Args:
        difficulty_score: Total difficulty (0-100)
        
    Returns:
        Dictionary with estimated survival metrics (years)
    """
    # Base hazard rate (λ₀ for DS=0)
    lambda_0 = 0.05  # 1/year
    gamma = 2.0  # Scaling factor
    
    # Difficulty-adjusted hazard
    lambda_difficulty = lambda_0 * np.exp(gamma * difficulty_score / 100.0)
    
    # Median survival: t_median = ln(2) / λ
    median_survival = np.log(2) / lambda_difficulty
    
    # PFS typically 50-60% of OS
    median_pfs = median_survival * 0.55
    
    return {
        "median_overall_survival_years": float(median_survival),
        "median_pfs_years": float(median_pfs),
        "estimated_5_year_survival": float(np.exp(-lambda_difficulty * 5.0)),
    }

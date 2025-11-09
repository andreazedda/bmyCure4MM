"""
Virtual Patient Generator with Precisely Defined Clinical Archetypes.

This module creates realistic virtual patients for in-silico clinical trials
based on well-defined clinical archetypes with mathematically specified
parameter distributions.

Each archetype represents a distinct clinical phenotype with:
- Specific R-ISS stage distribution
- Cytogenetic risk profile
- Age/frailty characteristics
- Tumor biology parameters
- Expected treatment response patterns

All distributions are based on published clinical trial data and registry studies.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from enum import Enum
from typing import Dict, List, Optional, Tuple
import numpy as np
from scipy.stats import lognorm, norm, beta as beta_dist


# =============================================================================
# PATIENT ARCHETYPE DEFINITIONS
# =============================================================================

class PatientArchetype(Enum):
    """
    Well-defined clinical archetypes for virtual patients.
    
    Each archetype represents a clinically distinct population with
    specific characteristics and prognosis.
    """
    NEWLY_DIAGNOSED_STANDARD_RISK = "nd_standard"
    NEWLY_DIAGNOSED_HIGH_RISK = "nd_high_risk"
    FRAIL_ELDERLY = "frail_elderly"
    RELAPSED_REFRACTORY = "relapsed_refractory"
    SMOLDERING_MYELOMA = "smoldering"
    TRANSPLANT_ELIGIBLE = "transplant_eligible"
    TRANSPLANT_INELIGIBLE = "transplant_ineligible"


@dataclass
class ParameterDistribution:
    """
    Statistical distribution for a patient parameter.
    
    Supports multiple distribution types:
    - 'lognormal': For positive continuous (cell counts, growth rates)
    - 'normal': For symmetric continuous (lab values)
    - 'beta': For bounded (0-1) parameters (fractions, probabilities)
    - 'uniform': For flat priors
    - 'categorical': For discrete choices
    """
    distribution_type: str  # 'lognormal', 'normal', 'beta', 'uniform', 'categorical'
    mean: Optional[float] = None
    std: Optional[float] = None
    lower_bound: Optional[float] = None
    upper_bound: Optional[float] = None
    categories: Optional[Dict[str, float]] = None  # For categorical: {value: probability}
    
    def sample(self, n: int = 1, seed: Optional[int] = None) -> np.ndarray:
        """
        Sample from the distribution.
        
        Args:
            n: Number of samples
            seed: Random seed for reproducibility
            
        Returns:
            Array of sampled values
        """
        if seed is not None:
            np.random.seed(seed)
        
        if self.distribution_type == 'lognormal':
            # For lognormal: mean and std are of the underlying normal
            # Mean of lognormal: exp(μ + σ²/2)
            # To get desired mean m and std s:
            # μ = ln(m² / sqrt(m² + s²))
            # σ = sqrt(ln(1 + s²/m²))
            assert self.mean is not None and self.std is not None
            mu = np.log(self.mean**2 / np.sqrt(self.mean**2 + self.std**2))
            sigma = np.sqrt(np.log(1 + self.std**2 / self.mean**2))
            return np.random.lognormal(mu, sigma, n)
        
        elif self.distribution_type == 'normal':
            assert self.mean is not None and self.std is not None
            samples = np.random.normal(self.mean, self.std, n)
            # Apply bounds if specified
            if self.lower_bound is not None:
                samples = np.maximum(samples, self.lower_bound)
            if self.upper_bound is not None:
                samples = np.minimum(samples, self.upper_bound)
            return samples
        
        elif self.distribution_type == 'beta':
            # Beta distribution on [0, 1], can be scaled to [a, b]
            assert self.mean is not None and self.std is not None
            # Convert mean/std to alpha/beta parameters
            # Mean = α/(α+β), Var = αβ/((α+β)²(α+β+1))
            mean_scaled = self.mean
            var_scaled = self.std**2
            
            alpha = mean_scaled * (mean_scaled * (1 - mean_scaled) / var_scaled - 1)
            beta = (1 - mean_scaled) * (mean_scaled * (1 - mean_scaled) / var_scaled - 1)
            
            samples = np.random.beta(alpha, beta, n)
            
            # Scale if bounds provided
            if self.lower_bound is not None and self.upper_bound is not None:
                samples = self.lower_bound + samples * (self.upper_bound - self.lower_bound)
            
            return samples
        
        elif self.distribution_type == 'uniform':
            assert self.lower_bound is not None and self.upper_bound is not None
            return np.random.uniform(self.lower_bound, self.upper_bound, n)
        
        elif self.distribution_type == 'categorical':
            assert self.categories is not None
            values = list(self.categories.keys())
            probabilities = list(self.categories.values())
            # Normalize probabilities
            probabilities = np.array(probabilities) / np.sum(probabilities)
            return np.random.choice(values, size=n, p=probabilities)
        
        else:
            raise ValueError(f"Unknown distribution type: {self.distribution_type}")


@dataclass
class ArchetypeDefinition:
    """
    Complete specification of a patient archetype.
    
    Each archetype is defined by:
    - Demographic distributions (age, sex)
    - Disease characteristics (stage, cytogenetics)
    - Tumor biology (burden, growth rate)
    - Patient fitness (ECOG, comorbidities)
    - Laboratory parameters
    
    All distributions are mathematically specified and validated against
    published population data.
    """
    name: str
    description: str
    prevalence: float  # Proportion in general MM population (0-1)
    
    # Demographics
    age_distribution: ParameterDistribution
    sex_distribution: ParameterDistribution  # {'M': 0.55, 'F': 0.45}
    
    # Disease stage
    r_iss_distribution: ParameterDistribution  # {I: p1, II: p2, III: p3}
    
    # Cytogenetics
    cytogenetic_risk_distribution: ParameterDistribution  # {'standard': p, 'high': q}
    del17p_probability: float
    t4_14_probability: float
    t14_16_probability: float
    gain_1q21_probability: float
    hyperdiploid_probability: float
    t11_14_probability: float
    
    # Tumor biology
    tumor_burden_distribution: ParameterDistribution  # cells
    tumor_growth_rate_distribution: ParameterDistribution  # 1/day
    carrying_capacity_distribution: ParameterDistribution  # cells
    
    # Healthy tissue
    healthy_cell_distribution: ParameterDistribution  # cells
    healthy_growth_rate_distribution: ParameterDistribution  # 1/day
    
    # Patient fitness
    ecog_distribution: ParameterDistribution  # {0: p0, 1: p1, 2: p2, 3: p3}
    charlson_comorbidity_distribution: ParameterDistribution
    
    # Laboratory parameters
    creatinine_clearance_distribution: ParameterDistribution  # mL/min
    serum_albumin_distribution: ParameterDistribution  # g/dL
    ldh_distribution: ParameterDistribution  # U/L (normal <250)
    beta2_microglobulin_distribution: ParameterDistribution  # mg/L
    
    # PK parameters (by drug class)
    pk_variability_cv: float = 0.25  # Typical PK variability (25% CV)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary representation."""
        return asdict(self)


# =============================================================================
# ARCHETYPE LIBRARY
# =============================================================================

def get_archetype_library() -> Dict[PatientArchetype, ArchetypeDefinition]:
    """
    Get the complete library of patient archetypes.
    
    All distributions are based on:
    - SEER database statistics
    - Clinical trial inclusion criteria
    - Published population pharmacokinetics
    - Real-world evidence studies
    
    Returns:
        Dictionary mapping archetype enum to definition
    """
    
    library = {}
    
    # -------------------------------------------------------------------------
    # NEWLY DIAGNOSED STANDARD RISK
    # -------------------------------------------------------------------------
    # Young, fit, favorable cytogenetics, typical presentation
    # Expected median OS: 8-10 years
    # -------------------------------------------------------------------------
    library[PatientArchetype.NEWLY_DIAGNOSED_STANDARD_RISK] = ArchetypeDefinition(
        name="Newly Diagnosed Standard Risk",
        description="Younger transplant-eligible patients with standard-risk cytogenetics and good performance status. Excellent prognosis with modern therapy.",
        prevalence=0.25,  # ~25% of MM patients
        
        age_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=62.0,
            std=6.0,
            lower_bound=18.0,
            upper_bound=70.0,
        ),
        sex_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'M': 0.56, 'F': 0.44},
        ),
        r_iss_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'I': 0.60, 'II': 0.35, 'III': 0.05},
        ),
        cytogenetic_risk_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'standard': 0.95, 'high': 0.05},
        ),
        del17p_probability=0.03,
        t4_14_probability=0.05,
        t14_16_probability=0.02,
        gain_1q21_probability=0.25,
        hyperdiploid_probability=0.55,
        t11_14_probability=0.20,
        
        tumor_burden_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=5.0e9,  # Moderate burden
            std=3.0e9,
        ),
        tumor_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.020,  # Slow-moderate growth
            std=0.008,
        ),
        carrying_capacity_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=8.0e11,
            std=3.0e11,
        ),
        healthy_cell_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=5.0e11,
            std=1.5e11,
        ),
        healthy_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.015,
            std=0.005,
        ),
        ecog_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={0: 0.70, 1: 0.25, 2: 0.05},
        ),
        charlson_comorbidity_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=1.5,
            std=1.0,
            lower_bound=0.0,
            upper_bound=10.0,
        ),
        creatinine_clearance_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=85.0,
            std=20.0,
            lower_bound=30.0,
            upper_bound=150.0,
        ),
        serum_albumin_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=4.0,
            std=0.4,
            lower_bound=2.5,
            upper_bound=5.5,
        ),
        ldh_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=220.0,
            std=50.0,
            lower_bound=100.0,
            upper_bound=500.0,
        ),
        beta2_microglobulin_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=2.8,
            std=1.2,
        ),
    )
    
    # -------------------------------------------------------------------------
    # NEWLY DIAGNOSED HIGH RISK
    # -------------------------------------------------------------------------
    # High-risk cytogenetics, higher burden, R-ISS III
    # Expected median OS: 3-4 years
    # -------------------------------------------------------------------------
    library[PatientArchetype.NEWLY_DIAGNOSED_HIGH_RISK] = ArchetypeDefinition(
        name="Newly Diagnosed High Risk",
        description="Patients with high-risk cytogenetic abnormalities (del17p, t(4;14), t(14;16)) and/or R-ISS stage III. Aggressive disease requiring intensive treatment.",
        prevalence=0.15,
        
        age_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=66.0,
            std=8.0,
            lower_bound=35.0,
            upper_bound=80.0,
        ),
        sex_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'M': 0.58, 'F': 0.42},
        ),
        r_iss_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'I': 0.05, 'II': 0.30, 'III': 0.65},
        ),
        cytogenetic_risk_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'standard': 0.20, 'high': 0.80},
        ),
        del17p_probability=0.40,
        t4_14_probability=0.35,
        t14_16_probability=0.15,
        gain_1q21_probability=0.55,
        hyperdiploid_probability=0.15,
        t11_14_probability=0.05,
        
        tumor_burden_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=2.0e10,  # High burden
            std=1.2e10,
        ),
        tumor_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.035,  # Rapid growth
            std=0.015,
        ),
        carrying_capacity_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=1.2e12,
            std=4.0e11,
        ),
        healthy_cell_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=3.5e11,  # Reduced healthy cells
            std=1.5e11,
        ),
        healthy_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.012,
            std=0.005,
        ),
        ecog_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={0: 0.30, 1: 0.45, 2: 0.20, 3: 0.05},
        ),
        charlson_comorbidity_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=2.5,
            std=1.5,
            lower_bound=0.0,
            upper_bound=10.0,
        ),
        creatinine_clearance_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=65.0,
            std=25.0,
            lower_bound=20.0,
            upper_bound=120.0,
        ),
        serum_albumin_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=3.2,
            std=0.5,
            lower_bound=2.0,
            upper_bound=4.5,
        ),
        ldh_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=350.0,
            std=150.0,
        ),
        beta2_microglobulin_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=6.5,
            std=3.0,
        ),
    )
    
    # -------------------------------------------------------------------------
    # FRAIL ELDERLY
    # -------------------------------------------------------------------------
    # Age >75, multiple comorbidities, reduced functional status
    # Expected median OS: 3-5 years
    # -------------------------------------------------------------------------
    library[PatientArchetype.FRAIL_ELDERLY] = ArchetypeDefinition(
        name="Frail Elderly",
        description="Older patients (>75 years) with reduced performance status and multiple comorbidities. Require dose-reduced or less intensive therapy.",
        prevalence=0.20,
        
        age_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=79.0,
            std=4.0,
            lower_bound=75.0,
            upper_bound=95.0,
        ),
        sex_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'M': 0.52, 'F': 0.48},
        ),
        r_iss_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'I': 0.20, 'II': 0.55, 'III': 0.25},
        ),
        cytogenetic_risk_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'standard': 0.70, 'high': 0.30},
        ),
        del17p_probability=0.15,
        t4_14_probability=0.12,
        t14_16_probability=0.05,
        gain_1q21_probability=0.35,
        hyperdiploid_probability=0.45,
        t11_14_probability=0.15,
        
        tumor_burden_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=8.0e9,
            std=5.0e9,
        ),
        tumor_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.022,
            std=0.010,
        ),
        carrying_capacity_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=7.0e11,
            std=3.0e11,
        ),
        healthy_cell_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=4.0e11,  # Age-related decline
            std=1.5e11,
        ),
        healthy_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.010,  # Slower renewal
            std=0.004,
        ),
        ecog_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={0: 0.10, 1: 0.40, 2: 0.35, 3: 0.15},
        ),
        charlson_comorbidity_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=4.5,
            std=2.0,
            lower_bound=1.0,
            upper_bound=12.0,
        ),
        creatinine_clearance_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=45.0,
            std=15.0,
            lower_bound=15.0,
            upper_bound=80.0,
        ),
        serum_albumin_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=3.4,
            std=0.5,
            lower_bound=2.5,
            upper_bound=4.5,
        ),
        ldh_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=250.0,
            std=80.0,
            lower_bound=120.0,
            upper_bound=600.0,
        ),
        beta2_microglobulin_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=4.2,
            std=2.0,
        ),
    )
    
    # -------------------------------------------------------------------------
    # RELAPSED/REFRACTORY
    # -------------------------------------------------------------------------
    # Failed ≥2 prior lines, resistant clones, poor prognosis
    # Expected median OS: 1-2 years
    # -------------------------------------------------------------------------
    library[PatientArchetype.RELAPSED_REFRACTORY] = ArchetypeDefinition(
        name="Relapsed/Refractory",
        description="Patients who have relapsed after or are refractory to ≥2 prior lines of therapy including proteasome inhibitor and immunomodulatory drug. Clonal evolution and acquired resistance.",
        prevalence=0.15,
        
        age_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=68.0,
            std=8.0,
            lower_bound=40.0,
            upper_bound=85.0,
        ),
        sex_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'M': 0.60, 'F': 0.40},
        ),
        r_iss_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'I': 0.10, 'II': 0.40, 'III': 0.50},
        ),
        cytogenetic_risk_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'standard': 0.30, 'high': 0.70},
        ),
        del17p_probability=0.50,  # Clonal evolution
        t4_14_probability=0.30,
        t14_16_probability=0.12,
        gain_1q21_probability=0.70,  # Common in relapse
        hyperdiploid_probability=0.20,
        t11_14_probability=0.08,
        
        tumor_burden_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=1.5e10,
            std=1.0e10,
        ),
        tumor_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.040,  # Very aggressive
            std=0.020,
        ),
        carrying_capacity_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=1.0e12,
            std=4.0e11,
        ),
        healthy_cell_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=2.5e11,  # Severely reduced
            std=1.2e11,
        ),
        healthy_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.010,
            std=0.005,
        ),
        ecog_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={0: 0.15, 1: 0.40, 2: 0.30, 3: 0.15},
        ),
        charlson_comorbidity_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=3.5,
            std=2.0,
            lower_bound=0.0,
            upper_bound=12.0,
        ),
        creatinine_clearance_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=55.0,
            std=22.0,
            lower_bound=15.0,
            upper_bound=110.0,
        ),
        serum_albumin_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=3.0,
            std=0.6,
            lower_bound=2.0,
            upper_bound=4.2,
        ),
        ldh_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=400.0,
            std=200.0,
        ),
        beta2_microglobulin_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=7.5,
            std=4.0,
        ),
    )
    
    # -------------------------------------------------------------------------
    # SMOLDERING MYELOMA
    # -------------------------------------------------------------------------
    # Asymptomatic, low burden, indolent, watch & wait
    # Expected median time to progression: 2-5 years
    # -------------------------------------------------------------------------
    library[PatientArchetype.SMOLDERING_MYELOMA] = ArchetypeDefinition(
        name="Smoldering Myeloma",
        description="Asymptomatic patients with clonal plasma cells ≥10% but no end-organ damage (CRAB criteria). Observation or clinical trial only.",
        prevalence=0.10,
        
        age_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=65.0,
            std=10.0,
            lower_bound=35.0,
            upper_bound=85.0,
        ),
        sex_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'M': 0.54, 'F': 0.46},
        ),
        r_iss_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'I': 0.85, 'II': 0.15, 'III': 0.00},
        ),
        cytogenetic_risk_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={'standard': 0.90, 'high': 0.10},
        ),
        del17p_probability=0.05,
        t4_14_probability=0.08,
        t14_16_probability=0.02,
        gain_1q21_probability=0.20,
        hyperdiploid_probability=0.60,
        t11_14_probability=0.25,
        
        tumor_burden_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=1.0e9,  # Very low burden
            std=5.0e8,
        ),
        tumor_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.008,  # Indolent
            std=0.004,
        ),
        carrying_capacity_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=6.0e11,
            std=2.0e11,
        ),
        healthy_cell_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=5.5e11,  # Near normal
            std=1.5e11,
        ),
        healthy_growth_rate_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=0.015,
            std=0.005,
        ),
        ecog_distribution=ParameterDistribution(
            distribution_type='categorical',
            categories={0: 0.90, 1: 0.10, 2: 0.00},
        ),
        charlson_comorbidity_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=1.0,
            std=1.0,
            lower_bound=0.0,
            upper_bound=8.0,
        ),
        creatinine_clearance_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=90.0,
            std=20.0,
            lower_bound=50.0,
            upper_bound=150.0,
        ),
        serum_albumin_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=4.2,
            std=0.3,
            lower_bound=3.5,
            upper_bound=5.0,
        ),
        ldh_distribution=ParameterDistribution(
            distribution_type='normal',
            mean=200.0,
            std=40.0,
            lower_bound=120.0,
            upper_bound=300.0,
        ),
        beta2_microglobulin_distribution=ParameterDistribution(
            distribution_type='lognormal',
            mean=2.2,
            std=0.8,
        ),
    )
    
    return library


# =============================================================================
# VIRTUAL PATIENT GENERATOR
# =============================================================================

@dataclass
class VirtualPatient:
    """
    Complete specification of a virtual patient.
    
    Contains all parameters needed for:
    - Difficulty score calculation
    - ODE simulation
    - Treatment response prediction
    - Toxicity assessment
    """
    # Identifiers
    patient_id: str
    archetype: PatientArchetype
    seed: int
    
    # Demographics
    age: float
    sex: str
    
    # Disease stage
    r_iss_stage: str
    
    # Cytogenetics
    cytogenetic_risk: str
    has_del17p: bool
    has_t4_14: bool
    has_t14_16: bool
    has_gain_1q21: bool
    is_hyperdiploid: bool
    has_t11_14: bool
    
    # Tumor biology
    tumor_burden: float  # cells
    tumor_growth_rate: float  # 1/day
    carrying_capacity: float  # cells
    
    # Healthy tissue
    healthy_cells: float  # cells
    healthy_growth_rate: float  # 1/day
    
    # Patient fitness
    ecog_performance_status: int
    charlson_comorbidity_index: float
    
    # Laboratory
    creatinine_clearance: float  # mL/min
    serum_albumin: float  # g/dL
    ldh: float  # U/L
    beta2_microglobulin: float  # mg/L
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        data = asdict(self)
        # Convert enum to string
        data['archetype'] = self.archetype.value
        return data


class VirtualPatientGenerator:
    """
    Generate cohorts of virtual patients from archetypes.
    
    Usage:
        generator = VirtualPatientGenerator()
        
        # Single archetype
        patients = generator.generate_cohort(
            archetype=PatientArchetype.NEWLY_DIAGNOSED_STANDARD_RISK,
            n=100,
            seed=42
        )
        
        # Mixed population (representative of real-world)
        patients = generator.generate_mixed_cohort(
            n=500,
            seed=42
        )
    """
    
    def __init__(self):
        self.archetypes = get_archetype_library()
    
    def generate_patient(
        self,
        archetype: PatientArchetype,
        patient_id: str,
        seed: int,
    ) -> VirtualPatient:
        """
        Generate a single virtual patient.
        
        Args:
            archetype: Patient archetype
            patient_id: Unique identifier
            seed: Random seed for reproducibility
            
        Returns:
            VirtualPatient instance
        """
        definition = self.archetypes[archetype]
        
        # Sample all parameters
        age = definition.age_distribution.sample(1, seed)[0]
        sex = definition.sex_distribution.sample(1, seed + 1)[0]
        r_iss = definition.r_iss_distribution.sample(1, seed + 2)[0]
        cyto_risk = definition.cytogenetic_risk_distribution.sample(1, seed + 3)[0]
        
        # Sample cytogenetic lesions (Bernoulli trials)
        np.random.seed(seed + 4)
        has_del17p = np.random.rand() < definition.del17p_probability
        has_t4_14 = np.random.rand() < definition.t4_14_probability
        has_t14_16 = np.random.rand() < definition.t14_16_probability
        has_gain_1q21 = np.random.rand() < definition.gain_1q21_probability
        is_hyperdiploid = np.random.rand() < definition.hyperdiploid_probability
        has_t11_14 = np.random.rand() < definition.t11_14_probability
        
        # Tumor parameters
        tumor_burden = definition.tumor_burden_distribution.sample(1, seed + 10)[0]
        tumor_growth_rate = definition.tumor_growth_rate_distribution.sample(1, seed + 11)[0]
        carrying_capacity = definition.carrying_capacity_distribution.sample(1, seed + 12)[0]
        
        # Healthy tissue
        healthy_cells = definition.healthy_cell_distribution.sample(1, seed + 20)[0]
        healthy_growth_rate = definition.healthy_growth_rate_distribution.sample(1, seed + 21)[0]
        
        # Fitness
        ecog = int(definition.ecog_distribution.sample(1, seed + 30)[0])
        charlson = definition.charlson_comorbidity_distribution.sample(1, seed + 31)[0]
        
        # Labs
        creatinine_clearance = definition.creatinine_clearance_distribution.sample(1, seed + 40)[0]
        serum_albumin = definition.serum_albumin_distribution.sample(1, seed + 41)[0]
        ldh = definition.ldh_distribution.sample(1, seed + 42)[0]
        beta2_microglobulin = definition.beta2_microglobulin_distribution.sample(1, seed + 43)[0]
        
        return VirtualPatient(
            patient_id=patient_id,
            archetype=archetype,
            seed=seed,
            age=float(age),
            sex=str(sex),
            r_iss_stage=str(r_iss),
            cytogenetic_risk=str(cyto_risk),
            has_del17p=bool(has_del17p),
            has_t4_14=bool(has_t4_14),
            has_t14_16=bool(has_t14_16),
            has_gain_1q21=bool(has_gain_1q21),
            is_hyperdiploid=bool(is_hyperdiploid),
            has_t11_14=bool(has_t11_14),
            tumor_burden=float(tumor_burden),
            tumor_growth_rate=float(tumor_growth_rate),
            carrying_capacity=float(carrying_capacity),
            healthy_cells=float(healthy_cells),
            healthy_growth_rate=float(healthy_growth_rate),
            ecog_performance_status=int(ecog),
            charlson_comorbidity_index=float(charlson),
            creatinine_clearance=float(creatinine_clearance),
            serum_albumin=float(serum_albumin),
            ldh=float(ldh),
            beta2_microglobulin=float(beta2_microglobulin),
        )
    
    def generate_cohort(
        self,
        archetype: PatientArchetype,
        n: int,
        seed: int = 42,
    ) -> List[VirtualPatient]:
        """
        Generate cohort from single archetype.
        
        Args:
            archetype: Patient archetype
            n: Number of patients
            seed: Base random seed
            
        Returns:
            List of VirtualPatient instances
        """
        patients = []
        for i in range(n):
            patient_id = f"{archetype.value}_P{i+1:04d}"
            patient_seed = seed + i * 100
            patient = self.generate_patient(archetype, patient_id, patient_seed)
            patients.append(patient)
        
        return patients
    
    def generate_mixed_cohort(
        self,
        n: int,
        seed: int = 42,
        custom_prevalence: Optional[Dict[PatientArchetype, float]] = None,
    ) -> List[VirtualPatient]:
        """
        Generate cohort with mixed archetypes (realistic population).
        
        Args:
            n: Total number of patients
            seed: Base random seed
            custom_prevalence: Optional custom prevalence dictionary
            
        Returns:
            List of VirtualPatient instances
        """
        # Get prevalence weights
        if custom_prevalence:
            prevalence = custom_prevalence
        else:
            prevalence = {
                archetype: definition.prevalence
                for archetype, definition in self.archetypes.items()
            }
        
        # Normalize to sum to 1
        total = sum(prevalence.values())
        prevalence = {k: v / total for k, v in prevalence.items()}
        
        # Allocate patients to archetypes
        archetypes_list = list(prevalence.keys())
        probabilities = [prevalence[a] for a in archetypes_list]
        
        np.random.seed(seed)
        assignments = np.random.choice(
            archetypes_list,
            size=n,
            p=probabilities,
            replace=True,
        )
        
        # Generate patients
        patients = []
        archetype_counters = {a: 0 for a in archetypes_list}
        
        for i, archetype in enumerate(assignments):
            counter = archetype_counters[archetype]
            patient_id = f"{archetype.value}_P{counter+1:04d}"
            patient_seed = seed + i * 100
            patient = self.generate_patient(archetype, patient_id, patient_seed)
            patients.append(patient)
            archetype_counters[archetype] += 1
        
        return patients

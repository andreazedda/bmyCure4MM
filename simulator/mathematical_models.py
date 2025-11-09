"""
Mathematical Models for Multiple Myeloma Treatment Simulation.

This module defines the precise mathematical equations and algorithms used
for tumor growth, drug pharmacokinetics/pharmacodynamics, and treatment response.

All models are based on published literature and clinical trial data.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Callable, Optional
import numpy as np
from scipy.integrate import odeint


# =============================================================================
# CORE MATHEMATICAL MODEL DEFINITIONS
# =============================================================================

@dataclass
class TumorGrowthModel:
    """
    Gompertzian tumor growth model with carrying capacity.
    
    Mathematical Formulation:
    -------------------------
    dT/dt = r_T * T * ln(K_T / T)
    
    Where:
        T = tumor cell count
        r_T = intrinsic tumor growth rate [1/day]
        K_T = carrying capacity (maximum tumor burden)
        
    Reference:
        Norton L. (1988). A Gompertzian model of human breast cancer growth.
        Cancer Research, 48(24 Pt 1), 7067-7071.
    
    Parameters:
        r_tumor: Intrinsic growth rate (typical range: 0.01-0.05 /day)
        K_tumor: Carrying capacity (typical: 1e12 cells)
        
    Notes:
        - Gompertzian growth accounts for nutrient/space limitations
        - More realistic than exponential growth for large tumors
        - Self-limiting at carrying capacity
    """
    r_tumor: float = 0.023  # /day - typical MM growth rate
    K_tumor: float = 1.0e12  # cells - maximum sustainable burden
    
    def growth_rate(self, tumor_cells: float) -> float:
        """
        Calculate instantaneous growth rate.
        
        Returns:
            dT/dt in cells/day
        """
        if tumor_cells <= 0 or tumor_cells >= self.K_tumor:
            return 0.0
        return self.r_tumor * tumor_cells * np.log(self.K_tumor / tumor_cells)


@dataclass
class HealthyCellModel:
    """
    Healthy plasma cell dynamics with homeostatic regulation.
    
    Mathematical Formulation:
    -------------------------
    dH/dt = r_H * H * (1 - H / K_H)
    
    Where:
        H = healthy plasma cell count
        r_H = healthy cell renewal rate [1/day]
        K_H = homeostatic equilibrium level
        
    Reference:
        Mackey MC. (1978). Unified hypothesis for the origin of aplastic anemia
        and periodic hematopoiesis. Blood, 51(5), 941-956.
    
    Parameters:
        r_healthy: Renewal rate (typical: 0.01-0.02 /day)
        K_healthy: Equilibrium level (typical: 5e11 cells)
    """
    r_healthy: float = 0.015  # /day - plasma cell renewal rate
    K_healthy: float = 5.0e11  # cells - normal plasma cell count
    
    def renewal_rate(self, healthy_cells: float) -> float:
        """
        Calculate instantaneous renewal rate.
        
        Returns:
            dH/dt in cells/day
        """
        if healthy_cells <= 0:
            return 0.0
        return self.r_healthy * healthy_cells * (1.0 - healthy_cells / self.K_healthy)


@dataclass
class PharmacokineticModel:
    """
    One-compartment PK model with first-order elimination.
    
    Mathematical Formulation:
    -------------------------
    dC/dt = -k_e * C + f_dose(t)
    
    Where:
        C = drug concentration [mg/L]
        k_e = elimination rate constant [1/hour]
        f_dose(t) = dosing function (bolus or infusion)
        
    Derived parameters:
        - Half-life: t_1/2 = ln(2) / k_e
        - Volume of distribution: Vd [L]
        - Clearance: CL = k_e * Vd [L/hour]
        
    Reference:
        Gabrielsson J, Weiner D. (2012). Pharmacokinetic and Pharmacodynamic
        Data Analysis: Concepts and Applications. 5th ed.
    
    Parameters:
        half_life: Drug half-life [hours]
        Vd: Volume of distribution [L]
    """
    half_life: float  # hours
    Vd: float  # liters
    
    @property
    def elimination_rate(self) -> float:
        """
        Calculate elimination rate constant.
        
        Returns:
            k_e in 1/hour
        """
        return np.log(2) / self.half_life
    
    @property
    def clearance(self) -> float:
        """
        Calculate total body clearance.
        
        Returns:
            CL in L/hour
        """
        return self.elimination_rate * self.Vd
    
    def concentration_at_time(
        self,
        dose: float,
        time: float,
        previous_concentration: float = 0.0,
    ) -> float:
        """
        Calculate concentration at time t after dose.
        
        For bolus administration:
            C(t) = (Dose / Vd) * exp(-k_e * t) + C_prev * exp(-k_e * t)
        
        Args:
            dose: Administered dose [mg]
            time: Time since dose [hours]
            previous_concentration: Baseline concentration [mg/L]
            
        Returns:
            Concentration in mg/L
        """
        k_e = self.elimination_rate
        C_dose = (dose / self.Vd) * np.exp(-k_e * time)
        C_residual = previous_concentration * np.exp(-k_e * time)
        return C_dose + C_residual


@dataclass
class PharmacodynamicModel:
    """
    Emax model for drug effect on tumor cells.
    
    Mathematical Formulation:
    -------------------------
    Effect = E_max * C^n / (EC50^n + C^n)
    
    Where:
        E_max = maximum achievable effect (0-1, dimensionless)
        C = drug concentration [mg/L]
        EC50 = concentration for 50% effect [mg/L]
        n = Hill coefficient (slope factor, typically 1-4)
        
    Cell kill rate:
        k_kill = Effect * λ_max
        
    Where:
        λ_max = maximum kill rate [1/day]
        
    Reference:
        Holford NH, Sheiner LB. (1981). Understanding the dose-effect relationship.
        Clinical Pharmacokinetics, 6(6), 429-453.
    
    Parameters:
        Emax: Maximum effect (0-1)
        EC50: Half-maximal concentration [mg/L]
        Hill: Hill coefficient (default 1)
        lambda_max: Maximum kill rate [1/day]
    """
    Emax: float  # 0-1 dimensionless
    EC50: float  # mg/L
    Hill: float = 1.0  # dimensionless
    lambda_max: float = 0.2  # 1/day - maximum kill rate
    
    def effect(self, concentration: float) -> float:
        """
        Calculate drug effect at given concentration.
        
        Args:
            concentration: Drug concentration [mg/L]
            
        Returns:
            Effect (0-1, dimensionless)
        """
        if concentration <= 0:
            return 0.0
        numerator = self.Emax * (concentration ** self.Hill)
        denominator = (self.EC50 ** self.Hill) + (concentration ** self.Hill)
        return numerator / denominator
    
    def kill_rate(self, concentration: float) -> float:
        """
        Calculate tumor cell kill rate.
        
        Args:
            concentration: Drug concentration [mg/L]
            
        Returns:
            Kill rate in 1/day
        """
        return self.effect(concentration) * self.lambda_max


@dataclass
class DrugInteractionModel:
    """
    Pharmacodynamic drug interaction model.
    
    Mathematical Formulation (Greco et al.):
    -----------------------------------------
    For two drugs A and B:
        Combined effect = E_A + E_B + α * E_A * E_B
        
    Where:
        E_A, E_B = individual drug effects (0-1)
        α = interaction parameter:
            α > 0: Synergistic interaction
            α = 0: Additive interaction
            α < 0: Antagonistic interaction
            
    For n drugs, pairwise interactions:
        E_total = Σ E_i + Σ Σ α_ij * E_i * E_j (i < j)
        
    Reference:
        Greco WR, Bravo G, Parsons JC. (1995). The search for synergy:
        a critical review from a response surface perspective.
        Pharmacological Reviews, 47(2), 331-385.
    
    Parameters:
        interaction_matrix: n×n matrix of interaction coefficients α_ij
    """
    interaction_matrix: np.ndarray
    
    def combined_effect(self, individual_effects: np.ndarray) -> float:
        """
        Calculate combined drug effect considering interactions.
        
        Args:
            individual_effects: Array of individual drug effects (0-1)
            
        Returns:
            Combined effect (0-1)
        """
        n = len(individual_effects)
        
        # Additive effects
        total_effect = float(np.sum(individual_effects))
        
        # Pairwise interactions
        for i in range(n):
            for j in range(i + 1, n):
                alpha_ij = self.interaction_matrix[i, j]
                interaction = alpha_ij * individual_effects[i] * individual_effects[j]
                total_effect += interaction
        
        # Cap at 1.0 (complete kill)
        return min(total_effect, 1.0)


# =============================================================================
# IMMUNE SYSTEM MODELING
# =============================================================================

@dataclass
class ImmuneResponseModel:
    """
    Simplified immune surveillance model.
    
    Mathematical Formulation:
    -------------------------
    Immune kill rate: k_immune = η * I * T / (T + K_I)
    
    Where:
        η = immune efficiency parameter
        I = immune competence index (0-1)
        T = tumor cell count
        K_I = tumor burden for half-maximal immune response
        
    Reference:
        de Pillis LG, Radunskaya AE. (2001). A mathematical tumor model
        with immune resistance and drug therapy. Computational and
        Mathematical Methods in Medicine, 3(2), 79-100.
    
    Parameters:
        immune_efficiency: η parameter (typical: 1e-9 to 1e-11)
        immune_competence: I index (0=no immunity, 1=fully competent)
        half_saturation: K_I (tumor burden for half-max response)
    """
    immune_efficiency: float = 1.0e-10  # 1/(cells·day)
    immune_competence: float = 1.0  # 0-1 dimensionless
    half_saturation: float = 1.0e9  # cells
    
    def kill_rate(self, tumor_cells: float) -> float:
        """
        Calculate immune-mediated tumor cell kill rate.
        
        Args:
            tumor_cells: Current tumor burden
            
        Returns:
            Kill rate in 1/day
        """
        if tumor_cells <= 0:
            return 0.0
        
        numerator = self.immune_efficiency * self.immune_competence * tumor_cells
        denominator = tumor_cells + self.half_saturation
        return numerator / denominator


# =============================================================================
# COMPLETE SYSTEM ODE
# =============================================================================

class MultipleMyelomaODE:
    """
    Complete ODE system for MM treatment simulation.
    
    Mathematical Formulation:
    -------------------------
    System of ODEs:
    
    1. Tumor cells:
       dT/dt = r_T * T * ln(K_T / T) 
               - Σ k_kill,i(C_i) * T
               - k_immune(T, I) * T
    
    2. Healthy plasma cells:
       dH/dt = r_H * H * (1 - H / K_H)
               - Σ ω_i * k_kill,i(C_i) * H
    
    3. Drug concentrations (for each drug i):
       dC_i/dt = -k_e,i * C_i + f_dose,i(t)
    
    Where:
        ω_i = toxicity weight (0-1, fraction of tumor effect on healthy cells)
        
    Boundary conditions:
        T(0) = T_0 (initial tumor burden)
        H(0) = H_0 (initial healthy cell count)
        C_i(0) = 0 (no drug at baseline)
        
    Parameters:
        All parameters from component models above
    """
    
    def __init__(
        self,
        tumor_model: TumorGrowthModel,
        healthy_model: HealthyCellModel,
        pk_models: Dict[str, PharmacokineticModel],
        pd_models: Dict[str, PharmacodynamicModel],
        interaction_model: DrugInteractionModel,
        immune_model: ImmuneResponseModel,
        toxicity_weights: Dict[str, float],
        dose_functions: Dict[str, Callable],
    ):
        self.tumor_model = tumor_model
        self.healthy_model = healthy_model
        self.pk_models = pk_models
        self.pd_models = pd_models
        self.interaction_model = interaction_model
        self.immune_model = immune_model
        self.toxicity_weights = toxicity_weights
        self.dose_functions = dose_functions
        
        # Build state vector index mapping
        self.idx = {"tumor": 0, "healthy": 1}
        drug_names = sorted(pk_models.keys())
        for i, drug in enumerate(drug_names):
            self.idx[drug] = 2 + i
        self.n_states = 2 + len(drug_names)
    
    def derivatives(self, state: np.ndarray, t: float) -> np.ndarray:
        """
        Calculate derivatives for ODE system.
        
        Args:
            state: Current state vector [T, H, C_1, ..., C_n]
            t: Current time [days]
            
        Returns:
            Derivatives dState/dt
        """
        # Extract state variables
        tumor = state[self.idx["tumor"]]
        healthy = state[self.idx["healthy"]]
        
        # Tumor growth
        dT_growth = self.tumor_model.growth_rate(tumor)
        
        # Healthy cell renewal
        dH_renewal = self.healthy_model.renewal_rate(healthy)
        
        # Drug effects
        drug_effects = []
        tumor_kill_total = 0.0
        healthy_kill_total = 0.0
        
        drug_names = sorted(self.pk_models.keys())
        dC_dt = {}
        
        for drug in drug_names:
            concentration = state[self.idx[drug]]
            
            # PD effect
            effect = self.pd_models[drug].effect(concentration)
            drug_effects.append(effect)
            
            # Individual kill rate
            kill_rate = self.pd_models[drug].kill_rate(concentration)
            tumor_kill_total += kill_rate
            
            # Healthy cell toxicity
            toxicity_weight = self.toxicity_weights.get(drug, 0.1)
            healthy_kill_total += toxicity_weight * kill_rate
            
            # PK derivative
            k_e = self.pk_models[drug].elimination_rate / 24.0  # Convert to 1/day
            dose_rate = self.dose_functions[drug](t) if drug in self.dose_functions else 0.0
            dC_dt[drug] = -k_e * concentration + dose_rate
        
        # Combined drug effect (with interactions)
        if len(drug_effects) > 1:
            drug_effects_array = np.array(drug_effects)
            combined_effect = self.interaction_model.combined_effect(drug_effects_array)
            # Adjust kill rates by interaction factor
            interaction_factor = combined_effect / max(sum(drug_effects), 1e-9)
            tumor_kill_total *= interaction_factor
            healthy_kill_total *= interaction_factor
        
        # Immune-mediated kill
        immune_kill = self.immune_model.kill_rate(tumor)
        
        # Assemble derivatives
        dT_dt = dT_growth - tumor_kill_total * tumor - immune_kill * tumor
        dH_dt = dH_renewal - healthy_kill_total * healthy
        
        # Build derivative vector
        derivatives = np.zeros(self.n_states)
        derivatives[self.idx["tumor"]] = dT_dt
        derivatives[self.idx["healthy"]] = dH_dt
        
        for drug in drug_names:
            derivatives[self.idx[drug]] = dC_dt[drug]
        
        return derivatives
    
    def simulate(
        self,
        initial_state: np.ndarray,
        time_points: np.ndarray,
    ) -> np.ndarray:
        """
        Simulate ODE system over time.
        
        Args:
            initial_state: Initial conditions [T_0, H_0, 0, ..., 0]
            time_points: Time points for solution [days]
            
        Returns:
            Solution array (n_timepoints × n_states)
        """
        solution = odeint(self.derivatives, initial_state, time_points)
        return solution


# =============================================================================
# CLINICAL OUTCOME METRICS
# =============================================================================

def calculate_response_category(tumor_reduction: float) -> str:
    """
    Classify treatment response per IMWG criteria.
    
    Response Categories:
    - CR (Complete Response): ≥95% reduction
    - VGPR (Very Good Partial Response): 90-95% reduction
    - PR (Partial Response): 50-90% reduction
    - SD (Stable Disease): -25% to +50% change
    - PD (Progressive Disease): >25% increase
    
    Reference:
        Kumar S, et al. (2016). International Myeloma Working Group
        consensus criteria for response. Lancet Oncol, 17(8), e328-e346.
    
    Args:
        tumor_reduction: Fractional tumor reduction (0-1)
        
    Returns:
        Response category code
    """
    if tumor_reduction >= 0.95:
        return "CR"
    elif tumor_reduction >= 0.90:
        return "VGPR"
    elif tumor_reduction >= 0.50:
        return "PR"
    elif tumor_reduction >= -0.25:
        return "SD"
    else:
        return "PD"


def calculate_toxicity_grade(healthy_loss: float) -> int:
    """
    Calculate toxicity grade based on healthy cell loss.
    
    CTCAE Grading (Common Terminology Criteria for Adverse Events):
    - Grade 0: <10% loss (none)
    - Grade 1: 10-25% loss (mild)
    - Grade 2: 25-50% loss (moderate)
    - Grade 3: 50-75% loss (severe)
    - Grade 4: >75% loss (life-threatening)
    
    Reference:
        NCI. (2017). Common Terminology Criteria for Adverse Events
        (CTCAE) Version 5.0.
    
    Args:
        healthy_loss: Fractional healthy cell loss (0-1)
        
    Returns:
        Toxicity grade (0-4)
    """
    if healthy_loss < 0.10:
        return 0
    elif healthy_loss < 0.25:
        return 1
    elif healthy_loss < 0.50:
        return 2
    elif healthy_loss < 0.75:
        return 3
    else:
        return 4


def calculate_progression_free_survival(
    trajectory: np.ndarray,
    time_points: np.ndarray,
    threshold: float = 1.25,
) -> Optional[float]:
    """
    Estimate progression-free survival (PFS) time.
    
    PFS Definition:
        Time from treatment start to disease progression (≥25% increase
        from nadir) or death.
        
    Args:
        trajectory: Tumor burden over time
        time_points: Corresponding time points [days]
        threshold: Progression threshold (1.25 = 25% increase)
        
    Returns:
        PFS time in days, or None if no progression observed
    """
    nadir_idx = np.argmin(trajectory)
    nadir_value = trajectory[nadir_idx]
    
    post_nadir = trajectory[nadir_idx + 1:]
    progression_threshold = nadir_value * threshold
    
    progression_indices = np.where(post_nadir >= progression_threshold)[0]
    
    if len(progression_indices) > 0:
        progression_idx = nadir_idx + 1 + progression_indices[0]
        return float(time_points[progression_idx])
    
    return None

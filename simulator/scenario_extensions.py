"""
Enhanced Scenario Model with Mathematical Framework Integration.

This module extends the Django Scenario model with:
- Automatic difficulty score calculation
- Virtual patient archetype assignment
- Mathematical model parameter storage
- Clinical outcome prediction

All additions maintain backward compatibility with existing models.
"""
from __future__ import annotations

from typing import Dict, Optional
from django.db import models
from django.core.exceptions import ValidationError

from .difficulty_scoring import (
    TumorBurdenScore,
    GrowthRateScore,
    CytogeneticScore,
    PatientFrailtyScore,
    StageScore,
    DifficultyScoreCalculator,
    RISSStagingSystem,
    estimate_response_probability,
    estimate_toxicity_risk,
    estimate_survival_metrics,
)
from .virtual_patients import (
    PatientArchetype,
    VirtualPatientGenerator,
    VirtualPatient,
)


# =============================================================================
# ENHANCED SCENARIO MODEL (Mixin Pattern)
# =============================================================================

class ScenarioMathematicalMixin(models.Model):
    """
    Mixin to add mathematical framework to Scenario model.
    
    This can be added to existing Scenario model via multi-table inheritance
    or as a separate model with OneToOneField relationship.
    
    Fields added:
    - difficulty_score: Calculated composite score (0-100)
    - difficulty_level: Human-readable level (Very Easy to Very Hard)
    - patient_archetype: Clinical phenotype category
    - mathematical_parameters: JSON storage for all simulation parameters
    - expected_outcomes: Predicted response rates, survival, toxicity
    """
    
    # Difficulty scoring
    difficulty_score = models.FloatField(
        null=True,
        blank=True,
        help_text="Composite difficulty score (0-100) based on clinical factors.",
    )
    difficulty_level = models.CharField(
        max_length=20,
        blank=True,
        choices=[
            ("Very Easy", "Very Easy"),
            ("Easy", "Easy"),
            ("Moderate", "Moderate"),
            ("Hard", "Hard"),
            ("Very Hard", "Very Hard"),
        ],
        help_text="Categorical difficulty level.",
    )
    
    # Patient archetype
    patient_archetype = models.CharField(
        max_length=50,
        blank=True,
        choices=[
            ("nd_standard", "Newly Diagnosed Standard Risk"),
            ("nd_high_risk", "Newly Diagnosed High Risk"),
            ("frail_elderly", "Frail Elderly"),
            ("relapsed_refractory", "Relapsed/Refractory"),
            ("smoldering", "Smoldering Myeloma"),
            ("transplant_eligible", "Transplant Eligible"),
            ("transplant_ineligible", "Transplant Ineligible"),
        ],
        help_text="Clinical archetype for virtual patient generation.",
    )
    
    # Mathematical parameters (JSON storage)
    mathematical_parameters = models.JSONField(
        default=dict,
        blank=True,
        help_text="Complete set of parameters for ODE simulation.",
    )
    
    # Expected outcomes (pre-calculated)
    expected_outcomes = models.JSONField(
        default=dict,
        blank=True,
        help_text="Predicted response rates, survival metrics, toxicity risk.",
    )
    
    # R-ISS staging
    r_iss_stage = models.CharField(
        max_length=10,
        blank=True,
        choices=[
            ("I", "R-ISS Stage I"),
            ("II", "R-ISS Stage II"),
            ("III", "R-ISS Stage III"),
        ],
        help_text="Revised International Staging System stage.",
    )
    
    # Cytogenetics
    has_del17p = models.BooleanField(default=False, help_text="Deletion 17p (TP53)")
    has_t4_14 = models.BooleanField(default=False, help_text="Translocation t(4;14)")
    has_t14_16 = models.BooleanField(default=False, help_text="Translocation t(14;16)")
    has_gain_1q21 = models.BooleanField(default=False, help_text="1q21 gain/amplification")
    is_hyperdiploid = models.BooleanField(default=False, help_text="Hyperdiploid karyotype")
    has_t11_14 = models.BooleanField(default=False, help_text="Translocation t(11;14)")
    
    # Patient characteristics
    patient_age = models.IntegerField(null=True, blank=True, help_text="Patient age (years)")
    ecog_performance_status = models.IntegerField(
        null=True,
        blank=True,
        choices=[(0, "0"), (1, "1"), (2, "2"), (3, "3"), (4, "4")],
        help_text="ECOG performance status (0=fully active, 4=bedbound)",
    )
    charlson_comorbidity_index = models.IntegerField(
        null=True,
        blank=True,
        help_text="Charlson Comorbidity Index",
    )
    
    # Laboratory values (for difficulty calculation)
    creatinine_clearance = models.FloatField(
        null=True,
        blank=True,
        help_text="Creatinine clearance (mL/min)",
    )
    serum_albumin = models.FloatField(
        null=True,
        blank=True,
        help_text="Serum albumin (g/dL)",
    )
    ldh = models.FloatField(
        null=True,
        blank=True,
        help_text="Lactate dehydrogenase (U/L)",
    )
    beta2_microglobulin = models.FloatField(
        null=True,
        blank=True,
        help_text="Beta-2 microglobulin (mg/L)",
    )
    
    # Tumor biology parameters
    tumor_cell_count = models.FloatField(
        null=True,
        blank=True,
        help_text="Initial tumor burden (cells)",
    )
    tumor_growth_rate = models.FloatField(
        null=True,
        blank=True,
        help_text="Tumor growth rate (1/day)",
    )
    
    class Meta:
        abstract = True
    
    def calculate_difficulty_score(self) -> float:
        """
        Calculate composite difficulty score from clinical parameters.
        
        Returns:
            Difficulty score (0-100)
        """
        # Validate required fields
        if not self._validate_difficulty_parameters():
            raise ValidationError(
                "Missing required parameters for difficulty calculation. "
                "Need: tumor_cell_count, tumor_growth_rate, r_iss_stage, "
                "patient_age, ecog_performance_status, charlson_comorbidity_index, "
                "creatinine_clearance, serum_albumin"
            )
        
        # Create component scores
        tumor_burden_score = TumorBurdenScore(tumor_cells=self.tumor_cell_count)
        growth_rate_score = GrowthRateScore(growth_rate=self.tumor_growth_rate)
        
        cytogenetic_score = CytogeneticScore(
            has_del17p=self.has_del17p,
            has_t4_14=self.has_t4_14,
            has_t14_16=self.has_t14_16,
            has_1q21_gain=self.has_gain_1q21,
            is_hyperdiploid=self.is_hyperdiploid,
            has_t11_14=self.has_t11_14,
        )
        
        frailty_score = PatientFrailtyScore(
            age=self.patient_age,
            ecog_performance_status=self.ecog_performance_status or 0,
            charlson_comorbidity_index=self.charlson_comorbidity_index or 0,
            creatinine_clearance=self.creatinine_clearance,
            serum_albumin=self.serum_albumin,
        )
        
        # Map R-ISS string to enum
        r_iss_map = {
            "I": RISSStagingSystem.STAGE_I,
            "II": RISSStagingSystem.STAGE_II,
            "III": RISSStagingSystem.STAGE_III,
        }
        r_iss_enum = r_iss_map.get(self.r_iss_stage, RISSStagingSystem.STAGE_II)
        stage_score = StageScore(r_iss_stage=r_iss_enum)
        
        # Calculate total
        calculator = DifficultyScoreCalculator(
            tumor_burden_score=tumor_burden_score,
            growth_rate_score=growth_rate_score,
            cytogenetic_score=cytogenetic_score,
            frailty_score=frailty_score,
            stage_score=stage_score,
        )
        
        total_score = calculator.compute_total()
        level = calculator.get_difficulty_level()
        
        # Store results
        self.difficulty_score = total_score
        self.difficulty_level = level
        
        # Store component breakdown
        breakdown = calculator.get_component_breakdown()
        if not self.mathematical_parameters:
            self.mathematical_parameters = {}
        self.mathematical_parameters["difficulty_breakdown"] = breakdown
        
        return total_score
    
    def calculate_expected_outcomes(self) -> Dict:
        """
        Predict expected clinical outcomes based on difficulty.
        
        Returns:
            Dictionary with response rates, survival, toxicity predictions
        """
        if self.difficulty_score is None:
            self.calculate_difficulty_score()
        
        frailty_component = self.mathematical_parameters.get(
            "difficulty_breakdown", {}
        ).get("frailty", 7.5)
        
        outcomes = {
            "response_probabilities": estimate_response_probability(self.difficulty_score),
            "toxicity_risk": estimate_toxicity_risk(self.difficulty_score, frailty_component),
            "survival_estimates": estimate_survival_metrics(self.difficulty_score),
        }
        
        self.expected_outcomes = outcomes
        return outcomes
    
    def generate_virtual_patient(self, seed: int = 42) -> VirtualPatient:
        """
        Generate a VirtualPatient instance matching this scenario.
        
        Args:
            seed: Random seed for reproducibility
            
        Returns:
            VirtualPatient with parameters from this scenario
        """
        if not self.patient_archetype:
            raise ValidationError("patient_archetype must be set to generate virtual patient")
        
        # Map string to enum
        archetype_map = {
            "nd_standard": PatientArchetype.NEWLY_DIAGNOSED_STANDARD_RISK,
            "nd_high_risk": PatientArchetype.NEWLY_DIAGNOSED_HIGH_RISK,
            "frail_elderly": PatientArchetype.FRAIL_ELDERLY,
            "relapsed_refractory": PatientArchetype.RELAPSED_REFRACTORY,
            "smoldering": PatientArchetype.SMOLDERING_MYELOMA,
            "transplant_eligible": PatientArchetype.TRANSPLANT_ELIGIBLE,
            "transplant_ineligible": PatientArchetype.TRANSPLANT_INELIGIBLE,
        }
        archetype_enum = archetype_map.get(self.patient_archetype)
        if archetype_enum is None:
            raise ValidationError(f"Unknown archetype: {self.patient_archetype}")
        
        generator = VirtualPatientGenerator()
        patient = generator.generate_patient(
            archetype=archetype_enum,
            patient_id=f"scenario_{self.pk}",
            seed=seed,
        )
        
        return patient
    
    def populate_from_virtual_patient(self, patient: VirtualPatient) -> None:
        """
        Populate scenario fields from a VirtualPatient instance.
        
        Args:
            patient: VirtualPatient to copy parameters from
        """
        # Demographics
        self.patient_age = int(patient.age)
        
        # Disease stage
        self.r_iss_stage = patient.r_iss_stage
        
        # Cytogenetics
        self.has_del17p = patient.has_del17p
        self.has_t4_14 = patient.has_t4_14
        self.has_t14_16 = patient.has_t14_16
        self.has_gain_1q21 = patient.has_gain_1q21
        self.is_hyperdiploid = patient.is_hyperdiploid
        self.has_t11_14 = patient.has_t11_14
        
        # Tumor biology
        self.tumor_cell_count = patient.tumor_burden
        self.tumor_growth_rate = patient.tumor_growth_rate
        
        # Patient fitness
        self.ecog_performance_status = patient.ecog_performance_status
        self.charlson_comorbidity_index = int(patient.charlson_comorbidity_index)
        
        # Laboratory
        self.creatinine_clearance = patient.creatinine_clearance
        self.serum_albumin = patient.serum_albumin
        self.ldh = patient.ldh
        self.beta2_microglobulin = patient.beta2_microglobulin
        
        # Set archetype
        self.patient_archetype = patient.archetype.value
        
        # Store all parameters
        self.mathematical_parameters = patient.to_dict()
    
    def _validate_difficulty_parameters(self) -> bool:
        """Check if all required parameters are present for difficulty calculation."""
        required = [
            self.tumor_cell_count,
            self.tumor_growth_rate,
            self.r_iss_stage,
            self.patient_age,
            self.creatinine_clearance,
            self.serum_albumin,
        ]
        return all(x is not None for x in required)
    
    def get_difficulty_display_html(self) -> str:
        """
        Generate HTML display of difficulty score with visual indicators.
        
        Returns:
            HTML string with styled difficulty display
        """
        if self.difficulty_score is None:
            return "<span class='badge badge-secondary'>Not Calculated</span>"
        
        score = self.difficulty_score
        level = self.difficulty_level or "Moderate"
        
        # Color coding
        color_map = {
            "Very Easy": "success",
            "Easy": "info",
            "Moderate": "warning",
            "Hard": "danger",
            "Very Hard": "dark",
        }
        badge_color = color_map.get(level, "secondary")
        
        html = f"""
        <div class="difficulty-display">
            <span class="badge badge-{badge_color}" style="font-size: 1.1em;">
                {level}
            </span>
            <span class="difficulty-score" style="margin-left: 10px;">
                Score: {score:.1f}/100
            </span>
        </div>
        """
        
        return html
    
    def get_difficulty_breakdown_table(self) -> str:
        """
        Generate HTML table showing component scores.
        
        Returns:
            HTML table with breakdown
        """
        if not self.mathematical_parameters.get("difficulty_breakdown"):
            return ""
        
        breakdown = self.mathematical_parameters["difficulty_breakdown"]
        
        rows = []
        component_names = {
            "tumor_burden": "Tumor Burden",
            "growth_rate": "Growth Rate",
            "cytogenetics": "Cytogenetics",
            "frailty": "Patient Frailty",
            "stage": "Disease Stage",
        }
        max_points = {
            "tumor_burden": 25,
            "growth_rate": 20,
            "cytogenetics": 25,
            "frailty": 15,
            "stage": 15,
        }
        
        for key, name in component_names.items():
            score = breakdown.get(key, 0)
            max_pts = max_points[key]
            percentage = (score / max_pts) * 100
            
            rows.append(f"""
                <tr>
                    <td>{name}</td>
                    <td>{score:.1f}</td>
                    <td>{max_pts}</td>
                    <td>
                        <div class="progress" style="height: 20px;">
                            <div class="progress-bar" role="progressbar" 
                                 style="width: {percentage}%">
                                {percentage:.0f}%
                            </div>
                        </div>
                    </td>
                </tr>
            """)
        
        table = f"""
        <table class="table table-sm">
            <thead>
                <tr>
                    <th>Component</th>
                    <th>Score</th>
                    <th>Max</th>
                    <th>Percentage</th>
                </tr>
            </thead>
            <tbody>
                {"".join(rows)}
            </tbody>
        </table>
        """
        
        return table


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def create_scenario_from_archetype(
    archetype: PatientArchetype,
    title: str,
    seed: int = 42,
    **kwargs,
) -> "Scenario":
    """
    Create a Scenario instance from a patient archetype.
    
    This is a factory function that generates a complete scenario with
    all mathematical parameters automatically populated.
    
    Args:
        archetype: Patient archetype to base scenario on
        title: Scenario title
        seed: Random seed for patient generation
        **kwargs: Additional fields to set on Scenario
        
    Returns:
        Scenario instance (not saved to database)
        
    Example:
        scenario = create_scenario_from_archetype(
            archetype=PatientArchetype.NEWLY_DIAGNOSED_HIGH_RISK,
            title="High-Risk NDMM Case",
            summary="68yo M with del(17p) and elevated LDH...",
            clinical_stage="newly_diagnosed",
        )
        scenario.save()
    """
    from .models import Scenario  # Avoid circular import
    
    # Generate virtual patient
    generator = VirtualPatientGenerator()
    patient = generator.generate_patient(
        archetype=archetype,
        patient_id=f"scenario_{title}",
        seed=seed,
    )
    
    # Create scenario
    scenario = Scenario(title=title, **kwargs)
    
    # Populate from virtual patient
    if hasattr(scenario, 'populate_from_virtual_patient'):
        scenario.populate_from_virtual_patient(patient)
    
    # Calculate difficulty
    if hasattr(scenario, 'calculate_difficulty_score'):
        scenario.calculate_difficulty_score()
        scenario.calculate_expected_outcomes()
    
    return scenario


def batch_calculate_difficulty(queryset) -> int:
    """
    Calculate difficulty scores for multiple scenarios.
    
    Args:
        queryset: QuerySet of Scenario instances
        
    Returns:
        Number of scenarios updated
    """
    count = 0
    for scenario in queryset:
        if hasattr(scenario, 'calculate_difficulty_score'):
            try:
                scenario.calculate_difficulty_score()
                scenario.calculate_expected_outcomes()
                scenario.save()
                count += 1
            except ValidationError:
                # Skip scenarios with missing data
                continue
    
    return count

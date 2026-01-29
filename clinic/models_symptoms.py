"""
SymptomAssessment model for preformatted symptom tracking.

This module provides structured symptom capture aligned with clinical scales:
- Bone pain: Numeric Rating Scale (NRS) 0-10
- Fatigue: FACIT-Fatigue scale (0-52, higher = less fatigue)
- Neuropathy: CTCAE v5.0 Grade (0-4)
- CRAB symptoms: Checkboxes for Calcium, Renal, Anemia, Bone
- Infections: Structured recording
"""

from __future__ import annotations

from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator


class SymptomAssessment(models.Model):
    """
    Preformatted symptom capture for MM patients.
    
    Uses validated clinical scales for objective tracking:
    - NRS for pain (0-10)
    - FACIT-F for fatigue (0-52)
    - CTCAE for neuropathy grading (0-4)
    - CRAB symptoms as boolean flags
    """
    
    # Relationship
    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="symptom_assessments"
    )
    assessment_date = models.DateField(auto_now_add=True)
    
    # ── Pain Assessment (NRS 0-10) ──────────────────────────────────────────
    bone_pain_nrs = models.PositiveSmallIntegerField(
        "Bone Pain (NRS 0-10)",
        validators=[MinValueValidator(0), MaxValueValidator(10)],
        null=True,
        blank=True,
        help_text="Numeric Rating Scale: 0=No pain, 10=Worst imaginable"
    )
    bone_pain_location = models.CharField(
        "Pain Location",
        max_length=255,
        blank=True,
        help_text="e.g., lumbar spine, ribs, pelvis"
    )
    
    # ── Fatigue (FACIT-Fatigue) ─────────────────────────────────────────────
    FACIT_CHOICES = [
        (0, "0 - Not at all"),
        (1, "1 - A little bit"),
        (2, "2 - Somewhat"),
        (3, "3 - Quite a bit"),
        (4, "4 - Very much"),
    ]
    
    # FACIT-F items (13 questions, each 0-4, sum = 0-52)
    # Higher score = LESS fatigue (better quality of life)
    fatigue_total = models.PositiveSmallIntegerField(
        "FACIT-Fatigue Score (0-52)",
        validators=[MinValueValidator(0), MaxValueValidator(52)],
        null=True,
        blank=True,
        help_text="Total FACIT-F score: 0=Severe fatigue, 52=No fatigue"
    )
    
    # Quick fatigue questions for simplified assessment
    fatigue_feel_tired = models.PositiveSmallIntegerField(
        "I feel fatigued",
        choices=FACIT_CHOICES,
        null=True,
        blank=True
    )
    fatigue_feel_weak = models.PositiveSmallIntegerField(
        "I feel weak all over",
        choices=FACIT_CHOICES,
        null=True,
        blank=True
    )
    fatigue_feel_listless = models.PositiveSmallIntegerField(
        "I feel listless (washed out)",
        choices=FACIT_CHOICES,
        null=True,
        blank=True
    )
    fatigue_trouble_starting = models.PositiveSmallIntegerField(
        "I have trouble starting things",
        choices=FACIT_CHOICES,
        null=True,
        blank=True
    )
    fatigue_trouble_finishing = models.PositiveSmallIntegerField(
        "I have trouble finishing things",
        choices=FACIT_CHOICES,
        null=True,
        blank=True
    )
    
    # ── Neuropathy (CTCAE v5.0 Grading) ─────────────────────────────────────
    CTCAE_NEUROPATHY_CHOICES = [
        (0, "Grade 0 - None"),
        (1, "Grade 1 - Asymptomatic; clinical or diagnostic findings only"),
        (2, "Grade 2 - Moderate; limiting instrumental ADL"),
        (3, "Grade 3 - Severe; limiting self-care ADL"),
        (4, "Grade 4 - Life-threatening; urgent intervention indicated"),
    ]
    
    neuropathy_sensory = models.PositiveSmallIntegerField(
        "Sensory Neuropathy (CTCAE)",
        choices=CTCAE_NEUROPATHY_CHOICES,
        default=0,
        help_text="Numbness, tingling, paresthesia"
    )
    neuropathy_motor = models.PositiveSmallIntegerField(
        "Motor Neuropathy (CTCAE)",
        choices=CTCAE_NEUROPATHY_CHOICES,
        default=0,
        help_text="Weakness, difficulty with fine motor skills"
    )
    neuropathy_pain = models.PositiveSmallIntegerField(
        "Neuropathic Pain (CTCAE)",
        choices=CTCAE_NEUROPATHY_CHOICES,
        default=0,
        help_text="Burning, shooting, electric-shock pain"
    )
    
    # ── CRAB Symptoms ───────────────────────────────────────────────────────
    # C - Calcium elevation
    crab_hypercalcemia = models.BooleanField(
        "Hypercalcemia (C)",
        default=False,
        help_text="Serum calcium >11 mg/dL or >1 mg/dL above ULN"
    )
    crab_hypercalcemia_symptoms = models.CharField(
        "Hypercalcemia Symptoms",
        max_length=255,
        blank=True,
        help_text="e.g., confusion, constipation, polyuria, nausea"
    )
    
    # R - Renal insufficiency
    crab_renal = models.BooleanField(
        "Renal Insufficiency (R)",
        default=False,
        help_text="Creatinine >2 mg/dL or CrCl <40 mL/min"
    )
    crab_renal_symptoms = models.CharField(
        "Renal Symptoms",
        max_length=255,
        blank=True,
        help_text="e.g., edema, decreased urine output, fatigue"
    )
    
    # A - Anemia
    crab_anemia = models.BooleanField(
        "Anemia (A)",
        default=False,
        help_text="Hemoglobin <10 g/dL or >2 g/dL below LLN"
    )
    crab_anemia_symptoms = models.CharField(
        "Anemia Symptoms",
        max_length=255,
        blank=True,
        help_text="e.g., fatigue, dyspnea, pallor, dizziness"
    )
    
    # B - Bone lesions
    crab_bone = models.BooleanField(
        "Bone Lesions (B)",
        default=False,
        help_text="Lytic lesions on imaging or osteoporosis with fractures"
    )
    crab_bone_symptoms = models.CharField(
        "Bone Symptoms",
        max_length=255,
        blank=True,
        help_text="e.g., fractures, vertebral compression, bone pain"
    )
    
    # ── Infections ──────────────────────────────────────────────────────────
    INFECTION_TYPE_CHOICES = [
        ("none", "No active infection"),
        ("respiratory", "Respiratory (pneumonia, bronchitis, URI)"),
        ("urinary", "Urinary tract infection"),
        ("skin", "Skin/soft tissue infection"),
        ("bloodstream", "Bloodstream infection (bacteremia/sepsis)"),
        ("viral", "Viral infection (COVID, herpes, shingles)"),
        ("fungal", "Fungal infection"),
        ("other", "Other infection"),
    ]
    
    infection_present = models.BooleanField(
        "Active Infection",
        default=False
    )
    infection_type = models.CharField(
        "Infection Type",
        max_length=32,
        choices=INFECTION_TYPE_CHOICES,
        default="none"
    )
    infection_organism = models.CharField(
        "Organism (if known)",
        max_length=128,
        blank=True,
        help_text="e.g., E. coli, S. pneumoniae, COVID-19"
    )
    infection_hospitalized = models.BooleanField(
        "Required Hospitalization",
        default=False
    )
    
    # ── Performance Status (ECOG) ───────────────────────────────────────────
    ECOG_CHOICES = [
        (0, "0 - Fully active, no restrictions"),
        (1, "1 - Restricted in strenuous activity, ambulatory"),
        (2, "2 - Ambulatory, capable of self-care, unable to work"),
        (3, "3 - Limited self-care, confined to bed/chair >50% of day"),
        (4, "4 - Completely disabled, confined to bed/chair"),
    ]
    
    ecog_status = models.PositiveSmallIntegerField(
        "ECOG Performance Status",
        choices=ECOG_CHOICES,
        default=0
    )
    
    # ── Additional Symptoms ─────────────────────────────────────────────────
    appetite_loss = models.BooleanField("Appetite Loss", default=False)
    weight_loss_pct = models.DecimalField(
        "Weight Loss (%)",
        max_digits=4,
        decimal_places=1,
        null=True,
        blank=True,
        help_text="Percentage weight loss in last 3 months"
    )
    night_sweats = models.BooleanField("Night Sweats", default=False)
    bruising_bleeding = models.BooleanField("Easy Bruising/Bleeding", default=False)
    
    # ── Notes ───────────────────────────────────────────────────────────────
    additional_notes = models.TextField(
        "Additional Notes",
        blank=True,
        help_text="Any other symptoms or observations"
    )
    
    class Meta:
        ordering = ["-assessment_date"]
        verbose_name = "Symptom Assessment"
        verbose_name_plural = "Symptom Assessments"
    
    def __str__(self) -> str:
        return f"Symptoms for {self.patient} on {self.assessment_date}"
    
    # ── Computed Properties ─────────────────────────────────────────────────
    
    @property
    def crab_count(self) -> int:
        """Count of active CRAB criteria."""
        return sum([
            self.crab_hypercalcemia,
            self.crab_renal,
            self.crab_anemia,
            self.crab_bone,
        ])
    
    @property
    def crab_string(self) -> str:
        """Active CRAB criteria as string (e.g., 'C-A-B')."""
        parts = []
        if self.crab_hypercalcemia:
            parts.append("C")
        if self.crab_renal:
            parts.append("R")
        if self.crab_anemia:
            parts.append("A")
        if self.crab_bone:
            parts.append("B")
        return "-".join(parts) if parts else "None"
    
    @property
    def max_neuropathy_grade(self) -> int:
        """Highest neuropathy grade across all types."""
        return max(
            self.neuropathy_sensory,
            self.neuropathy_motor,
            self.neuropathy_pain
        )
    
    @property
    def pain_severity(self) -> str:
        """Pain severity category based on NRS."""
        if self.bone_pain_nrs is None:
            return "Not assessed"
        if self.bone_pain_nrs == 0:
            return "None"
        if self.bone_pain_nrs <= 3:
            return "Mild"
        if self.bone_pain_nrs <= 6:
            return "Moderate"
        return "Severe"
    
    @property
    def fatigue_severity(self) -> str:
        """Fatigue severity based on FACIT-F score."""
        if self.fatigue_total is None:
            return "Not assessed"
        # Higher FACIT-F = less fatigue
        if self.fatigue_total >= 40:
            return "Minimal"
        if self.fatigue_total >= 30:
            return "Mild"
        if self.fatigue_total >= 20:
            return "Moderate"
        return "Severe"
    
    def calculate_facit_from_items(self) -> int | None:
        """Calculate FACIT-F total from individual items if available."""
        items = [
            self.fatigue_feel_tired,
            self.fatigue_feel_weak,
            self.fatigue_feel_listless,
            self.fatigue_trouble_starting,
            self.fatigue_trouble_finishing,
        ]
        if all(i is not None for i in items):
            # Reverse scoring: 4 - score for negative items
            # This is simplified; full FACIT has 13 items
            raw_sum = sum(items)
            # Scale to 0-52 range (simplified)
            return 52 - int((raw_sum / 20) * 52)
        return None
    
    def get_symptom_summary(self) -> dict:
        """Generate summary dict for API/display."""
        return {
            "date": self.assessment_date.isoformat() if self.assessment_date else None,
            "pain": {
                "nrs": self.bone_pain_nrs,
                "severity": self.pain_severity,
                "location": self.bone_pain_location,
            },
            "fatigue": {
                "facit_score": self.fatigue_total,
                "severity": self.fatigue_severity,
            },
            "neuropathy": {
                "sensory_grade": self.neuropathy_sensory,
                "motor_grade": self.neuropathy_motor,
                "pain_grade": self.neuropathy_pain,
                "max_grade": self.max_neuropathy_grade,
            },
            "crab": {
                "count": self.crab_count,
                "string": self.crab_string,
                "hypercalcemia": self.crab_hypercalcemia,
                "renal": self.crab_renal,
                "anemia": self.crab_anemia,
                "bone": self.crab_bone,
            },
            "infection": {
                "present": self.infection_present,
                "type": self.infection_type,
                "hospitalized": self.infection_hospitalized,
            },
            "ecog": self.ecog_status,
        }

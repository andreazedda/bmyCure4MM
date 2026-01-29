"""
Forms for SymptomAssessment preformatted input.

Provides user-friendly forms with validation and bilingual help text.
"""

from __future__ import annotations

from django import forms
from django.core.validators import MinValueValidator, MaxValueValidator

from mmportal.forms_mixins import BootstrapValidationMixin
from .models_symptoms import SymptomAssessment


class SymptomAssessmentForm(BootstrapValidationMixin, forms.ModelForm):
    """
    Comprehensive symptom assessment form with clinical scales.
    
    Uses validated instruments:
    - NRS for pain
    - FACIT-F for fatigue
    - CTCAE for neuropathy
    - CRAB criteria checkboxes
    """
    
    class Meta:
        model = SymptomAssessment
        exclude = ["patient", "assessment_date"]
        
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # ── Pain Section ────────────────────────────────────────────────────
        self.fields["bone_pain_nrs"].widget = forms.NumberInput(attrs={
            "type": "range",
            "min": "0",
            "max": "10",
            "step": "1",
            "class": "form-range",
            "oninput": "this.nextElementSibling.value = this.value",
        })
        self.fields["bone_pain_nrs"].help_text = (
            '<span class="t-en">Numeric Rating Scale: 0 = No pain, 10 = Worst pain imaginable</span>'
            '<span class="t-it">Scala Numerica: 0 = Nessun dolore, 10 = Peggior dolore immaginabile</span>'
        )
        
        self.fields["bone_pain_location"].widget = forms.TextInput(attrs={
            "class": "form-control",
            "placeholder": "e.g., lumbar spine, ribs, pelvis / es. colonna lombare, costole, bacino",
        })
        
        # ── Fatigue Section ─────────────────────────────────────────────────
        self.fields["fatigue_total"].widget = forms.NumberInput(attrs={
            "class": "form-control",
            "min": "0",
            "max": "52",
            "placeholder": "0-52",
        })
        self.fields["fatigue_total"].help_text = (
            '<span class="t-en">FACIT-Fatigue total score (0-52). Higher = LESS fatigue.</span>'
            '<span class="t-it">Punteggio totale FACIT-Fatigue (0-52). Più alto = MENO fatica.</span>'
        )
        
        # Quick fatigue questions
        for field_name in ["fatigue_feel_tired", "fatigue_feel_weak", "fatigue_feel_listless", 
                          "fatigue_trouble_starting", "fatigue_trouble_finishing"]:
            self.fields[field_name].widget = forms.Select(attrs={"class": "form-select"})
        
        # ── Neuropathy Section ──────────────────────────────────────────────
        for field_name in ["neuropathy_sensory", "neuropathy_motor", "neuropathy_pain"]:
            self.fields[field_name].widget = forms.Select(attrs={"class": "form-select"})
        
        self.fields["neuropathy_sensory"].help_text = (
            '<span class="t-en">Numbness, tingling, paresthesia</span>'
            '<span class="t-it">Intorpidimento, formicolio, parestesia</span>'
        )
        self.fields["neuropathy_motor"].help_text = (
            '<span class="t-en">Weakness, difficulty with fine motor skills</span>'
            '<span class="t-it">Debolezza, difficoltà con movimenti fini</span>'
        )
        self.fields["neuropathy_pain"].help_text = (
            '<span class="t-en">Burning, shooting, electric-shock pain</span>'
            '<span class="t-it">Dolore bruciante, lancinante, a scossa elettrica</span>'
        )
        
        # ── CRAB Section ────────────────────────────────────────────────────
        crab_fields = [
            ("crab_hypercalcemia", "Hypercalcemia (C) / Ipercalcemia", 
             "Calcium >11 mg/dL / Calcio >11 mg/dL"),
            ("crab_renal", "Renal Insufficiency (R) / Insufficienza Renale",
             "Creatinine >2 mg/dL / Creatinina >2 mg/dL"),
            ("crab_anemia", "Anemia (A)",
             "Hemoglobin <10 g/dL / Emoglobina <10 g/dL"),
            ("crab_bone", "Bone Lesions (B) / Lesioni Ossee",
             "Lytic lesions or fractures / Lesioni litiche o fratture"),
        ]
        
        for field_name, label, help_text in crab_fields:
            self.fields[field_name].widget = forms.CheckboxInput(attrs={
                "class": "form-check-input",
            })
            self.fields[field_name].label = label
            self.fields[field_name].help_text = help_text
        
        # CRAB symptom text fields
        for field_name in ["crab_hypercalcemia_symptoms", "crab_renal_symptoms", 
                          "crab_anemia_symptoms", "crab_bone_symptoms"]:
            self.fields[field_name].widget = forms.TextInput(attrs={
                "class": "form-control form-control-sm",
                "placeholder": "Specific symptoms / Sintomi specifici",
            })
        
        # ── Infection Section ───────────────────────────────────────────────
        self.fields["infection_present"].widget = forms.CheckboxInput(attrs={
            "class": "form-check-input",
            "onchange": "toggleInfectionFields(this.checked)",
        })
        
        self.fields["infection_type"].widget = forms.Select(attrs={
            "class": "form-select",
        })
        
        self.fields["infection_organism"].widget = forms.TextInput(attrs={
            "class": "form-control",
            "placeholder": "e.g., E. coli, S. pneumoniae",
        })
        
        self.fields["infection_hospitalized"].widget = forms.CheckboxInput(attrs={
            "class": "form-check-input",
        })
        
        # ── Performance Status ──────────────────────────────────────────────
        self.fields["ecog_status"].widget = forms.Select(attrs={
            "class": "form-select",
        })
        self.fields["ecog_status"].help_text = (
            '<span class="t-en">ECOG Performance Status (0-4)</span>'
            '<span class="t-it">Status di Performance ECOG (0-4)</span>'
        )
        
        # ── Additional Symptoms ─────────────────────────────────────────────
        for field_name in ["appetite_loss", "night_sweats", "bruising_bleeding"]:
            self.fields[field_name].widget = forms.CheckboxInput(attrs={
                "class": "form-check-input",
            })
        
        self.fields["weight_loss_pct"].widget = forms.NumberInput(attrs={
            "class": "form-control",
            "step": "0.1",
            "min": "0",
            "max": "100",
            "placeholder": "%",
        })
        self.fields["weight_loss_pct"].help_text = (
            '<span class="t-en">Percentage weight loss in last 3 months</span>'
            '<span class="t-it">Percentuale perdita peso ultimi 3 mesi</span>'
        )
        
        # ── Notes ───────────────────────────────────────────────────────────
        self.fields["additional_notes"].widget = forms.Textarea(attrs={
            "class": "form-control",
            "rows": "3",
            "placeholder": "Any other symptoms or observations / Altri sintomi o osservazioni",
        })


class QuickSymptomForm(BootstrapValidationMixin, forms.ModelForm):
    """
    Simplified symptom form for quick assessments.
    Captures only essential data points.
    """
    
    class Meta:
        model = SymptomAssessment
        fields = [
            "bone_pain_nrs",
            "fatigue_total",
            "neuropathy_sensory",
            "ecog_status",
            "crab_hypercalcemia",
            "crab_renal",
            "crab_anemia",
            "crab_bone",
            "infection_present",
        ]
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Pain slider
        self.fields["bone_pain_nrs"].widget = forms.NumberInput(attrs={
            "type": "range",
            "min": "0",
            "max": "10",
            "class": "form-range",
        })
        
        # Fatigue
        self.fields["fatigue_total"].widget = forms.NumberInput(attrs={
            "class": "form-control",
            "min": "0",
            "max": "52",
        })
        
        # Neuropathy
        self.fields["neuropathy_sensory"].widget = forms.Select(attrs={
            "class": "form-select",
        })
        
        # ECOG
        self.fields["ecog_status"].widget = forms.Select(attrs={
            "class": "form-select",
        })
        
        # CRAB checkboxes
        for field_name in ["crab_hypercalcemia", "crab_renal", "crab_anemia", 
                          "crab_bone", "infection_present"]:
            self.fields[field_name].widget = forms.CheckboxInput(attrs={
                "class": "form-check-input",
            })

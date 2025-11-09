from __future__ import annotations

from datetime import date

from django import forms
from django.core.exceptions import ValidationError

from mmportal.forms_mixins import BootstrapValidationMixin

from . import models


class PatientForm(BootstrapValidationMixin, forms.ModelForm):
    """Form for creating and editing patients."""

    mrn = forms.CharField(
        label="MRN (Medical Record Number)",
        max_length=50,
        help_text=(
            "Unique identifier. Examples: MM-2024-001, P12345, 9876543. "
            "Will be automatically converted to uppercase."
        ),
        widget=forms.TextInput(attrs={
            "placeholder": "e.g., MM-2024-001",
            "class": "form-control"
        })
    )
    
    first_name = forms.CharField(
        label="First Name / Nome",
        max_length=100,
        help_text="Patient's given name. Example: John, Maria, Ahmed",
        widget=forms.TextInput(attrs={
            "placeholder": "e.g., John",
            "class": "form-control"
        })
    )
    
    last_name = forms.CharField(
        label="Last Name / Cognome",
        max_length=100,
        help_text="Patient's family name. Example: Smith, Rossi, García",
        widget=forms.TextInput(attrs={
            "placeholder": "e.g., Smith",
            "class": "form-control"
        })
    )

    birth_date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date", "max": date.today().isoformat()}),
        help_text=(
            "Patient's date of birth (YYYY-MM-DD). "
            "Used to calculate current age. Must be in the past."
        )
    )
    
    diagnosis_date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date", "max": date.today().isoformat()}),
        help_text=(
            "Date when Multiple Myeloma was diagnosed (YYYY-MM-DD). "
            "Must be after birth date and not in the future."
        )
    )
    
    sex = forms.ChoiceField(
        label="Sex / Sesso",
        choices=[("", "Select sex / Seleziona sesso")] + list(models.Patient.SEX_CHOICES),
        help_text="Biological sex (M=Male/Maschio, F=Female/Femmina)",
        widget=forms.Select(attrs={"class": "form-select"})
    )
    
    notes = forms.CharField(
        label="Clinical Notes / Note Cliniche",
        required=False,
        help_text=(
            "Optional additional information. "
            "Example: Comorbidities, allergies, special considerations"
        ),
        widget=forms.Textarea(attrs={
            "rows": "3",
            "placeholder": "e.g., Diabetes mellitus type 2, allergic to penicillin...",
            "class": "form-control"
        })
    )

    class Meta:
        model = models.Patient
        fields = [
            "mrn",
            "first_name", 
            "last_name",
            "birth_date",
            "sex",
            "diagnosis_date",
            "notes",
        ]
        
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Bootstrap styling already handled by field definitions

    def clean_mrn(self):
        mrn = self.cleaned_data.get("mrn")
        if mrn:
            mrn = mrn.strip().upper()
            # Check for uniqueness if this is a new patient or MRN changed
            if self.instance.pk is None or self.instance.mrn != mrn:
                if models.Patient.objects.filter(mrn=mrn).exists():
                    raise forms.ValidationError(
                        "⚠️ A patient with this MRN already exists. "
                        "Each patient must have a unique Medical Record Number."
                    )
        return mrn
    
    def clean_birth_date(self):
        birth_date = self.cleaned_data.get("birth_date")
        if birth_date:
            if birth_date > date.today():
                raise ValidationError("⚠️ Birth date cannot be in the future.")
            # Check age is reasonable (e.g., not more than 120 years old)
            age_years = (date.today() - birth_date).days / 365.25
            if age_years > 120:
                raise ValidationError(
                    f"⚠️ Birth date results in age {int(age_years)} years. Please verify."
                )
            if age_years < 18:
                raise ValidationError(
                    "⚠️ Patient must be at least 18 years old (Multiple Myeloma is rare in children)."
                )
        return birth_date
    
    def clean_diagnosis_date(self):
        diagnosis_date = self.cleaned_data.get("diagnosis_date")
        if diagnosis_date:
            if diagnosis_date > date.today():
                raise ValidationError("⚠️ Diagnosis date cannot be in the future.")
        return diagnosis_date
    
    def clean(self):
        cleaned = super().clean()
        birth_date = cleaned.get("birth_date")
        diagnosis_date = cleaned.get("diagnosis_date")
        
        if birth_date and diagnosis_date:
            if diagnosis_date < birth_date:
                raise ValidationError(
                    "⚠️ Diagnosis date must be after birth date."
                )
            
            # Check diagnosis age is reasonable
            age_at_diagnosis = (diagnosis_date - birth_date).days / 365.25
            if age_at_diagnosis < 18:
                self.add_error("diagnosis_date", 
                    f"⚠️ Patient was {int(age_at_diagnosis)} years old at diagnosis. "
                    "Multiple Myeloma is extremely rare in patients under 18."
                )
        
        return cleaned


class AssessmentForm(BootstrapValidationMixin, forms.ModelForm):
    """Quick entry form for laboratory assessment with normal ranges and validation."""

    date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date", "max": date.today().isoformat()}),
        help_text="Assessment date (cannot be in future)"
    )
    
    m_protein_g_dl = forms.DecimalField(
        label="M-Protein (g/dL)",
        required=False,
        decimal_places=2,
        min_value=0,
        max_value=20,
        help_text="Normal: 0 g/dL. Lower is better. Range: 0-20 g/dL",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 2.5",
            "step": "0.01",
            "class": "form-control"
        })
    )
    
    flc_ratio = forms.DecimalField(
        label="FLC Ratio (κ/λ)",
        required=False,
        decimal_places=3,
        min_value=0.001,
        max_value=1000,
        help_text="Normal range: 0.26-1.65. Outside this range indicates clonal disease.",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 0.85",
            "step": "0.001",
            "class": "form-control"
        })
    )
    
    hemoglobin_g_dl = forms.DecimalField(
        label="Hemoglobin (g/dL)",
        required=False,
        decimal_places=1,
        min_value=0,
        max_value=25,
        help_text="Normal: M:13-17, F:12-15 g/dL. Low indicates anemia.",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 12.5",
            "step": "0.1",
            "class": "form-control"
        })
    )
    
    calcium_mg_dl = forms.DecimalField(
        label="Calcium (mg/dL)",
        required=False,
        decimal_places=1,
        min_value=0,
        max_value=20,
        help_text="Normal: 8.5-10.5 mg/dL. High (>11) indicates hypercalcemia.",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 9.2",
            "step": "0.1",
            "class": "form-control"
        })
    )
    
    creatinine_mg_dl = forms.DecimalField(
        label="Creatinine (mg/dL)",
        required=False,
        decimal_places=2,
        min_value=0,
        max_value=15,
        help_text="Normal: 0.6-1.2 mg/dL. High indicates kidney damage.",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 0.9",
            "step": "0.01",
            "class": "form-control"
        })
    )
    
    beta2m_mg_l = forms.DecimalField(
        label="β2-Microglobulin (mg/L)",
        required=False,
        decimal_places=2,
        min_value=0,
        max_value=50,
        help_text="Normal: <2.5 mg/L. Used for R-ISS staging. 2.5-5.5: R-ISS I/II, >5.5: R-ISS III",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 3.2",
            "step": "0.01",
            "class": "form-control"
        })
    )
    
    albumin_g_dl = forms.DecimalField(
        label="Albumin (g/dL)",
        required=False,
        decimal_places=1,
        min_value=0,
        max_value=10,
        help_text="Normal: 3.5-5.5 g/dL. Low indicates malnutrition or disease.",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 4.0",
            "step": "0.1",
            "class": "form-control"
        })
    )
    
    ldH_u_l = forms.DecimalField(
        label="LDH (U/L)",
        required=False,
        decimal_places=0,
        min_value=0,
        max_value=10000,
        help_text="Normal: 120-240 U/L. High indicates tissue damage or disease progression.",
        widget=forms.NumberInput(attrs={
            "placeholder": "e.g., 180",
            "step": "1",
            "class": "form-control"
        })
    )
    
    r_iss = forms.ChoiceField(
        label="R-ISS Stage",
        required=False,
        choices=[("", "Not assessed")] + list(models.Assessment.R_ISS_CHOICES),
        help_text="Revised International Staging System: I (best), II (intermediate), III (worst)",
        widget=forms.Select(attrs={"class": "form-select"})
    )
    
    response = forms.ChoiceField(
        label="Treatment Response",
        required=False,
        choices=[("", "Not assessed")] + list(models.Assessment.RESPONSE_CHOICES),
        help_text="CR=Complete, VGPR=Very Good Partial, PR=Partial, SD=Stable, PD=Progressive",
        widget=forms.Select(attrs={"class": "form-select"})
    )

    class Meta:
        model = models.Assessment
        fields = [
            "date",
            "m_protein_g_dl",
            "kappa_mg_l",
            "lambda_mg_l",
            "flc_ratio",
            "hemoglobin_g_dl",
            "calcium_mg_dl",
            "creatinine_mg_dl",
            "beta2m_mg_l",
            "albumin_g_dl",
            "ldH_u_l",
            "r_iss",
            "response",
            "notes",
        ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name, field in self.fields.items():
            if name not in ['m_protein_g_dl', 'flc_ratio', 'hemoglobin_g_dl', 
                           'calcium_mg_dl', 'creatinine_mg_dl', 'beta2m_mg_l',
                           'albumin_g_dl', 'ldH_u_l', 'r_iss', 'response', 'date']:
                widget = field.widget
                if isinstance(widget, forms.Select):
                    widget.attrs.setdefault("class", "form-select")
                elif isinstance(widget, forms.CheckboxInput):
                    widget.attrs.setdefault("class", "form-check-input")
                else:
                    widget.attrs.setdefault("class", "form-control")
    
    def clean_date(self):
        assessment_date = self.cleaned_data.get("date")
        if assessment_date and assessment_date > date.today():
            raise ValidationError("⚠️ Assessment date cannot be in the future.")
        return assessment_date

    def clean(self):
        cleaned = super().clean()
        
        # Validate M-Protein range
        m_protein = cleaned.get("m_protein_g_dl")
        if m_protein is not None:
            if m_protein > 10:
                self.add_error("m_protein_g_dl", 
                    "⚠️ Very high M-Protein (>10 g/dL). Please verify measurement.")
        
        # Validate FLC ratio
        flc_ratio = cleaned.get("flc_ratio")
        if flc_ratio is not None:
            if flc_ratio < 0.26 or flc_ratio > 1.65:
                # This is expected in MM, but flag if extremely abnormal
                if flc_ratio < 0.01 or flc_ratio > 100:
                    self.add_error("flc_ratio",
                        f"⚠️ Extremely abnormal FLC ratio ({flc_ratio}). Normal: 0.26-1.65. Please verify.")
        
        # Validate calcium
        calcium = cleaned.get("calcium_mg_dl")
        if calcium is not None and calcium > 14:
            self.add_error("calcium_mg_dl",
                "⚠️ Severe hypercalcemia (>14 mg/dL). Requires immediate intervention!")
        
        # Validate creatinine
        creatinine = cleaned.get("creatinine_mg_dl")
        if creatinine is not None and creatinine > 10:
            self.add_error("creatinine_mg_dl",
                "⚠️ Severe renal failure (>10 mg/dL). Patient may require dialysis.")
        
        return cleaned


class PatientTherapyForm(BootstrapValidationMixin, forms.ModelForm):
    """Manage therapy course entries."""

    start_date = forms.DateField(widget=forms.DateInput(attrs={"type": "date"}))
    end_date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date"}),
        required=False,
    )

    class Meta:
        model = models.PatientTherapy
        fields = [
            "regimen",
            "start_date",
            "end_date",
            "outcome",
            "adverse_events",
        ]

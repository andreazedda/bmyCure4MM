from __future__ import annotations

from django import forms
from django.core.exceptions import ValidationError

from mmportal.forms_mixins import BootstrapValidationMixin

from clinic.models import Assessment, Regimen

from . import models
from .presets import PRESETS


class RegimenForm(BootstrapValidationMixin, forms.ModelForm):
    """Editor form for maintaining regimens."""

    INTENT_CHOICES = [
        ("", "Select intent..."),
        ("curative", "Curative / Induction"),
        ("maintenance", "Maintenance"),
        ("consolidation", "Consolidation"),
        ("salvage", "Salvage / Relapsed"),
        ("palliative", "Palliative / Supportive"),
    ]

    intent = forms.ChoiceField(
        label="Treatment intent",
        required=False,
        choices=INTENT_CHOICES,
        widget=forms.Select(attrs={"class": "form-select"}),
    )
    components = forms.CharField(
        label="Components",
        widget=forms.Textarea(attrs={"class": "form-control", "rows": 3}),
    )
    notes = forms.CharField(
        required=False,
        widget=forms.Textarea(attrs={"class": "form-control", "rows": 3}),
    )

    class Meta:
        model = Regimen
        fields = ["name", "line", "components", "intent", "notes"]
        widgets = {
            "name": forms.TextInput(attrs={"class": "form-control"}),
            "line": forms.TextInput(attrs={"class": "form-control"}),
        }


class ScenarioForm(BootstrapValidationMixin, forms.ModelForm):
    """Editor form for scenarios."""

    recommended_regimens = forms.ModelMultipleChoiceField(
        queryset=Regimen.objects.order_by("name"),
        required=False,
        widget=forms.SelectMultiple(attrs={"class": "form-select", "size": "8"}),
    )

    class Meta:
        model = models.Scenario
        fields = [
            "title",
            "clinical_stage",
            "risk_stratification",
            "summary",
            "guideline_notes",
            "recommended_regimens",
            "expected_response",
            "active",
        ]
        widgets = {
            "title": forms.TextInput(attrs={"class": "form-control"}),
            "clinical_stage": forms.Select(attrs={"class": "form-select"}),
            "risk_stratification": forms.TextInput(attrs={"class": "form-control"}),
            "summary": forms.Textarea(attrs={"class": "form-control", "rows": 4}),
            "guideline_notes": forms.Textarea(attrs={"class": "form-control", "rows": 4}),
            "active": forms.CheckboxInput(attrs={"class": "form-check-input"}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        expected_field = self.fields["expected_response"]
        expected_field.choices = [("", "Select response...")] + list(Assessment.RESPONSE_CHOICES)
        expected_field.widget.attrs.setdefault("class", "form-select")


class SimulationParameterForm(BootstrapValidationMixin, forms.Form):
    """Guardrailed three-step parameter capture for simulations."""

    PRESET_CHOICES = [(key, value["label"]) for key, value in PRESETS.items()]

    preset = forms.ChoiceField(choices=PRESET_CHOICES, required=True, label="Treatment preset")

    creatinine_clearance = forms.FloatField(
        min_value=5,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "1",
                "min": "5",
                "max": "150",
                "inputmode": "decimal",
            }
        ),
        label="Creatinine clearance (ml/min)",
        help_text="Renal function estimate driving lenalidomide dose adjustments (thresholds at 60 and 30 ml/min).",
    )
    neuropathy_grade = forms.ChoiceField(
        choices=[(str(i), f"Grade {i}") for i in range(0, 4)],
        label="Peripheral neuropathy (grade)",
        widget=forms.Select(attrs={"class": "form-select"}),
        help_text="CTCAE sensory neuropathy grade (0–3). Grades ≥2 restrict bortezomib to ≤1.0 mg/m².",
    )
    anc = forms.FloatField(
        min_value=0.1,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0.1",
                "max": "10",
                "inputmode": "decimal",
            }
        ),
        label="Absolute neutrophil count (×10⁹/L)",
        help_text="ANC <1.0 ×10⁹/L indicates clinically significant neutropenia and blocks simulation.",
    )
    platelets = forms.FloatField(
        min_value=10,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "5",
                "min": "10",
                "max": "400",
                "inputmode": "decimal",
            }
        ),
        label="Platelets (×10⁹/L)",
        help_text="Platelets <75 ×10⁹/L represent high bleeding risk and block simulation.",
    )
    pregnancy = forms.BooleanField(
        required=False,
        widget=forms.CheckboxInput(attrs={"class": "form-check-input"}),
        label="Pregnancy flagged",
        help_text="Pregnancy is a contraindication for IMiDs; simulation runs are blocked when flagged.",
    )

    baseline_tumor_cells = forms.FloatField(
        min_value=1e6,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "1000000",
                "min": "1000000",
                "max": "1000000000000",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Initial tumor cells",
        help_text="Estimated malignant plasma cell burden (cells). Newly diagnosed cases often present with 10⁹–10¹² cells.",
    )
    baseline_healthy_cells = forms.FloatField(
        min_value=1e8,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "100000000",
                "min": "100000000",
                "max": "10000000000000",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Initial healthy plasma cells",
        help_text="Approximate pool of normal plasma cells (cells). Typical range ≈10¹¹–10¹².",
    )
    lenalidomide_dose = forms.FloatField(
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0",
                "max": "50",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Lenalidomide daily dose (mg)",
        help_text="Standard induction dose: 25 mg/day (21 days on, 7 off). Preset tweaks limited to ±20%.",
    )
    bortezomib_dose = forms.FloatField(
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0",
                "max": "2",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Bortezomib weekly dose (mg/m²)",
        help_text="Typical SC dosing 1.3 mg/m² on specified days; slider allows controlled ±20% adjustment.",
    )
    daratumumab_dose = forms.FloatField(
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0",
                "max": "20",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Daratumumab dose (mg/kg)",
        help_text="Anti-CD38 monoclonal antibody loading dose 16 mg/kg; slider extends within investigational bounds.",
    )
    time_horizon = forms.FloatField(
        min_value=7.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "1",
                "min": "7",
                "max": "365",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Simulation horizon (days)",
        help_text="Length of the virtual treatment window (days). Limits ensure numerical stability and clinical realism.",
    )
    tumor_growth_rate = forms.FloatField(
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.001",
                "min": "0",
                "max": "0.1",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Tumor growth rate (day⁻¹)",
        help_text="Logistic growth rate of malignant plasma cells. Values >0.05 day⁻¹ are rare and capped.",
    )
    healthy_growth_rate = forms.FloatField(
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.001",
                "min": "0",
                "max": "0.05",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Healthy growth rate (day⁻¹)",
        help_text="Marrow recovery kinetics for healthy plasma cells. Typical range 0.01–0.02 day⁻¹.",
    )
    interaction_strength = forms.FloatField(
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.01",
                "min": "0",
                "max": "0.2",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Drug interaction strength",
        help_text="Synergy/toxicity coupling between agents (dimensionless). Values >0.2 are unsupported.",
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.warnings: list[str] = []
        self.preset_key = self.data.get("preset") or self.initial.get("preset") or self.PRESET_CHOICES[0][0]
        if self.preset_key not in PRESETS:
            self.preset_key = self.PRESET_CHOICES[0][0]
        self.preset_config = PRESETS[self.preset_key]
        defaults = self.preset_config["default_params"]

        slider_fields = [
            "baseline_tumor_cells",
            "baseline_healthy_cells",
            "lenalidomide_dose",
            "bortezomib_dose",
            "daratumumab_dose",
            "time_horizon",
            "tumor_growth_rate",
            "healthy_growth_rate",
            "interaction_strength",
        ]

        for field in slider_fields:
            if field in defaults and field not in self.data:
                self.fields[field].initial = defaults[field]

        if "creatinine_clearance" not in self.data:
            self.fields["creatinine_clearance"].initial = 90
        if "neuropathy_grade" not in self.data:
            self.fields["neuropathy_grade"].initial = "0"
        if "anc" not in self.data:
            self.fields["anc"].initial = 2.0
        if "platelets" not in self.data:
            self.fields["platelets"].initial = 150

        self.slider_bounds = self._compute_slider_bounds(defaults, slider_fields)

    def _compute_slider_bounds(self, defaults: dict[str, float], slider_fields: list[str]) -> dict[str, dict[str, float]]:
        bounds = {}
        pct = float(self.preset_config.get("bounds_pct", 20)) / 100.0
        for field in slider_fields:
            base_val = float(self.fields[field].initial or defaults.get(field, 0.0))
            bounds[field] = {
                "base": base_val,
                "min": base_val * (1.0 - pct),
                "max": base_val * (1.0 + pct),
            }
        return bounds

    def get_slider_payload(self) -> dict[str, dict[str, float]]:
        return self.slider_bounds

    def get_presets_payload(self) -> dict[str, dict[str, object]]:
        payload = {}
        for key, config in PRESETS.items():
            defaults = config["default_params"]
            pct = float(config.get("bounds_pct", 20)) / 100.0
            bounds = {}
            for field, value in defaults.items():
                val = float(value)
                bounds[field] = {
                    "base": val,
                    "min": val * (1.0 - pct),
                    "max": val * (1.0 + pct),
                }
            payload[key] = {
                "label": config["label"],
                "bounds": bounds,
                "schedule": config.get("schedule", {}),
            }
        return payload

    def _clean_numeric(self, cleaned: dict[str, float], key: str) -> float:
        value = cleaned.get(key)
        if value is None:
            initial = self.fields[key].initial
            return float(initial or 0.0)
        return float(value)

    def clean(self):
        cleaned = super().clean()
        cleaned["preset"] = self.preset_key
        cleaned["schedule"] = self.preset_config.get("schedule", {})

        len_dose = self._clean_numeric(cleaned, "lenalidomide_dose")
        bor_dose = self._clean_numeric(cleaned, "bortezomib_dose")
        dara_dose = self._clean_numeric(cleaned, "daratumumab_dose")
        time_horizon = self._clean_numeric(cleaned, "time_horizon")
        tumor_growth = self._clean_numeric(cleaned, "tumor_growth_rate")
        healthy_growth = self._clean_numeric(cleaned, "healthy_growth_rate")
        baseline_tumor = self._clean_numeric(cleaned, "baseline_tumor_cells")
        baseline_healthy = self._clean_numeric(cleaned, "baseline_healthy_cells")
        interaction_strength = self._clean_numeric(cleaned, "interaction_strength")
        creatinine = cleaned.get("creatinine_clearance")
        neuropathy = int(cleaned.get("neuropathy_grade") or 0)
        anc = cleaned.get("anc")
        platelets = cleaned.get("platelets")
        pregnancy = cleaned.get("pregnancy")

        errors: dict[str, ValidationError] = {}

        if pregnancy:
            raise ValidationError("Simulation blocked: pregnancy flagged.")

        if baseline_tumor > 1_000_000_000_000:
            errors["baseline_tumor_cells"] = ValidationError(
                "Baseline tumor burden exceeds modeled capacity (≤ 1×10¹² cells)."
            )
        if baseline_healthy > 10_000_000_000_000:
            errors["baseline_healthy_cells"] = ValidationError(
                "Healthy plasma cell pool must be ≤ 1×10¹³ cells."
            )
        if len_dose > 50:
            errors["lenalidomide_dose"] = ValidationError("Lenalidomide dose must be ≤ 50 mg/day.")
        if bor_dose > 2:
            errors["bortezomib_dose"] = ValidationError("Bortezomib dose must be ≤ 2 mg/m² per week.")
        if dara_dose > 20:
            errors["daratumumab_dose"] = ValidationError("Daratumumab dose must be ≤ 20 mg/kg.")
        if time_horizon > 365:
            errors["time_horizon"] = ValidationError("Simulation horizon is limited to 365 days.")
        if tumor_growth > 0.1:
            errors["tumor_growth_rate"] = ValidationError("Tumor growth rate must be ≤ 0.10 day⁻¹.")
        if healthy_growth > 0.05:
            errors["healthy_growth_rate"] = ValidationError("Healthy growth rate must be ≤ 0.05 day⁻¹.")
        if interaction_strength > 0.2:
            errors["interaction_strength"] = ValidationError("Interaction strength must be ≤ 0.2.")

        if creatinine is not None:
            if creatinine < 30 and len_dose > 10:
                errors["lenalidomide_dose"] = ValidationError(
                    "Creatinine clearance <30 ml/min requires lenalidomide dose ≤10 mg."
                )
            elif creatinine < 60 and len_dose > 15:
                errors["lenalidomide_dose"] = ValidationError(
                    "Creatinine clearance <60 ml/min requires lenalidomide dose ≤15 mg."
                )

        if neuropathy >= 2 and bor_dose > 1.0:
            errors["bortezomib_dose"] = ValidationError(
                "Peripheral neuropathy grade ≥2 mandates bortezomib dose ≤1.0 mg/m²."
            )

        if anc is not None and anc < 1.0:
            raise ValidationError("Simulation blocked: absolute neutrophil count <1.0 ×10⁹/L (myelosuppression).")
        if platelets is not None and platelets < 75:
            raise ValidationError("Simulation blocked: platelets <75 ×10⁹/L (thrombocytopenia).")

        if len_dose > 40 and bor_dose > 1.5:
            raise ValidationError(
                "Combined high doses of lenalidomide (>40 mg) and bortezomib (>1.5 mg/m²) may exceed safe toxicity bounds."
            )

        if interaction_strength > 0.15:
            self.warnings.append("Interaction strength above 0.15 implies potent synergy—monitor toxicity closely.")
        if 35 < len_dose <= 40:
            self.warnings.append("Lenalidomide dose is near the upper investigational range (>35 mg).")
        if 1.3 < bor_dose <= 1.5:
            self.warnings.append("Bortezomib dose exceeds standard 1.3 mg/m²—evaluate neuropathy risk.")
        if 16 < dara_dose <= 20:
            self.warnings.append("Daratumumab dose beyond typical loading (16 mg/kg) may increase immunosuppression.")
        if time_horizon >= 300:
            self.warnings.append("Long simulation horizon (>300 days) may magnify numerical stiffness; consider shorter intervals.")

        if errors:
            raise ValidationError(errors)

        return cleaned
class SimulationAttemptForm(BootstrapValidationMixin, forms.ModelForm):
    """Form used for clinician responses within a scenario."""

    selected_regimen = forms.ModelChoiceField(
        label="Chosen regimen",
        queryset=Regimen.objects.all(),
        required=False,
        help_text="Select the primary regimen you would start for this case.",
        widget=forms.Select(attrs={"class": "form-select"}),
    )
    predicted_response = forms.ChoiceField(
        label="Expected response",
        choices=[("", "Select response...")] + list(Assessment.RESPONSE_CHOICES),
        required=False,
        widget=forms.Select(attrs={"class": "form-select"}),
    )
    confidence = forms.IntegerField(
        min_value=0,
        max_value=100,
        initial=70,
        help_text="How confident are you in this plan? (%)",
        widget=forms.NumberInput(attrs={"class": "form-control", "min": "0", "max": "100"}),
    )
    notes = forms.CharField(
        widget=forms.Textarea(attrs={"rows": 4, "class": "form-control"}),
        required=False,
        label="Clinical reasoning",
    )

    class Meta:
        model = models.SimulationAttempt
        fields = ["selected_regimen", "predicted_response", "confidence", "notes"]

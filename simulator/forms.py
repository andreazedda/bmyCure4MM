from __future__ import annotations

from django import forms
from django.core.exceptions import ValidationError

from mmportal.forms_mixins import BootstrapValidationMixin

from clinic.models import Assessment, Regimen

from . import models
from .presets import PRESETS

SIMULATION_FORM_HELP_TEXT_EN = {
    "preset": "Clinical preset with guardrailed defaults for common regimens.",
    "creatinine_clearance": "Renal function estimate driving lenalidomide adjustments (thresholds at 60 and 30 ml/min).",
    "neuropathy_grade": "CTCAE sensory neuropathy grade (0–3). Grades ≥2 restrict bortezomib to ≤1.0 mg/m².",
    "anc": "ANC <1.0 ×10⁹/L blocks the simulation because of neutropenia.",
    "platelets": "Platelets <75 ×10⁹/L block the simulation (bleeding risk).",
    "pregnancy": "Pregnancy is a contraindication for IMiDs; simulations are blocked when flagged.",
    "baseline_tumor_cells": "Estimated malignant plasma cell burden (cells).",
    "baseline_healthy_cells": "Approximate pool of normal plasma cells (cells).",
    "lenalidomide_dose": "Standard induction dose: 25 mg/day (21 days on, 7 off). Adjustments are limited to preset bounds.",
    "bortezomib_dose": "Typical SC dosing 1.3 mg/m² on scheduled days; stay within ±20% to avoid toxicity.",
    "daratumumab_dose": "Loading dose 16 mg/kg; higher exposure raises immunosuppression risk.",
    "time_horizon": "Length of the virtual treatment window (days). Longer horizons increase numerical stiffness.",
    "tumor_growth_rate": "Logistic growth rate of malignant plasma cells. Values >0.05 day⁻¹ are rare.",
    "healthy_growth_rate": "Marrow recovery kinetics for healthy plasma cells (≈0.01–0.02 day⁻¹).",
    "interaction_strength": "Synergy/toxicity coupling between agents (dimensionless).",
    "cohort_size": "Number of virtual subjects sampled for uncertainty bands.",
    "use_twin": "Enable the Patient Twin to auto-derive biology from the latest labs.",
    "seed": "Optional seed for reproducible virtual cohorts and solver noise.",
}

SIMULATION_FORM_HELP_TEXT_IT = {
    "preset": "Preset clinico con default protetti per i regimi più comuni.",
    "creatinine_clearance": "Stima della funzione renale che guida gli aggiustamenti di lenalidomide (soglie a 60 e 30 ml/min).",
    "neuropathy_grade": "Grado di neuropatia sensitiva CTCAE (0–3). Gradi ≥2 limitano bortezomib a ≤1.0 mg/m².",
    "anc": "ANC <1.0 ×10⁹/L indica neutropenia clinicamente significativa e blocca la simulazione.",
    "platelets": "Piastrine <75 ×10⁹/L rappresentano un alto rischio emorragico e bloccano la simulazione.",
    "pregnancy": "La gravidanza è una controindicazione per gli IMiD; le prove virtuali vengono bloccate quando selezionata.",
    "baseline_tumor_cells": "Stima del carico di cellule plasmatiche maligne (cellule).",
    "baseline_healthy_cells": "Pool approssimativo di cellule plasmatiche sane (cellule).",
    "lenalidomide_dose": "Dose standard: 25 mg/die (21 giorni su 28). Il preset limita l’aggiustamento controllato.",
    "bortezomib_dose": "Dose SC tipica 1.3 mg/m² nei giorni programmati; restare entro ±20% evita tossicità.",
    "daratumumab_dose": "Dose di carico 16 mg/kg; variazioni maggiori aumentano il rischio d’immunosoppressione.",
    "time_horizon": "Durata della finestra terapeutica virtuale (giorni). Orizzonti lunghi aumentano l’incertezza.",
    "tumor_growth_rate": "Velocità logistica di crescita del tumore. Valori >0.05 day⁻¹ sono rari.",
    "healthy_growth_rate": "Cinetica di recupero per le cellule sane. Tipicamente 0.01–0.02 day⁻¹.",
    "interaction_strength": "Ampiezza della sinergia/tossicità tra i farmaci (senza unità).",
    "cohort_size": "Numero di pazienti virtuali campionati per stimare le bande di incertezza.",
    "use_twin": "Attiva il Gemello Paziente per riempire automaticamente i parametri biologici dai laboratori disponibili.",
    "seed": "Seed opzionale per rendere ripetibili coorti virtuali e solver.",
}


class RegimenForm(BootstrapValidationMixin, forms.ModelForm):
    """
    Enhanced editor form for maintaining regimens with dosing validation.
    
    This form enforces meaningful data entry by:
    - Validating drug dosing within safe clinical ranges
    - Checking for contraindications and interactions
    - Providing educational warnings with clinical rationale
    - Suggesting standard dosing regimens
    
    Components field format: "Drug1: dose units, Drug2: dose units, ..."
    Examples:
    - "Lenalidomide: 25 mg days 1-21, Dexamethasone: 40 mg weekly"
    - "Bortezomib: 1.3 mg/m² days 1,4,8,11, Dexamethasone: 20 mg days 1,2,4,5,8,9,11,12"
    """

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
        widget=forms.Textarea(attrs={
            "class": "form-control",
            "rows": 5,
            "placeholder": "e.g., Lenalidomide: 25 mg days 1-21 every 28 days, Dexamethasone: 40 mg weekly"
        }),
        help_text=(
            "List drugs with dosing. Format: 'Drug: dose schedule'. "
            "Standard dosing examples: Lenalidomide 25mg, Bortezomib 1.3mg/m², "
            "Daratumumab 16mg/kg, Carfilzomib 20-56mg/m², Dexamethasone 20-40mg"
        )
    )
    notes = forms.CharField(
        required=False,
        widget=forms.Textarea(attrs={
            "class": "form-control",
            "rows": 3,
            "placeholder": "Additional notes, adjustments, or special considerations..."
        }),
    )

    class Meta:
        model = Regimen
        fields = ["name", "line", "components", "intent", "notes"]
        widgets = {
            "name": forms.TextInput(attrs={
                "class": "form-control",
                "placeholder": "e.g., VRd (Bortezomib-Lenalidomide-Dexamethasone)"
            }),
            "line": forms.TextInput(attrs={
                "class": "form-control",
                "placeholder": "e.g., First-line, Second-line"
            }),
        }
    
    # Standard dosing ranges (mg or mg/m²)
    DRUG_DOSING_GUIDELINES = {
        "lenalidomide": {
            "min": 5,
            "max": 25,
            "standard": 25,
            "unit": "mg",
            "warnings": {
                "max_exceeded": "⚠️ Lenalidomide >25mg/day significantly increases neutropenia and VTE risk. Standard dose is 25mg days 1-21.",
                "adjust_renal": "💡 Lenalidomide requires dose reduction for CrCl <60: CrCl 30-60 → 10mg, CrCl <30 → 15mg every other day.",
            }
        },
        "bortezomib": {
            "min": 0.7,
            "max": 1.3,
            "standard": 1.3,
            "unit": "mg/m²",
            "warnings": {
                "max_exceeded": "⚠️ Bortezomib >1.3mg/m² increases peripheral neuropathy risk. Standard is 1.3mg/m² twice weekly.",
                "neuropathy": "💡 For neuropathy grade ≥2: reduce to 1.0mg/m² or switch to weekly dosing.",
                "subcutaneous": "💡 Subcutaneous administration reduces neuropathy vs IV (preferred route).",
            }
        },
        "daratumumab": {
            "min": 8,
            "max": 16,
            "standard": 16,
            "unit": "mg/kg",
            "warnings": {
                "max_exceeded": "⚠️ Daratumumab >16mg/kg not recommended. Standard is 16mg/kg weekly ×8, then Q2W ×16, then Q4W.",
                "infusion_reactions": "💡 Premedicate with antihistamines, acetaminophen, corticosteroids. Monitor for infusion reactions (first dose).",
            }
        },
        "carfilzomib": {
            "min": 20,
            "max": 56,
            "standard": 27,  # 20/27 or 20/56 dosing
            "unit": "mg/m²",
            "warnings": {
                "max_exceeded": "⚠️ Carfilzomib >56mg/m² exceeds FDA-approved dosing. Standard: 20mg/m² cycle 1 day 1-2, then 27 or 56mg/m².",
                "cardiac": "⚠️ Monitor cardiac function. Carfilzomib increases cardiac events (HF, HTN). Baseline ECHO recommended.",
                "hydration": "💡 IV hydration (250-500mL pre/post) reduces renal toxicity.",
            }
        },
        "dexamethasone": {
            "min": 4,
            "max": 40,
            "standard": 40,
            "unit": "mg",
            "warnings": {
                "max_exceeded": "⚠️ Dexamethasone >40mg/week increases infection, hyperglycemia, psychiatric effects. Consider 20mg for elderly/frail.",
                "elderly": "💡 For age >75 or frail: reduce to 20mg weekly to minimize toxicity.",
                "monitoring": "💡 Monitor glucose, blood pressure, mood. Prophylaxis for PCP if prolonged high-dose use.",
            }
        },
        "pomalidomide": {
            "min": 2,
            "max": 4,
            "standard": 4,
            "unit": "mg",
            "warnings": {
                "max_exceeded": "⚠️ Pomalidomide >4mg/day not recommended. Standard: 4mg days 1-21 every 28 days.",
                "adjust_renal": "💡 No dose reduction for renal impairment, but monitor closely if CrCl <45.",
            }
        },
        "cyclophosphamide": {
            "min": 50,
            "max": 500,
            "standard": 300,
            "unit": "mg/m²",
            "warnings": {
                "max_exceeded": "⚠️ Cyclophosphamide >500mg/m²/week increases myelosuppression. Typical: 300mg/m² weekly.",
                "hydration": "💡 Ensure adequate hydration to prevent hemorrhagic cystitis. Consider mesna for high doses.",
            }
        },
        "melphalan": {
            "min": 5,
            "max": 200,
            "standard": 9,  # oral dosing
            "unit": "mg/m²",
            "warnings": {
                "high_dose": "⚠️ Melphalan >40mg/m² is high-dose (HDM) requiring stem cell support. Typical oral: 9mg/m² days 1-4.",
                "myelosuppression": "⚠️ Melphalan causes profound myelosuppression. Monitor CBC closely.",
            }
        },
    }
    
    def clean_components(self):
        """Validate drug components and dosing."""
        components = self.cleaned_data.get("components", "")
        if not components:
            raise forms.ValidationError(
                "💡 Please specify regimen components with dosing. "
                "Example: 'Lenalidomide: 25 mg days 1-21, Dexamethasone: 40 mg weekly'"
            )
        
        # Parse components for known drugs
        components_lower = components.lower()
        warnings = []
        
        for drug, guidelines in self.DRUG_DOSING_GUIDELINES.items():
            if drug in components_lower:
                # Extract numeric values
                import re
                # Match patterns like "25 mg", "1.3 mg/m2", "16 mg/kg"
                dose_pattern = r'(\d+(?:\.\d+)?)\s*(?:mg|mg/m²|mg/m2|mg/kg)'
                matches = re.findall(dose_pattern, components, re.IGNORECASE)
                
                if matches:
                    doses = [float(m) for m in matches]
                    max_dose = max(doses)
                    
                    if max_dose > guidelines["max"]:
                        warnings.append(guidelines["warnings"]["max_exceeded"])
                    
                    # Add context-specific warnings
                    if drug == "lenalidomide":
                        warnings.append(guidelines["warnings"]["adjust_renal"])
                    elif drug == "bortezomib":
                        warnings.append(guidelines["warnings"]["neuropathy"])
                        warnings.append(guidelines["warnings"]["subcutaneous"])
                    elif drug == "daratumumab":
                        warnings.append(guidelines["warnings"]["infusion_reactions"])
                    elif drug == "carfilzomib":
                        warnings.append(guidelines["warnings"]["cardiac"])
                        warnings.append(guidelines["warnings"]["hydration"])
                    elif drug == "dexamethasone":
                        warnings.append(guidelines["warnings"]["elderly"])
                        warnings.append(guidelines["warnings"]["monitoring"])
                    elif drug == "cyclophosphamide":
                        warnings.append(guidelines["warnings"]["hydration"])
                    elif drug == "melphalan" and max_dose > 40:
                        warnings.append(guidelines["warnings"]["high_dose"])
        
        # Check for drug combinations
        has_imid = any(drug in components_lower for drug in ["lenalidomide", "pomalidomide", "thalidomide"])
        has_pi = any(drug in components_lower for drug in ["bortezomib", "carfilzomib", "ixazomib"])
        has_dara = "daratumumab" in components_lower
        has_dex = "dexamethasone" in components_lower
        
        if has_imid and not has_dex:
            warnings.append(
                "💡 IMiD regimens typically include dexamethasone for synergy. "
                "Consider adding dexamethasone 40mg weekly (20mg for elderly)."
            )
        
        if has_imid and has_pi and has_dara:
            warnings.append(
                "⚠️ Quad-therapy (IMiD + PI + Daratumumab + Dex) is intensive. "
                "Reserve for high-risk or aggressive disease. Monitor toxicity closely."
            )
        
        if has_imid:
            warnings.append(
                "⚠️ IMiDs (lenalidomide, pomalidomide) increase VTE risk. "
                "Thromboprophylaxis required: aspirin 81-325mg or LMWH for high-risk patients."
            )
        
        # Store warnings for display
        if warnings:
            self._component_warnings = warnings
        
        return components
    
    def clean(self):
        """Cross-field validation and contraindication checking."""
        cleaned_data = super().clean()
        
        intent = cleaned_data.get("intent")
        components = cleaned_data.get("components", "").lower()
        
        # Validate intent matches components
        if intent == "maintenance" and "daratumumab" in components:
            self.add_warning(
                "⚠️ Daratumumab maintenance is not standard. Consider lenalidomide maintenance (10-15mg days 1-21)."
            )
        
        if intent == "maintenance" and "bortezomib" in components:
            self.add_warning(
                "💡 Bortezomib maintenance is typically weekly SC at reduced dose (1.3mg/m² every 2 weeks)."
            )
        
        if intent == "curative" and "cyclophosphamide" in components and "melphalan" not in components:
            self.add_warning(
                "💡 For transplant-eligible patients, avoid melphalan before stem cell collection (stem cell toxicity)."
            )
        
        # Add component warnings to form errors
        if hasattr(self, "_component_warnings"):
            for warning in self._component_warnings:
                self.add_warning(warning)
        
        return cleaned_data
    
    def add_warning(self, message: str):
        """Add a non-blocking warning message."""
        if not hasattr(self, "_warnings"):
            self._warnings = []
        self._warnings.append(message)


class ScenarioForm(BootstrapValidationMixin, forms.ModelForm):
    """
    Enhanced editor form for scenarios with comprehensive clinical validation.
    
    This form enforces meaningful data entry by:
    - Validating physiological ranges for all laboratory values
    - Enforcing logical relationships between clinical parameters
    - Providing educational error messages with clinical rationale
    - Calculating difficulty scores automatically
    """

    recommended_regimens = forms.ModelMultipleChoiceField(
        queryset=Regimen.objects.order_by("name"),
        required=False,
        widget=forms.SelectMultiple(attrs={"class": "form-select", "size": "8"}),
    )

    class Meta:
        model = models.Scenario
        fields = [
            # Basic fields
            "title",
            "clinical_stage",
            "summary",
            "guideline_notes",
            "recommended_regimens",
            "expected_response",
            "active",
            # Patient characteristics
            "patient_age",
            "ecog_performance_status",
            "charlson_comorbidity_index",
            "patient_archetype",
            # Cytogenetics
            "del17p",
            "t_4_14",
            "t_14_16",
            "gain_1q21",
            "hyperdiploid",
            "t_11_14",
            # Tumor biology
            "tumor_cell_count",
            "tumor_growth_rate",
            "carrying_capacity",
            # Laboratory values
            "creatinine_clearance",
            "albumin",
            "beta2_microglobulin",
            "ldh",
            "hemoglobin",
            "calcium",
            # Risk stratification
            "riss_stage",
            # Calculated fields (read-only in form)
            "difficulty_score",
            "difficulty_level",
        ]
        widgets = {
            "title": forms.TextInput(attrs={"class": "form-control", "placeholder": "e.g., Newly Diagnosed High-Risk MM"}),
            "clinical_stage": forms.Select(attrs={"class": "form-select"}),
            "summary": forms.Textarea(attrs={"class": "form-control", "rows": 4, "placeholder": "Brief clinical description..."}),
            "guideline_notes": forms.Textarea(attrs={"class": "form-control", "rows": 4}),
            "active": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            # Patient characteristics
            "patient_age": forms.NumberInput(attrs={"class": "form-control", "min": "18", "max": "120"}),
            "ecog_performance_status": forms.Select(attrs={"class": "form-select"}),
            "charlson_comorbidity_index": forms.NumberInput(attrs={"class": "form-control", "min": "0", "max": "10"}),
            "patient_archetype": forms.Select(attrs={"class": "form-select"}),
            # Cytogenetics (checkboxes)
            "del17p": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            "t_4_14": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            "t_14_16": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            "gain_1q21": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            "hyperdiploid": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            "t_11_14": forms.CheckboxInput(attrs={"class": "form-check-input"}),
            # Tumor biology
            "tumor_cell_count": forms.NumberInput(attrs={"class": "form-control", "step": "any", "placeholder": "e.g., 1e10"}),
            "tumor_growth_rate": forms.NumberInput(attrs={"class": "form-control", "step": "0.001", "min": "0.001", "max": "0.1"}),
            "carrying_capacity": forms.NumberInput(attrs={"class": "form-control", "step": "any", "placeholder": "e.g., 1e12"}),
            # Laboratory values
            "creatinine_clearance": forms.NumberInput(attrs={"class": "form-control", "step": "1", "min": "5", "max": "150"}),
            "albumin": forms.NumberInput(attrs={"class": "form-control", "step": "0.1", "min": "1.0", "max": "6.0"}),
            "beta2_microglobulin": forms.NumberInput(attrs={"class": "form-control", "step": "0.1", "min": "0.5", "max": "20.0"}),
            "ldh": forms.NumberInput(attrs={"class": "form-control", "step": "1", "min": "50", "max": "1000"}),
            "hemoglobin": forms.NumberInput(attrs={"class": "form-control", "step": "0.1", "min": "4.0", "max": "20.0"}),
            "calcium": forms.NumberInput(attrs={"class": "form-control", "step": "0.1", "min": "6.0", "max": "18.0"}),
            # Risk stratification
            "riss_stage": forms.Select(attrs={"class": "form-select"}),
            # Calculated fields (disabled)
            "difficulty_score": forms.NumberInput(attrs={"class": "form-control", "readonly": "readonly", "disabled": "disabled"}),
            "difficulty_level": forms.Select(attrs={"class": "form-select", "disabled": "disabled"}),
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        expected_field = self.fields["expected_response"]
        expected_field.choices = [("", "Select response...")] + list(Assessment.RESPONSE_CHOICES)
        expected_field.widget.attrs.setdefault("class", "form-select")
        
        # Make calculated fields display-only
        self.fields["difficulty_score"].required = False
        self.fields["difficulty_level"].required = False
    
    def clean_patient_age(self):
        """Validate patient age is in reasonable range."""
        age = self.cleaned_data.get("patient_age")
        if age is not None:
            if age < 18:
                raise forms.ValidationError(
                    "⚠️ Age must be at least 18 years. Multiple myeloma typically affects adults aged 65+."
                )
            if age > 120:
                raise forms.ValidationError(
                    "⚠️ Age must be 120 years or less. Please verify patient demographics."
                )
            if age > 100:
                self.add_warning(
                    "⚠️ Age >100 is extremely rare. Please double-check this value."
                )
        return age
    
    def clean_tumor_cell_count(self):
        """Validate tumor cell count is in reasonable range."""
        count = self.cleaned_data.get("tumor_cell_count")
        if count is not None:
            if count < 1e6:
                raise forms.ValidationError(
                    "💡 Tumor cell count must be at least 1×10⁶ cells. "
                    "Try: 1e6 for minimal disease, 1e9 for smoldering MM, 1e10 for newly diagnosed MM, 1e11+ for advanced disease."
                )
            if count > 1e13:
                raise forms.ValidationError(
                    "⚠️ Tumor cell count exceeds physiological maximum (~10¹³ cells). "
                    "This would exceed total body cell count. Please verify."
                )
            if count > 1e12:
                self.add_warning(
                    "⚠️ Very high tumor burden (>10¹² cells). Consider plasma cell leukemia or aggressive disease."
                )
        return count
    
    def clean_tumor_growth_rate(self):
        """Validate tumor growth rate."""
        rate = self.cleaned_data.get("tumor_growth_rate")
        if rate is not None:
            if rate < 0.001:
                raise forms.ValidationError(
                    "💡 Growth rate too low (<0.001/day). "
                    "Try: 0.005 for indolent, 0.01 for typical MM, 0.05 for aggressive disease."
                )
            if rate > 0.1:
                raise forms.ValidationError(
                    "⚠️ Growth rate exceeds maximum observed in myeloma (>0.1/day). "
                    "Even aggressive disease rarely exceeds 0.05/day."
                )
        return rate
    
    def clean_creatinine_clearance(self):
        """Validate creatinine clearance."""
        crcl = self.cleaned_data.get("creatinine_clearance")
        if crcl is not None:
            if crcl < 5:
                raise forms.ValidationError(
                    "⚠️ Creatinine clearance <5 mL/min indicates end-stage renal disease requiring dialysis. "
                    "IMiD dosing must be significantly adjusted."
                )
            if crcl < 30:
                self.add_warning(
                    "💡 Severe renal impairment (CrCl <30). Lenalidomide dose should not exceed 15mg every other day."
                )
            elif crcl < 60:
                self.add_warning(
                    "💡 Moderate renal impairment (CrCl 30-60). Consider lenalidomide dose reduction to 10-15mg daily."
                )
        return crcl
    
    def clean_albumin(self):
        """Validate serum albumin."""
        albumin = self.cleaned_data.get("albumin")
        if albumin is not None:
            if albumin < 1.0:
                raise forms.ValidationError(
                    "⚠️ Severe hypoalbuminemia (<1.0 g/dL) indicates critical illness or synthetic liver dysfunction."
                )
            if albumin < 3.0:
                self.add_warning(
                    "💡 Low albumin (<3.0 g/dL) affects ISS staging and prognosis. Consider nutritional support."
                )
        return albumin
    
    def clean_beta2_microglobulin(self):
        """Validate β2-microglobulin."""
        b2m = self.cleaned_data.get("beta2_microglobulin")
        if b2m is not None:
            if b2m > 10.0:
                self.add_warning(
                    "⚠️ Very elevated β2M (>10 mg/L) indicates ISS stage III and poor prognosis."
                )
            elif b2m > 5.5:
                self.add_warning(
                    "💡 Elevated β2M (>5.5 mg/L) contributes to R-ISS III classification."
                )
        return b2m
    
    def clean_hemoglobin(self):
        """Validate hemoglobin."""
        hgb = self.cleaned_data.get("hemoglobin")
        if hgb is not None:
            if hgb < 6.0:
                self.add_warning(
                    "⚠️ Severe anemia (Hgb <6 g/dL). Transfusion support typically required."
                )
            elif hgb < 10.0:
                self.add_warning(
                    "💡 Anemia (Hgb <10 g/dL) is common in MM. Consider erythropoietin if symptomatic."
                )
        return hgb
    
    def clean_calcium(self):
        """Validate serum calcium."""
        ca = self.cleaned_data.get("calcium")
        if ca is not None:
            if ca > 14.0:
                self.add_warning(
                    "⚠️ Severe hypercalcemia (>14 mg/dL) requires immediate intervention: IV fluids, bisphosphonates, calcitonin."
                )
            elif ca > 11.5:
                self.add_warning(
                    "💡 Moderate hypercalcemia (>11.5 mg/dL). Ensure adequate hydration and consider bisphosphonates."
                )
        return ca
    
    def clean(self):
        """Cross-field validation and educational warnings."""
        cleaned_data = super().clean()
        
        # Check for high-risk cytogenetics
        del17p = cleaned_data.get("del17p")
        t_4_14 = cleaned_data.get("t_4_14")
        t_14_16 = cleaned_data.get("t_14_16")
        gain_1q21 = cleaned_data.get("gain_1q21")
        
        high_risk_count = sum([del17p, t_4_14, t_14_16, gain_1q21])
        
        if high_risk_count >= 2:
            self.add_warning(
                "⚠️ Multiple high-risk cytogenetic abnormalities detected. "
                "Consider quad-therapy (e.g., Dara-VRd) or clinical trial enrollment."
            )
        
        # Validate R-ISS staging consistency with labs
        riss_stage = cleaned_data.get("riss_stage")
        albumin = cleaned_data.get("albumin")
        b2m = cleaned_data.get("beta2_microglobulin")
        ldh = cleaned_data.get("ldh")
        
        if riss_stage and albumin and b2m:
            # R-ISS I: ISS I + standard-risk CA + normal LDH
            # ISS I: β2M <3.5 AND albumin ≥3.5
            iss_1 = b2m < 3.5 and albumin >= 3.5
            iss_3 = b2m >= 5.5
            
            if riss_stage == "I" and iss_3:
                self.add_warning(
                    "💡 R-ISS I staging inconsistent with β2M ≥5.5 (suggests ISS III). "
                    "Please verify staging or laboratory values."
                )
            
            if riss_stage == "III" and iss_1:
                self.add_warning(
                    "💡 R-ISS III staging inconsistent with ISS I parameters. "
                    "R-ISS III requires ISS III OR high-risk CA OR elevated LDH."
                )
        
        # Warn about ECOG and Charlson mismatch
        ecog = cleaned_data.get("ecog_performance_status")
        charlson = cleaned_data.get("charlson_comorbidity_index")
        
        if ecog is not None and charlson is not None:
            if ecog >= 3 and charlson <= 1:
                self.add_warning(
                    "⚠️ High ECOG score (≥3) but low comorbidity index (≤1) is unusual. "
                    "Consider if MM-related disability vs comorbidities."
                )
        
        # Calculate difficulty score if tumor parameters are present
        tumor_count = cleaned_data.get("tumor_cell_count")
        growth_rate = cleaned_data.get("tumor_growth_rate")
        
        if tumor_count and growth_rate:
            # Import difficulty scoring here to avoid circular imports
            try:
                from .difficulty_scoring import (
                    TumorBurdenScore,
                    GrowthRateScore,
                    CytogeneticRiskScore,
                    FrailtyScore,
                    StageScore,
                    DifficultyCalculator
                )
                
                # Calculate component scores
                tumor_burden = TumorBurdenScore(tumor_cells=tumor_count)
                growth_score = GrowthRateScore(growth_rate=growth_rate)
                
                # Cytogenetic risk
                cyto_score = CytogeneticRiskScore(
                    del17p=bool(del17p),
                    t_4_14=bool(t_4_14),
                    t_14_16=bool(t_14_16),
                    gain_1q21=bool(gain_1q21),
                )
                
                # Frailty from ECOG + Charlson
                frailty_score = FrailtyScore(
                    age=cleaned_data.get("patient_age") or 65,
                    ecog=ecog or 0,
                    charlson=charlson or 0,
                )
                
                # Stage score from R-ISS
                stage_map = {"I": 1, "II": 2, "III": 3, "": 2}
                stage_score = StageScore(stage=stage_map.get(riss_stage, 2))
                
                # Calculate total difficulty
                calculator = DifficultyCalculator(
                    tumor_burden=tumor_burden,
                    growth_rate=growth_score,
                    cytogenetics=cyto_score,
                    frailty=frailty_score,
                    stage=stage_score,
                )
                
                total_score = calculator.calculate_total_score()
                difficulty_level = calculator.categorize_difficulty(total_score)
                
                # Update cleaned data
                cleaned_data["difficulty_score"] = total_score
                cleaned_data["difficulty_level"] = difficulty_level
                
            except ImportError:
                pass  # Difficulty scoring not available
        
        return cleaned_data
    
    def add_warning(self, message: str):
        """Add a non-blocking warning message."""
        if not hasattr(self, "_warnings"):
            self._warnings = []
        self._warnings.append(message)


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
    cohort_size = forms.TypedChoiceField(
        choices=[(1, "1"), (10, "10"), (50, "50"), (200, "200")],
        coerce=int,
        initial=1,
        label="Virtual cohort size",
        widget=forms.Select(attrs={"class": "form-select"}),
        help_text="Number of virtual subjects sampled for uncertainty bands (larger takes longer).",
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
    use_twin = forms.BooleanField(
        required=False,
        initial=True,
        widget=forms.CheckboxInput(attrs={"class": "form-check-input"}),
        label="Use Patient Twin",
        help_text="Auto-derive biology from linked laboratory assessments when available.",
    )
    seed = forms.IntegerField(
        required=False,
        min_value=0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "1",
                "min": "0",
                "inputmode": "numeric",
            }
        ),
        label="Random seed",
        help_text="Leave blank for stochastic runs; set a value to reproduce cohorts/optimizations.",
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
        # 🎓 Educational error messages with explanations
        if len_dose > 50:
            errors["lenalidomide_dose"] = ValidationError(
                "Lenalidomide dose must be ≤ 50 mg/day. "
                "💡 Why? Doses >50mg cause severe neutropenia (low white blood cells) in most patients. "
                "Try the standard dose of 25mg/day used in clinical trials."
            )
        if bor_dose > 2:
            errors["bortezomib_dose"] = ValidationError(
                "Bortezomib dose must be ≤ 2 mg/m² per week. "
                "💡 Why? Higher doses dramatically increase peripheral neuropathy risk (nerve damage). "
                "Standard dosing is 1.3 mg/m² which balances efficacy and safety."
            )
        if dara_dose > 20:
            errors["daratumumab_dose"] = ValidationError(
                "Daratumumab dose must be ≤ 20 mg/kg. "
                "💡 Why? Doses beyond 20mg/kg don't improve outcomes but increase infusion reactions. "
                "Clinical trials use 16 mg/kg loading dose."
            )
        if time_horizon > 365:
            errors["time_horizon"] = ValidationError(
                "Simulation horizon is limited to 365 days. "
                "💡 Why? Long-term predictions become unreliable beyond 1 year due to disease evolution. "
                "Try 180 days (6 months) for a realistic treatment cycle."
            )
        if tumor_growth > 0.1:
            errors["tumor_growth_rate"] = ValidationError(
                "Tumor growth rate must be ≤ 0.10 day⁻¹. "
                "💡 Why? This represents doubling time <7 days, extremely aggressive even for advanced myeloma. "
                "Typical range is 0.015-0.03 day⁻¹ (doubling every 23-46 days)."
            )
        if healthy_growth > 0.05:
            errors["healthy_growth_rate"] = ValidationError(
                "Healthy growth rate must be ≤ 0.05 day⁻¹. "
                "💡 Why? Bone marrow regenerates slower than tumor cells. "
                "Use 0.01-0.02 day⁻¹ for realistic marrow recovery."
            )
        if interaction_strength > 0.2:
            errors["interaction_strength"] = ValidationError(
                "Interaction strength must be ≤ 0.2. "
                "💡 Why? Values >0.2 imply unrealistic competition between cell types. "
                "Try 0.05-0.15 for moderate interaction effects."
            )

        if creatinine is not None:
            if creatinine < 30 and len_dose > 10:
                errors["lenalidomide_dose"] = ValidationError(
                    "Creatinine clearance <30 ml/min requires lenalidomide dose ≤10 mg. "
                    "💡 Why? Kidneys eliminate lenalidomide—poor kidney function causes drug accumulation and toxicity. "
                    "Dose reduction prevents life-threatening side effects."
                )
            elif creatinine < 60 and len_dose > 15:
                errors["lenalidomide_dose"] = ValidationError(
                    "Creatinine clearance <60 ml/min requires lenalidomide dose ≤15 mg. "
                    "💡 Why? Moderately impaired kidneys need dose adjustment to avoid myelosuppression."
                )

        if neuropathy >= 2 and bor_dose > 1.0:
            errors["bortezomib_dose"] = ValidationError(
                "Peripheral neuropathy grade ≥2 mandates bortezomib dose ≤1.0 mg/m². "
                "💡 Why? Bortezomib causes nerve damage—continuing high doses with existing neuropathy leads to irreversible disability. "
                "Dose reduction or drug holiday preserves quality of life."
            )

        if anc is not None and anc < 1.0:
            raise ValidationError(
                "Simulation blocked: absolute neutrophil count <1.0 ×10⁹/L (myelosuppression). "
                "💡 Why? Dangerously low white blood cells mean infection risk is too high to start treatment. "
                "Use growth factors (G-CSF) to raise counts first."
            )
        if platelets is not None and platelets < 75:
            raise ValidationError(
                "Simulation blocked: platelets <75 ×10⁹/L (thrombocytopenia). "
                "💡 Why? Low platelets cause bleeding risk—treatment would worsen this. "
                "Wait for platelet recovery >75 or use platelet transfusion."
            )

        if len_dose > 40 and bor_dose > 1.5:
            raise ValidationError(
                "Combined high doses of lenalidomide (>40 mg) and bortezomib (>1.5 mg/m²) may exceed safe toxicity bounds. "
                "💡 Why? Both drugs suppress bone marrow—combining high doses compounds neutropenia and thrombocytopenia risk. "
                "Reduce at least one drug to stay within safe margins."
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

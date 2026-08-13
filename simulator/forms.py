from __future__ import annotations

from django import forms
from django.core.exceptions import ValidationError

from mmportal.forms_mixins import BootstrapValidationMixin

from clinic.models import Assessment, Regimen
from django.utils.text import slugify

from . import models
from .presets import PRESETS

SIMULATION_FORM_HELP_TEXT_EN = {
    "preset": "Research preset containing hypothetical model inputs and schedules.",
    "creatinine_clearance": "User-provided renal-state input used by configured model rules.",
    "neuropathy_grade": "User-provided CTCAE neuropathy grade used by configured model constraints.",
    "anc": "User-provided ANC; configured low-value constraints may block a research run.",
    "platelets": "User-provided platelet count; configured low-value constraints may block a research run.",
    "pregnancy": "User-provided flag used by configured research constraints.",
    "baseline_tumor_cells": "Estimated malignant plasma cell burden (cells).",
    "baseline_healthy_cells": "Approximate pool of normal plasma cells (cells).",
    "lenalidomide_dose": "Hypothetical lenalidomide model input constrained to the configured research range.",
    "bortezomib_dose": "Hypothetical bortezomib model input constrained to the configured research range.",
    "daratumumab_dose": "Hypothetical daratumumab model input constrained to the configured research range.",
    "carfilzomib_dose": "Hypothetical carfilzomib model input constrained to the configured research range.",
    "time_horizon": "Length of the virtual treatment window (days). Longer horizons increase numerical stiffness.",
    "tumor_growth_rate": "Logistic growth rate of malignant plasma cells. Values >0.05 day⁻¹ are rare.",
    "healthy_growth_rate": "Marrow recovery kinetics for healthy plasma cells (≈0.01–0.02 day⁻¹).",
    "interaction_strength": "Synergy/toxicity coupling between agents (dimensionless).",
    "cohort_size": "Number of virtual subjects sampled for uncertainty bands.",
    "use_twin": "Enable the Patient Twin to auto-derive biology from the latest labs.",
    "seed": "Optional seed for reproducible virtual cohorts and solver noise.",
    "custom_drug_enabled": "Enable an experimental custom agent using manual PK/PD parameters (editor-only).",
    "custom_drug_name": "Display name for the experimental agent. Used only for labeling outputs.",
    "custom_drug_dose": "Total dose amount (arbitrary units). The simulator spreads exposure over the horizon unless a preset schedule exists.",
    "custom_pk_half_life": "PK half-life in hours. Used to compute elimination rate.",
    "custom_pk_vd": "Volume of distribution (arbitrary units). Used to scale initial concentration.",
    "custom_pd_emax": "Maximum effect (0–1). Higher values imply stronger tumor kill and toxicity coupling.",
    "custom_pd_ec50": "Concentration giving half-maximal effect (must be >0).",
}

SIMULATION_FORM_HELP_TEXT_IT = {
    "preset": "Preset di ricerca con input e calendari ipotetici del modello.",
    "creatinine_clearance": "Input sullo stato renale fornito dall'utente e usato dalle regole del modello.",
    "neuropathy_grade": "Grado CTCAE fornito dall'utente e usato dai vincoli del modello.",
    "anc": "ANC fornito dall'utente; vincoli configurati possono bloccare la simulazione.",
    "platelets": "Piastrine fornite dall'utente; vincoli configurati possono bloccare la simulazione.",
    "pregnancy": "Flag fornito dall'utente e usato dai vincoli di ricerca configurati.",
    "baseline_tumor_cells": "Stima del carico di cellule plasmatiche maligne (cellule).",
    "baseline_healthy_cells": "Pool approssimativo di cellule plasmatiche sane (cellule).",
    "lenalidomide_dose": "Input ipotetico di lenalidomide limitato all'intervallo di ricerca configurato.",
    "bortezomib_dose": "Input ipotetico di bortezomib limitato all'intervallo di ricerca configurato.",
    "daratumumab_dose": "Input ipotetico di daratumumab limitato all'intervallo di ricerca configurato.",
    "carfilzomib_dose": "Input ipotetico di carfilzomib limitato all'intervallo di ricerca configurato.",
    "time_horizon": "Durata della finestra terapeutica virtuale (giorni). Orizzonti lunghi aumentano l’incertezza.",
    "tumor_growth_rate": "Velocità logistica di crescita del tumore. Valori >0.05 day⁻¹ sono rari.",
    "healthy_growth_rate": "Cinetica di recupero per le cellule sane. Tipicamente 0.01–0.02 day⁻¹.",
    "interaction_strength": "Ampiezza della sinergia/tossicità tra i farmaci (senza unità).",
    "cohort_size": "Numero di pazienti virtuali campionati per stimare le bande di incertezza.",
    "use_twin": "Attiva il Gemello Paziente per riempire automaticamente i parametri biologici dai laboratori disponibili.",
    "seed": "Seed opzionale per rendere ripetibili coorti virtuali e solver.",
    "custom_drug_enabled": "Abilita un agente sperimentale personalizzato con parametri PK/PD manuali (solo editor).",
    "custom_drug_name": "Nome mostrato per l’agente sperimentale. Usato solo per etichettare gli output.",
    "custom_drug_dose": "Dose totale (unità arbitrarie). Il simulatore distribuisce l’esposizione sull’orizzonte salvo schedula preset.",
    "custom_pk_half_life": "Emivita PK in ore. Usata per calcolare la velocità di eliminazione.",
    "custom_pk_vd": "Volume di distribuzione (unità arbitrarie). Serve a scalare la concentrazione iniziale.",
    "custom_pd_emax": "Effetto massimo (0–1). Valori più alti implicano maggiore efficacia e coupling di tossicità.",
    "custom_pd_ec50": "Concentrazione a metà effetto massimo (deve essere >0).",
}


class RegimenForm(BootstrapValidationMixin, forms.ModelForm):
    """
    Editor form for recording research regimen metadata.
    
    Numeric checks below enforce configured research-input ranges only. They do
    not validate clinical safety or provide patient-specific instructions.
    
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
            "Example research inputs: Lenalidomide 25mg, Bortezomib 1.3mg/m², "
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
                "max_exceeded": "Research input is above the configured lenalidomide range; the value is not clinically validated.",
                "adjust_renal": "Renal-state assumptions and units must be reported with this hypothetical input.",
            }
        },
        "bortezomib": {
            "min": 0.7,
            "max": 1.3,
            "standard": 1.3,
            "unit": "mg/m²",
            "warnings": {
                "max_exceeded": "Research input is above the configured bortezomib range; the value is not clinically validated.",
                "neuropathy": "Neuropathy is a user-provided constraint input; no patient action is inferred.",
                "subcutaneous": "Administration route is not identified by this research form.",
            }
        },
        "daratumumab": {
            "min": 8,
            "max": 16,
            "standard": 16,
            "unit": "mg/kg",
            "warnings": {
                "max_exceeded": "Research input is above the configured daratumumab range; the value is not clinically validated.",
                "infusion_reactions": "Infusion management is outside this research prototype's intended use.",
            }
        },
        "carfilzomib": {
            "min": 20,
            "max": 56,
            "standard": 27,  # 20/27 or 20/56 dosing
            "unit": "mg/m²",
            "warnings": {
                "max_exceeded": "Research input is above the configured carfilzomib range; the value is not clinically validated.",
                "cardiac": "Cardiac state is not sufficiently represented for a clinical conclusion.",
                "hydration": "Supportive-care instructions are outside this research prototype's intended use.",
            }
        },
        "dexamethasone": {
            "min": 4,
            "max": 40,
            "standard": 40,
            "unit": "mg",
            "warnings": {
                "max_exceeded": "Research input is above the configured dexamethasone range; the value is not clinically validated.",
                "elderly": "Age and frailty are incompletely represented research covariates.",
                "monitoring": "Monitoring instructions are outside this research prototype's intended use.",
            }
        },
        "pomalidomide": {
            "min": 2,
            "max": 4,
            "standard": 4,
            "unit": "mg",
            "warnings": {
                "max_exceeded": "Research input is above the configured pomalidomide range; the value is not clinically validated.",
                "adjust_renal": "Renal-state assumptions must be reported with this hypothetical input.",
            }
        },
        "cyclophosphamide": {
            "min": 50,
            "max": 500,
            "standard": 300,
            "unit": "mg/m²",
            "warnings": {
                "max_exceeded": "Research input is above the configured cyclophosphamide range; the value is not clinically validated.",
                "hydration": "Supportive-care instructions are outside this research prototype's intended use.",
            }
        },
        "melphalan": {
            "min": 5,
            "max": 200,
            "standard": 9,  # oral dosing
            "unit": "mg/m²",
            "warnings": {
                "high_dose": "Research input is above the configured melphalan exploratory threshold.",
                "myelosuppression": "Clinical monitoring instructions are outside this research prototype's intended use.",
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
                "The research regimen metadata omits a component used by the legacy heuristic rule."
            )
        
        if has_imid and has_pi and has_dara:
            warnings.append(
                "This multi-component research scenario requires explicit assumptions and uncertainty reporting."
            )
        
        if has_imid:
            warnings.append(
                "Supportive-care and patient-management instructions are outside this research prototype."
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
                "The selected intent/component combination is an unvalidated research scenario."
            )
        
        if intent == "maintenance" and "bortezomib" in components:
            self.add_warning(
                "The selected intent/component combination is an unvalidated research scenario."
            )
        
        if intent == "curative" and "cyclophosphamide" in components and "melphalan" not in components:
            self.add_warning(
                "Transplant sequencing is outside this research prototype's intended use."
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
    Editor form for synthetic and research scenarios.
    
    This form enforces meaningful data entry by:
    - Validating physiological ranges for all laboratory values
    - Enforcing logical relationships between clinical parameters
    - Reporting model-input validation without clinical instructions
    - Calculating difficulty scores automatically
    """

    recommended_regimens = forms.ModelMultipleChoiceField(
        label="Associated exploratory regimen references",
        queryset=Regimen.objects.order_by("name"),
        required=False,
        help_text="Catalog links for research comparison; not patient-specific selections.",
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
        expected_field.label = "Stored synthetic response label"
        expected_field.help_text = "Optional scenario metadata; not a patient prediction."
        
        # Make calculated fields display-only
        self.fields["difficulty_score"].required = False
        self.fields["difficulty_level"].required = False
    
    def clean_patient_age(self):
        """Validate patient age is in reasonable range."""
        age = self.cleaned_data.get("patient_age")
        if age is not None:
            if age < 18:
                raise forms.ValidationError(
                    "Age must be at least 18 years for this configured synthetic-scenario domain."
                )
            if age > 120:
                raise forms.ValidationError(
                    "Age must be 120 years or less for this configured synthetic-scenario domain."
                )
            if age > 100:
                self.add_warning(
                    "Age >100 activates a synthetic-input verification flag."
                )
        return age
    
    def clean_tumor_cell_count(self):
        """Validate tumor cell count is in reasonable range."""
        count = self.cleaned_data.get("tumor_cell_count")
        if count is not None:
            if count < 1e6:
                raise forms.ValidationError(
                    "Tumor-state input must be at least 1×10⁶ cells in this model."
                )
            if count > 1e13:
                raise forms.ValidationError(
                    "Tumor-state input exceeds the configured model maximum (10¹³ cells)."
                )
            if count > 1e12:
                self.add_warning(
                    "Very high model input (>10¹² cells); record the assumption and test sensitivity."
                )
        return count
    
    def clean_tumor_growth_rate(self):
        """Validate tumor growth rate."""
        rate = self.cleaned_data.get("tumor_growth_rate")
        if rate is not None:
            if rate < 0.001:
                raise forms.ValidationError(
                    "Growth-rate input is below the configured model domain (0.001/day)."
                )
            if rate > 0.1:
                raise forms.ValidationError(
                    "Growth-rate input exceeds the configured model domain (0.1/day)."
                )
        return rate
    
    def clean_creatinine_clearance(self):
        """Validate creatinine clearance."""
        crcl = self.cleaned_data.get("creatinine_clearance")
        if crcl is not None:
            if crcl < 5:
                raise forms.ValidationError(
                    "Creatinine clearance <5 mL/min is outside the configured research domain."
                )
            if crcl < 30:
                self.add_warning(
                    "CrCl <30 activates an unvalidated renal-state constraint flag."
                )
            elif crcl < 60:
                self.add_warning(
                    "CrCl 30-60 activates an unvalidated renal-state watch flag."
                )
        return crcl
    
    def clean_albumin(self):
        """Validate serum albumin."""
        albumin = self.cleaned_data.get("albumin")
        if albumin is not None:
            if albumin < 1.0:
                raise forms.ValidationError(
                    "Albumin input below 1.0 g/dL is outside the configured research domain."
                )
            if albumin < 3.0:
                self.add_warning(
                    "Albumin <3.0 g/dL changes the model's derived risk context."
                )
        return albumin
    
    def clean_beta2_microglobulin(self):
        """Validate β2-microglobulin."""
        b2m = self.cleaned_data.get("beta2_microglobulin")
        if b2m is not None:
            if b2m > 10.0:
                self.add_warning(
                    "β2M >10 mg/L activates a high-value input verification flag."
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
                    "Hgb <6 g/dL is outside the calibrated domain of this research prototype."
                )
            elif hgb < 10.0:
                self.add_warning(
                    "Hgb <10 g/dL activates a descriptive host-state flag."
                )
        return hgb
    
    def clean_calcium(self):
        """Validate serum calcium."""
        ca = self.cleaned_data.get("calcium")
        if ca is not None:
            if ca > 14.0:
                self.add_warning(
                    "Calcium >14 mg/dL is outside the configured research domain and requires source verification before simulation."
                )
            elif ca > 11.5:
                self.add_warning(
                    "Calcium >11.5 mg/dL activates a descriptive host-state flag."
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
                "This combination activates a high-risk research-context flag."
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
                    "The model cannot distinguish disease-related disability from comorbidity."
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

    preset = forms.ChoiceField(choices=PRESET_CHOICES, required=False, label="Hypothetical model preset")

    creatinine_clearance = forms.FloatField(
        required=False,
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
        help_text="User-provided renal-state input used by configured model rules.",
    )
    neuropathy_grade = forms.ChoiceField(
        required=False,
        choices=[(str(i), f"Grade {i}") for i in range(0, 4)],
        label="Peripheral neuropathy (grade)",
        widget=forms.Select(attrs={"class": "form-select"}),
        help_text="User-provided CTCAE grade (0–3); configured combinations may be outside the research domain.",
    )
    anc = forms.FloatField(
        required=False,
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
        help_text="ANC <1.0 ×10⁹/L is outside the configured research domain and blocks simulation.",
    )
    platelets = forms.FloatField(
        required=False,
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
        help_text="Platelets <75 ×10⁹/L are outside the configured research domain and block simulation.",
    )
    pregnancy = forms.BooleanField(
        required=False,
        widget=forms.CheckboxInput(attrs={"class": "form-check-input"}),
        label="Pregnancy flagged",
        help_text="The configured research protocol excludes this input state and blocks simulation when flagged.",
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
        help_text="Heuristic initial malignant-cell model state (cells); not a measured patient count.",
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
        help_text="Heuristic initial healthy-cell model state (cells); not a validated marrow endpoint.",
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
        help_text="Hypothetical exposure input constrained to the configured research range.",
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
        help_text="Hypothetical bortezomib input constrained to the configured research range.",
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
        help_text="Hypothetical exposure input constrained to the configured research range.",
    )
    carfilzomib_dose = forms.FloatField(
        required=False,
        min_value=0.0,
        initial=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0",
                "max": "70",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Carfilzomib dose (mg/m²)",
        help_text="Hypothetical exposure input constrained to the configured research range.",
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
        help_text="Length of the virtual model window (days); bounds support numerical stability only.",
    )
    cohort_size = forms.TypedChoiceField(
        choices=[(1, "1"), (10, "10"), (50, "50"), (200, "200")],
        coerce=int,
        required=False,
        empty_value=1,
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
        help_text="Heuristic recovery rate for the simplified healthy-cell state.",
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

    custom_drug_enabled = forms.BooleanField(
        required=False,
        initial=False,
        widget=forms.CheckboxInput(attrs={"class": "form-check-input"}),
        label="Enable custom drug (experimental)",
        help_text="Editor-only: simulate an experimental agent via manual PK/PD parameters.",
    )
    custom_drug_name = forms.CharField(
        required=False,
        max_length=40,
        widget=forms.TextInput(
            attrs={
                "class": "form-control",
                "placeholder": "e.g., ExperimentalAgent",
            }
        ),
        label="Custom drug name",
        help_text="Used for labeling; does not create a preset.",
    )
    custom_drug_dose = forms.FloatField(
        required=False,
        min_value=0.0,
        initial=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0",
                "max": "1000",
                "inputmode": "decimal",
            }
        ),
        label="Custom drug dose (total)",
        help_text="Total administered amount across the horizon (arbitrary units).",
    )
    custom_pk_half_life = forms.FloatField(
        required=False,
        min_value=0.1,
        initial=24.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0.1",
                "max": "2000",
                "inputmode": "decimal",
            }
        ),
        label="Custom PK: half-life (hours)",
        help_text="PK elimination half-life (hours).",
    )
    custom_pk_vd = forms.FloatField(
        required=False,
        min_value=0.1,
        initial=30.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.1",
                "min": "0.1",
                "max": "10000",
                "inputmode": "decimal",
            }
        ),
        label="Custom PK: Vd",
        help_text="Volume of distribution (arbitrary units).",
    )
    custom_pd_emax = forms.FloatField(
        required=False,
        min_value=0.0,
        max_value=1.0,
        initial=0.5,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.01",
                "min": "0",
                "max": "1",
                "inputmode": "decimal",
            }
        ),
        label="Custom PD: Emax (0–1)",
        help_text="Max effect (0–1).",
    )
    custom_pd_ec50 = forms.FloatField(
        required=False,
        min_value=1e-6,
        initial=1.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.01",
                "min": "0",
                "max": "10000",
                "inputmode": "decimal",
            }
        ),
        label="Custom PD: EC50",
        help_text="Half-maximal concentration (must be >0).",
    )

    twin_assessment_id = forms.TypedChoiceField(
        required=False,
        coerce=int,
        choices=[("", "— Select assessment —")],
        widget=forms.Select(attrs={"class": "form-select"}),
        label="Patient Twin assessment",
        help_text="Select the lab snapshot (Assessment) used to derive Patient Twin biology.",
    )
    twin_biology_mode = forms.ChoiceField(
        required=False,
        initial="auto",
        choices=[("auto", "Auto"), ("manual", "Manual")],
        widget=forms.Select(attrs={"class": "form-select"}),
        label="Twin-driven biology",
        help_text="Auto uses Patient Twin-derived biology. Manual keeps your entered values.",
    )

    guided_tumor_aggressiveness = forms.ChoiceField(
        required=False,
        choices=[
            ("", "— No guided adjustment —"),
            ("lower", "Lower"),
            ("typical", "Typical"),
            ("higher", "Higher"),
        ],
        widget=forms.Select(attrs={"class": "form-select"}),
        label="Guided choice: tumor aggressiveness",
        help_text="Optional. Applies a small adjustment to tumor growth rate without manual typing.",
    )
    guided_immune_status = forms.ChoiceField(
        required=False,
        choices=[
            ("", "— No guided adjustment —"),
            ("better", "Better"),
            ("typical", "Typical"),
            ("worse", "Worse"),
        ],
        widget=forms.Select(attrs={"class": "form-select"}),
        label="Guided choice: immune status",
        help_text="Optional. Applies a small adjustment to immune compromise index without manual typing.",
    )
    carrying_capacity_tumor = forms.FloatField(
        required=False,
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "1",
                "min": "0",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Tumor carrying capacity (cells)",
        help_text="Maximum tumor burden for logistic growth. Leave on Auto for Twin-driven value.",
    )
    carrying_capacity_healthy = forms.FloatField(
        required=False,
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "1",
                "min": "0",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Healthy carrying capacity (cells)",
        help_text="Maximum healthy plasma cell pool for logistic growth. Leave on Auto for Twin-driven value.",
    )
    immune_compromise_index = forms.FloatField(
        required=False,
        min_value=0.0,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control",
                "step": "0.01",
                "min": "0",
                "max": "5",
                "inputmode": "decimal",
                "readonly": "readonly",
            }
        ),
        label="Immune compromise index",
        help_text="Dimensionless immune impairment multiplier. Leave on Auto for Twin-derived value.",
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
        self.user = kwargs.pop("user", None)
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
            "carfilzomib_dose",
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

        if "carrying_capacity_tumor" not in self.data:
            self.fields["carrying_capacity_tumor"].initial = 1.0e9
        if "carrying_capacity_healthy" not in self.data:
            self.fields["carrying_capacity_healthy"].initial = 5.0e11
        if "immune_compromise_index" not in self.data:
            self.fields["immune_compromise_index"].initial = 1.0

        self._init_twin_assessment_choices()

        self.slider_bounds = self._compute_slider_bounds(defaults, slider_fields)

    def _init_twin_assessment_choices(self) -> None:
        """Populate assessment choices.

        Ownership policy (PR1.1):
        - user sees assessments for patients they own
        - staff/editor sees all
        """
        if not self.user or not getattr(self.user, "is_authenticated", False):
            self.fields["twin_assessment_id"].choices = [("", "— Select assessment —")]
            self.fields["twin_assessment_id"].disabled = True
            return

        from simulator.permissions import is_editor

        is_privileged = getattr(self.user, "is_staff", False) or is_editor(self.user)

        try:
            from clinic.models import Assessment
            from django.db.models import Q

            qs = Assessment.objects.select_related("patient").order_by("-date")
            if not is_privileged:
                qs = qs.filter(Q(patient__owner=self.user) | Q(patient__mrn__startswith="DEMO"))
        except Exception:
            self.fields["twin_assessment_id"].choices = [("", "— Select assessment —")]
            return

        if not is_privileged and not qs.exists():
            self.fields["twin_assessment_id"].choices = [("", "— Select assessment —")]
            self.fields["twin_assessment_id"].disabled = True
            return

        choices: list[tuple[object, str]] = [("", "— Select assessment —")]
        for assessment in qs[:200]:
            first = (assessment.patient.first_name or "").strip()
            last = (assessment.patient.last_name or "").strip()
            full_name = (f"{first} {last}").strip() or last or first or "(unknown)"
            label = f"{assessment.patient.mrn} · {full_name} · {assessment.date}"
            choices.append((assessment.pk, label))
        self.fields["twin_assessment_id"].choices = choices

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

        use_twin = bool(cleaned.get("use_twin"))
        assessment_id = cleaned.get("twin_assessment_id")
        raw_assessment_id = (self.data.get("twin_assessment_id") or "").strip() if self.is_bound else ""
        biology_mode = (cleaned.get("twin_biology_mode") or "auto").lower()
        cleaned["twin_biology_mode"] = biology_mode

        if use_twin and not assessment_id and not raw_assessment_id:
            # Beginner-friendly behavior: don't hard-block runs when Twin is unavailable.
            # Twin can only work with a selected Assessment snapshot.
            cleaned["use_twin"] = False
            use_twin = False
            self.warnings.append(
                "Patient Twin disabled: select an Assessment snapshot to enable Twin-driven biology."
            )

        twin_keys = (
            "tumor_growth_rate",
            "healthy_growth_rate",
            "carrying_capacity_tumor",
            "carrying_capacity_healthy",
            "immune_compromise_index",
        )
        if use_twin and biology_mode == "auto":
            for key in twin_keys:
                cleaned[key] = None

        len_dose = self._clean_numeric(cleaned, "lenalidomide_dose")
        bor_dose = self._clean_numeric(cleaned, "bortezomib_dose")
        dara_dose = self._clean_numeric(cleaned, "daratumumab_dose")
        carf_dose = self._clean_numeric(cleaned, "carfilzomib_dose")
        time_horizon = self._clean_numeric(cleaned, "time_horizon")
        tumor_growth = self._clean_numeric(cleaned, "tumor_growth_rate")
        healthy_growth = self._clean_numeric(cleaned, "healthy_growth_rate")
        baseline_tumor = self._clean_numeric(cleaned, "baseline_tumor_cells")
        baseline_healthy = self._clean_numeric(cleaned, "baseline_healthy_cells")
        interaction_strength = self._clean_numeric(cleaned, "interaction_strength")
        creatinine = self._clean_numeric(cleaned, "creatinine_clearance")
        neuropathy = int(cleaned.get("neuropathy_grade") or (self.fields["neuropathy_grade"].initial or 0))
        anc = self._clean_numeric(cleaned, "anc")
        platelets = self._clean_numeric(cleaned, "platelets")
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
                "Lenalidomide input exceeds the configured research maximum (50 mg/day)."
            )
        if bor_dose > 2:
            errors["bortezomib_dose"] = ValidationError(
                "Bortezomib input exceeds the configured research maximum (2 mg/m²)."
            )
        if dara_dose > 20:
            errors["daratumumab_dose"] = ValidationError(
                "Daratumumab input exceeds the configured research maximum (20 mg/kg)."
            )
        if carf_dose > 70:
            errors["carfilzomib_dose"] = ValidationError(
                "Carfilzomib input exceeds the configured research maximum (70 mg/m²)."
            )
        if time_horizon > 365:
            errors["time_horizon"] = ValidationError(
                "Simulation horizon exceeds the configured research maximum (365 days)."
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
                    "The input combination crosses an unvalidated renal-state research constraint."
                )
            elif creatinine < 60 and len_dose > 15:
                errors["lenalidomide_dose"] = ValidationError(
                    "The input combination crosses an unvalidated renal-state research constraint."
                )

        if neuropathy >= 2 and bor_dose > 1.0:
            errors["bortezomib_dose"] = ValidationError(
                "The input combination crosses an unvalidated neuropathy research constraint."
            )

        if anc is not None and anc < 1.0:
            raise ValidationError(
                "Simulation blocked: ANC input is below the configured research domain."
            )
        if platelets is not None and platelets < 75:
            raise ValidationError(
                "Simulation blocked: platelet input is below the configured research domain."
            )

        if len_dose > 40 and bor_dose > 1.5:
            raise ValidationError(
                "Combined inputs exceed the configured research domain."
            )


        if interaction_strength > 0.15:
            self.warnings.append("Interaction strength above 0.15 activates an unvalidated model constraint flag.")
        if 35 < len_dose <= 40:
            self.warnings.append("Lenalidomide dose is near the upper investigational range (>35 mg).")
        if 1.3 < bor_dose <= 1.5:
            self.warnings.append("Bortezomib input exceeds the configured reference value.")
        if 16 < dara_dose <= 20:
            self.warnings.append("Daratumumab input exceeds the configured reference value.")
        if 56 < carf_dose <= 70:
            self.warnings.append("Carfilzomib input exceeds the configured reference value.")
        if time_horizon >= 300:
            self.warnings.append("Long horizon (>300 days) may magnify numerical stiffness; report solver sensitivity.")

        if errors:
            raise ValidationError(errors)

        # Optional experimental custom drug support (editor-only).
        custom_enabled = bool(cleaned.get("custom_drug_enabled"))
        if custom_enabled:
            from simulator.permissions import is_editor

            privileged = bool(
                self.user
                and getattr(self.user, "is_authenticated", False)
                and (getattr(self.user, "is_staff", False) or is_editor(self.user))
            )
            if not privileged:
                cleaned["custom_drug_enabled"] = False
                self.warnings.append(
                    "Custom drug disabled: editor permissions required for experimental agents."
                )
            else:
                name = (cleaned.get("custom_drug_name") or "").strip()
                if not name:
                    errors["custom_drug_name"] = ValidationError(
                        "Custom drug name is required when enabled."
                    )
                dose = float(cleaned.get("custom_drug_dose") or 0.0)
                if dose <= 0.0:
                    errors["custom_drug_dose"] = ValidationError(
                        "Custom drug dose must be > 0 when enabled."
                    )
                # Sanitize a stable key used in outputs/columns.
                key = slugify(name).replace("-", "_")
                key = "".join(ch for ch in key if ch.isalnum() or ch == "_")
                key = (key or "custom_drug")[:24]
                cleaned["custom_drug_key"] = f"custom_{key}" if not key.startswith("custom_") else key

                # Keep hard bounds conservative.
                half_life = float(cleaned.get("custom_pk_half_life") or 0.0)
                vd = float(cleaned.get("custom_pk_vd") or 0.0)
                emax = float(cleaned.get("custom_pd_emax") or 0.0)
                ec50 = float(cleaned.get("custom_pd_ec50") or 0.0)
                if half_life <= 0:
                    errors["custom_pk_half_life"] = ValidationError("Half-life must be > 0.")
                if vd <= 0:
                    errors["custom_pk_vd"] = ValidationError("Vd must be > 0.")
                if not (0.0 <= emax <= 1.0):
                    errors["custom_pd_emax"] = ValidationError("Emax must be within 0–1.")
                if ec50 <= 0:
                    errors["custom_pd_ec50"] = ValidationError("EC50 must be > 0.")

        if errors:
            raise ValidationError(errors)

        return cleaned
class SimulationAttemptForm(BootstrapValidationMixin, forms.ModelForm):
    """Form used to record a learner's research hypothesis."""

    selected_regimen = forms.ModelChoiceField(
        label="Exploratory regimen hypothesis",
        queryset=Regimen.objects.all(),
        required=False,
        help_text="Optional catalog reference for a synthetic learning comparison.",
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
        help_text="Confidence in this research hypothesis (%).",
        widget=forms.NumberInput(attrs={"class": "form-control", "min": "0", "max": "100"}),
    )
    notes = forms.CharField(
        widget=forms.Textarea(attrs={"rows": 4, "class": "form-control"}),
        required=False,
        label="Research reasoning",
    )

    class Meta:
        model = models.SimulationAttempt
        fields = ["selected_regimen", "predicted_response", "confidence", "notes"]

    def clean(self):
        cleaned = super().clean()
        selected_regimen = cleaned.get("selected_regimen")
        predicted_response = (cleaned.get("predicted_response") or "").strip()
        notes = (cleaned.get("notes") or "").strip()

        if not selected_regimen and not predicted_response and not notes:
            raise ValidationError(
                "Please select a regimen or expected response, or add a short note."
            )

        return cleaned

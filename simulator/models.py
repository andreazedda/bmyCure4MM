from __future__ import annotations

import json
import logging
from pathlib import Path

from django.conf import settings
from django.db import models
from django.utils import timezone

import numpy as np
import pandas as pd
from plotly import io as pio
from plotly.subplots import make_subplots
import plotly.graph_objects as go

IMWG_RESPONSE_CHOICES = [
    ("sCR", "Stringent CR"),
    ("CR", "Complete Response"),
    ("VGPR", "Very Good Partial Response"),
    ("PR", "Partial Response"),
    ("SD", "Stable Disease"),
    ("PD", "Progressive Disease"),
]

from .models_simulation import MathematicalModel
from .pharmaco import registry as pharmaco_registry
from .permissions import is_editor
from .twin import build_patient_twin

logger = logging.getLogger(__name__)

AUTO_SENTINELS = {"auto", "AUTO", "Auto", ""}
DEFAULT_PK_PARAMS = {
    "lenalidomide": {"half_life": 3.0, "Vd": 50.0},
    "bortezomib": {"half_life": 40.0, "Vd": 20.0},
    "daratumumab": {"half_life": 480.0, "Vd": 60.0},
}
DEFAULT_PD_PARAMS = {
    "lenalidomide": {"Emax": 0.8, "EC50": 10.0},
    "bortezomib": {"Emax": 0.9, "EC50": 0.2},
    "daratumumab": {"Emax": 0.95, "EC50": 50.0},
}


class Scenario(models.Model):
    """
    Clinical scenarios used to practice MM treatment decision making.
    
    Enhanced with comprehensive clinical parameters for realistic case scenarios:
    - Cytogenetics: del(17p), t(4;14), t(14;16), 1q21 gain, hyperdiploid, t(11;14)
    - Tumor biology: cell count, growth rate, carrying capacity
    - Patient characteristics: age, ECOG, Charlson comorbidity index
    - Laboratory values: renal function, albumin, β2M, LDH, hemoglobin, calcium
    - Risk stratification: R-ISS stage, patient archetype
    - Difficulty scoring: auto-calculated from clinical parameters
    """

    STAGE_CHOICES = [
        ("newly_diagnosed", "Newly Diagnosed"),
        ("relapsed_refractory", "Relapsed/Refractory"),
        ("maintenance", "Maintenance"),
        ("supportive", "Supportive Care"),
    ]

    RISS_CHOICES = [
        ("", "Not specified"),
        ("I", "R-ISS I (Low risk)"),
        ("II", "R-ISS II (Intermediate risk)"),
        ("III", "R-ISS III (High risk)"),
    ]

    ECOG_CHOICES = [
        (0, "0 - Fully active"),
        (1, "1 - Restricted in strenuous activity"),
        (2, "2 - Ambulatory, self-care"),
        (3, "3 - Limited self-care"),
        (4, "4 - Completely disabled"),
    ]

    ARCHETYPE_CHOICES = [
        ("", "Custom/Mixed"),
        ("nd_standard", "Newly Diagnosed - Standard Risk"),
        ("nd_high_risk", "Newly Diagnosed - High Risk"),
        ("frail_elderly", "Frail Elderly (Age >75)"),
        ("renal_impaired", "Renal Impairment (CrCl <60)"),
        ("rr_early", "Relapsed/Refractory - Early (1-2 lines)"),
        ("rr_late", "Relapsed/Refractory - Late (≥3 lines)"),
        ("aggressive", "Aggressive/Plasma Cell Leukemia"),
    ]

    DIFFICULTY_CHOICES = [
        ("", "Not calculated"),
        ("easy", "Easy (0-30)"),
        ("moderate", "Moderate (30-50)"),
        ("hard", "Hard (50-70)"),
        ("very_hard", "Very Hard (70-85)"),
        ("expert", "Expert (85-100)"),
    ]

    # Basic scenario fields
    title = models.CharField(max_length=200)
    clinical_stage = models.CharField(
        max_length=32,
        choices=STAGE_CHOICES,
        default="newly_diagnosed",
    )
    summary = models.TextField(
        help_text="Concise description of the case including age, stage and key features."
    )
    risk_stratification = models.CharField(
        max_length=128,
        blank=True,
        help_text="e.g., High-risk cytogenetics, R-ISS II.",
    )
    lab_snapshot = models.JSONField(
        blank=True,
        default=dict,
        help_text="Dictionary of key laboratory values shown to the learner.",
    )
    recommended_regimens = models.ManyToManyField(
        "clinic.Regimen",
        related_name="training_scenarios",
        blank=True,
    )
    guideline_notes = models.TextField(
        blank=True,
        help_text="Rationale, supporting evidence or pearls displayed after completion.",
    )
    expected_response = models.CharField(
        max_length=4,
        choices=IMWG_RESPONSE_CHOICES,
        blank=True,
        help_text="Target IMWG response. Optional; used for feedback.",
    )

    # Patient characteristics
    patient_age = models.IntegerField(
        null=True,
        blank=True,
        help_text="Patient age in years. Range: 18-120"
    )
    ecog_performance_status = models.IntegerField(
        null=True,
        blank=True,
        choices=ECOG_CHOICES,
        help_text="ECOG Performance Status (0-4)"
    )
    charlson_comorbidity_index = models.IntegerField(
        null=True,
        blank=True,
        help_text="Charlson Comorbidity Index. Range: 0-10. Low: 0-1, Moderate: 2-3, High: ≥4"
    )
    patient_archetype = models.CharField(
        max_length=32,
        blank=True,
        default="",
        choices=ARCHETYPE_CHOICES,
        help_text="Patient archetype from virtual patient generator"
    )

    # Cytogenetics
    del17p = models.BooleanField(
        default=False,
        help_text="TP53 deletion - high risk",
        verbose_name="del(17p)"
    )
    t_4_14 = models.BooleanField(
        default=False,
        help_text="Translocation (4;14) - high risk",
        verbose_name="t(4;14)"
    )
    t_14_16 = models.BooleanField(
        default=False,
        help_text="Translocation (14;16) - very high risk",
        verbose_name="t(14;16)"
    )
    gain_1q21 = models.BooleanField(
        default=False,
        help_text="1q21 gain/amplification - proliferation advantage",
        verbose_name="1q21 gain"
    )
    hyperdiploid = models.BooleanField(
        default=False,
        help_text="Hyperdiploid karyotype - standard risk",
        verbose_name="Hyperdiploid"
    )
    t_11_14 = models.BooleanField(
        default=False,
        help_text="Translocation (11;14) - standard risk, better prognosis",
        verbose_name="t(11;14)"
    )

    # Tumor biology
    tumor_cell_count = models.FloatField(
        null=True,
        blank=True,
        help_text="Tumor cell count (cells). Range: 1e6-1e12. Typical newly diagnosed: 1e10"
    )
    tumor_growth_rate = models.FloatField(
        null=True,
        blank=True,
        help_text="Tumor growth rate (per day). Range: 0.001-0.1. Typical: 0.01"
    )
    carrying_capacity = models.FloatField(
        null=True,
        blank=True,
        help_text="Maximum tumor burden (cells). Range: 1e11-1e13. Default: 1e12"
    )

    # Laboratory values
    creatinine_clearance = models.FloatField(
        null=True,
        blank=True,
        help_text="Creatinine clearance (mL/min). Normal: >60, Moderate renal impairment: 30-60, Severe: <30"
    )
    albumin = models.FloatField(
        null=True,
        blank=True,
        help_text="Serum albumin (g/dL). Normal: 3.5-5.0"
    )
    beta2_microglobulin = models.FloatField(
        null=True,
        blank=True,
        help_text="β2-microglobulin (mg/L). Normal: <2, Elevated: 2-5.5, Very high: >5.5",
        verbose_name="β2-microglobulin"
    )
    ldh = models.FloatField(
        null=True,
        blank=True,
        help_text="Lactate dehydrogenase (U/L). Normal: 140-280, Elevated: >280",
        verbose_name="LDH"
    )
    hemoglobin = models.FloatField(
        null=True,
        blank=True,
        help_text="Hemoglobin (g/dL). Normal: M 13-17, F 12-15. Anemia: <10"
    )
    calcium = models.FloatField(
        null=True,
        blank=True,
        help_text="Serum calcium (mg/dL). Normal: 8.5-10.2, Hypercalcemia: >10.5, Severe: >14"
    )

    # Risk stratification
    riss_stage = models.CharField(
        max_length=16,
        blank=True,
        default="",
        choices=RISS_CHOICES,
        help_text="Revised International Staging System stage",
        verbose_name="R-ISS Stage"
    )

    # Calculated fields
    difficulty_score = models.FloatField(
        null=True,
        blank=True,
        help_text="Calculated difficulty score (0-100). Auto-populated from difficulty scoring system."
    )
    difficulty_level = models.CharField(
        max_length=16,
        blank=True,
        default="",
        choices=DIFFICULTY_CHOICES,
        help_text="Difficulty level category"
    )

    # Metadata
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
    active = models.BooleanField(default=True)

    class Meta:
        ordering = ["title"]

    def __str__(self) -> str:
        return self.title

    def lab_items(self) -> list[tuple[str, str]]:
        snapshot = self.lab_snapshot or {}
        return list(snapshot.items())


class SimulationAttempt(models.Model):
    """Captures a clinician's decisions for a given scenario."""

    scenario = models.ForeignKey(
        Scenario,
        on_delete=models.CASCADE,
        related_name="attempts",
    )
    user = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    selected_regimen = models.ForeignKey(
        "clinic.Regimen",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    predicted_response = models.CharField(
        max_length=4,
        choices=IMWG_RESPONSE_CHOICES,
        blank=True,
    )
    confidence = models.PositiveSmallIntegerField(
        default=70,
        help_text="Confidence in chosen plan (0-100%).",
    )
    notes = models.TextField(blank=True)
    submitted = models.DateTimeField(auto_now_add=True)
    is_guideline_aligned = models.BooleanField(default=False)
    parameters = models.JSONField(default=dict, blank=True)
    results = models.JSONField(default=dict, blank=True)
    results_summary = models.JSONField(default=dict, blank=True)
    seed = models.IntegerField(null=True, blank=True)
    cohort_size = models.PositiveIntegerField(default=1)
    artifacts = models.JSONField(default=dict, blank=True)

    class Meta:
        ordering = ["-submitted"]

    def __str__(self) -> str:
        user_display = self.user.get_username() if self.user else "anonymous"
        return f"{self.scenario} attempt by {user_display} @ {self.submitted:%Y-%m-%d}"

    def run_model(self) -> dict:
        """Instantiate and execute the mathematical model, persisting outputs."""
        params = dict(self.parameters or {})
        resolved_params = dict(params)
        twin_payload: dict[str, float] = {}
        twin_assessment = None

        use_predlab = getattr(settings, "PREDLAB_V2", False)
        use_twin = resolved_params.get("use_twin", True)
        twin_mode = (resolved_params.get("twin_biology_mode") or "auto").strip().lower()
        if use_twin:
            assessment_id = resolved_params.get("twin_assessment_id") or resolved_params.get("assessment_id")
            if assessment_id:
                try:
                    from clinic.models import Assessment
                    from .access import get_accessible_assessment_by_id

                    if self.user and getattr(self.user, "is_authenticated", False):
                        is_privileged = self.user.is_staff or is_editor(self.user)
                        if is_privileged:
                            twin_assessment = Assessment.objects.filter(pk=assessment_id).first()
                        else:
                            twin_assessment = get_accessible_assessment_by_id(self.user, int(assessment_id))
                    else:
                        twin_assessment = None
                except Exception:
                    twin_assessment = None
            if twin_assessment:
                twin_payload = build_patient_twin(twin_assessment)
                if twin_mode == "auto":
                    resolved_params = self._merge_twin_parameters(resolved_params, twin_payload)

        try:
            solver_inputs = self._resolve_solver_inputs(resolved_params)
        except ValueError as exc:
            logger.warning(
                "solver-input-resolution-failed attempt_id=%s user_id=%s twin_mode=%s error=%s",
                self.pk,
                self.user_id,
                (resolved_params.get("twin_biology_mode") or "auto"),
                str(exc),
            )
            raise
        baseline_tumor = solver_inputs["baseline_tumor_cells"]
        baseline_healthy = solver_inputs["baseline_healthy_cells"]
        time_horizon = solver_inputs["time_horizon"]
        growth_rates = solver_inputs["growth_rates"]
        drug_doses = solver_inputs["drug_doses"]
        default_pk = {drug: values.copy() for drug, values in DEFAULT_PK_PARAMS.items()}
        default_pd = {drug: values.copy() for drug, values in DEFAULT_PD_PARAMS.items()}
        if use_predlab:
            pk_params, pd_params, dose_functions = pharmaco_registry.resolve(
                drug_doses,
                time_horizon,
                default_pk,
                default_pd,
            )
        else:
            pk_params = default_pk
            pd_params = default_pd
            dose_functions = {}
        interaction_strength = solver_inputs["interaction_strength"]
        interaction_matrix = np.eye(len(drug_doses), dtype=float) * interaction_strength
        immune_index = solver_inputs["immune_compromise_index"]
        carrying_t_override = solver_inputs["carrying_capacity_tumor"]
        carrying_h_override = solver_inputs["carrying_capacity_healthy"]

        model = MathematicalModel(
            baseline_tumor_cells=baseline_tumor,
            baseline_healthy_cells=baseline_healthy,
            drug_doses=drug_doses,
            pk_params=pk_params,
            pd_params=pd_params,
            growth_rates=growth_rates,
            interaction_matrix=interaction_matrix,
            time_span=(0.0, time_horizon),
            carrying_capacity_tumor=carrying_t_override,
            carrying_capacity_healthy=carrying_h_override,
            immune_compromise_index=immune_index,
            dose_functions=dose_functions,
        )
        trajectory = model.simulate()

        tumor_series = trajectory["tumor_cells"]
        healthy_series = trajectory["healthy_cells"]
        time_series = trajectory["time"]
        tumor_start = tumor_series.iloc[0]
        tumor_end = tumor_series.iloc[-1]
        healthy_start = healthy_series.iloc[0]
        healthy_end = healthy_series.iloc[-1]
        tumor_reduction = float(1.0 - tumor_end / max(tumor_start, 1e-9))
        healthy_loss = float(1.0 - healthy_end / max(healthy_start, 1e-9))
        nadir_idx = tumor_series.idxmin()
        post_nadir = trajectory.loc[nadir_idx + 1 :] if nadir_idx + 1 < len(trajectory) else pd.DataFrame()
        recurrence_threshold = 0.5 * tumor_start
        time_to_recurrence = None
        if not post_nadir.empty:
            recurrence = post_nadir[post_nadir["tumor_cells"] > recurrence_threshold]
            if not recurrence.empty:
                time_to_recurrence = float(recurrence["time"].iloc[0])
        auc_drug = {}
        for drug in drug_doses.keys():
            series = trajectory[f"{drug}_concentration"]
            auc_drug[drug] = float(np.trapz(series, time_series))
        effective = bool(tumor_reduction > 0.9 and healthy_loss < 0.2)
        summary = {
            "tumor_reduction": tumor_reduction,
            "healthy_loss": healthy_loss,
            "time_to_recurrence": time_to_recurrence,
            "auc": auc_drug,
            "effective": effective,
        }

        media_root = Path(settings.MEDIA_ROOT)
        plot_dir = media_root / "sim_plots"
        plot_dir.mkdir(parents=True, exist_ok=True)
        data_dir = media_root / "sim_data"
        data_dir.mkdir(parents=True, exist_ok=True)

        csv_filename = f"attempt_{self.pk}.csv"
        csv_path = data_dir / csv_filename
        trajectory.to_csv(csv_path, index=False)

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.6, 0.4])
        fig.add_trace(
            go.Scatter(x=time_series, y=tumor_series, name="Tumor burden", mode="lines"),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=time_series, y=healthy_series, name="Healthy plasma cells", mode="lines"),
            row=1,
            col=1,
        )
        for drug in drug_doses.keys():
            fig.add_trace(
                go.Scatter(
                    x=time_series,
                    y=trajectory[f"{drug}_concentration"],
                    name=f"{drug.title()} concentration",
                    mode="lines",
                ),
                row=2,
                col=1,
            )
        fig.update_layout(
            title=f"Simulation Outcome – {self.scenario.title}",
            xaxis_title="Time (days)",
            yaxis_title="Cell count",
            legend={"orientation": "h"},
            template="plotly_white",
        )
        plot_html = pio.to_html(fig, full_html=False, include_plotlyjs="cdn")
        plot_filename = f"attempt_{self.pk}.html"
        plot_path = plot_dir / plot_filename
        plot_path.write_text(plot_html, encoding="utf-8")

        twin_filename = None
        if twin_payload:
            twin_filename = f"attempt_{self.pk}_twin.json"
            twin_path = data_dir / twin_filename
            serialized = {
                **twin_payload,
                "assessment_id": twin_assessment.pk if twin_assessment else None,
            }
            twin_path.write_text(json.dumps(serialized, indent=2), encoding="utf-8")

        media_url = settings.MEDIA_URL.rstrip("/")
        result_payload = {
            "csv": f"{media_url}/sim_data/{csv_filename}",
            "plot": f"{media_url}/sim_plots/{plot_filename}",
            "generated_at": timezone.now().isoformat(),
        }
        if twin_filename:
            result_payload["twin_params.json"] = f"{media_url}/sim_data/{twin_filename}"

        self.results = result_payload
        self.results_summary = summary
        artifacts = {key: value for key, value in result_payload.items() if key != "generated_at"}
        if self.seed is not None:
            artifacts["seed"] = self.seed
        self.artifacts = artifacts
        self.save(update_fields=["results", "results_summary", "artifacts"])
        return summary

    @classmethod
    def _resolve_solver_inputs(cls, resolved_params: dict) -> dict:
        """Centralized conversion layer before the solver.

        Contract: no string (including 'auto') may reach MathematicalModel.
        """
        baseline_tumor = cls._float_or_strict(
            resolved_params.get("baseline_tumor_cells"), 1.0e9, "baseline_tumor_cells"
        )
        baseline_healthy = cls._float_or_strict(
            resolved_params.get("baseline_healthy_cells"), 5.0e11, "baseline_healthy_cells"
        )
        time_horizon = cls._float_or_strict(resolved_params.get("time_horizon"), 90.0, "time_horizon")
        growth_rates = {
            "tumor": cls._float_or_strict(
                resolved_params.get("tumor_growth_rate"), 0.023, "growth_rates.tumor"
            ),
            "healthy": cls._float_or_strict(
                resolved_params.get("healthy_growth_rate"), 0.015, "growth_rates.healthy"
            ),
        }
        drug_doses = {
            "lenalidomide": cls._float_or_strict(
                resolved_params.get("lenalidomide_dose"), 25.0, "drug_doses.lenalidomide"
            ),
            "bortezomib": cls._float_or_strict(
                resolved_params.get("bortezomib_dose"), 1.3, "drug_doses.bortezomib"
            ),
            "daratumumab": cls._float_or_strict(
                resolved_params.get("daratumumab_dose"), 16.0, "drug_doses.daratumumab"
            ),
        }
        interaction_strength = cls._float_or_strict(
            resolved_params.get("interaction_strength"), 0.05, "interaction_strength"
        )
        immune_index = cls._float_or_strict(
            resolved_params.get("immune_compromise_index"), 1.0, "immune_compromise_index"
        )
        carrying_t_override = cls._float_or_strict(
            resolved_params.get("carrying_capacity_tumor"), None, "carrying_capacity_tumor"
        )
        carrying_h_override = cls._float_or_strict(
            resolved_params.get("carrying_capacity_healthy"), None, "carrying_capacity_healthy"
        )

        return {
            "baseline_tumor_cells": baseline_tumor,
            "baseline_healthy_cells": baseline_healthy,
            "time_horizon": time_horizon,
            "growth_rates": growth_rates,
            "drug_doses": drug_doses,
            "interaction_strength": interaction_strength,
            "immune_compromise_index": immune_index,
            "carrying_capacity_tumor": carrying_t_override,
            "carrying_capacity_healthy": carrying_h_override,
        }

    @staticmethod
    def _float_or_strict(value, default, name: str):
        if value in AUTO_SENTINELS or value is None:
            return default
        if isinstance(value, str) and value.strip().lower() in {"auto", "", "none", "null"}:
            return default
        try:
            return float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid solver input for {name}: {value!r}") from exc

    @staticmethod
    def _merge_twin_parameters(params: dict, twin_payload: dict) -> dict:
        merged = dict(params)
        allowed = {
            "tumor_growth_rate",
            "healthy_growth_rate",
            "carrying_capacity_tumor",
            "carrying_capacity_healthy",
            "immune_compromise_index",
        }

        def _is_auto(value) -> bool:
            if value is None:
                return True
            if value in AUTO_SENTINELS:
                return True
            if isinstance(value, str) and value.strip().lower() in {"auto", "", "none", "null"}:
                return True
            return False

        for key in allowed:
            twin_value = twin_payload.get(key)
            if twin_value is None:
                continue
            current = merged.get(key)
            if key not in merged or _is_auto(current):
                merged[key] = twin_value
        return merged

    @staticmethod
    def _float_or(value, default):
        if value in AUTO_SENTINELS or value is None:
            return default
        if isinstance(value, str) and value.strip().lower() in {"auto", "", "none", "null"}:
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

from __future__ import annotations

import json
from pathlib import Path

from django.conf import settings
from django.db import models
from django.utils import timezone

import numpy as np
import pandas as pd
from plotly import io as pio
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from clinic.models import Assessment, Regimen

from .models_simulation import MathematicalModel
from .pharmaco import registry as pharmaco_registry
from .twin import build_patient_twin

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
    """Clinical scenarios used to practice MM treatment decision making."""

    STAGE_CHOICES = [
        ("newly_diagnosed", "Newly Diagnosed"),
        ("relapsed_refractory", "Relapsed/Refractory"),
        ("maintenance", "Maintenance"),
        ("supportive", "Supportive Care"),
    ]

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
        Regimen,
        related_name="training_scenarios",
        blank=True,
    )
    guideline_notes = models.TextField(
        blank=True,
        help_text="Rationale, supporting evidence or pearls displayed after completion.",
    )
    expected_response = models.CharField(
        max_length=4,
        choices=Assessment.RESPONSE_CHOICES,
        blank=True,
        help_text="Target IMWG response. Optional; used for feedback.",
    )
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
        Regimen,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    predicted_response = models.CharField(
        max_length=4,
        choices=Assessment.RESPONSE_CHOICES,
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
        if use_predlab and resolved_params.get("use_twin", True):
            assessment_id = resolved_params.get("twin_assessment_id") or resolved_params.get("assessment_id")
            if assessment_id:
                try:
                    twin_assessment = Assessment.objects.get(pk=assessment_id)
                except Assessment.DoesNotExist:
                    twin_assessment = None
            if twin_assessment:
                twin_payload = build_patient_twin(twin_assessment)
                resolved_params = self._merge_twin_parameters(resolved_params, twin_payload)

        baseline_tumor = self._float_or(resolved_params.get("baseline_tumor_cells"), 1.0e9)
        baseline_healthy = self._float_or(resolved_params.get("baseline_healthy_cells"), 5.0e11)
        time_horizon = self._float_or(resolved_params.get("time_horizon"), 90.0)
        growth_rates = {
            "tumor": self._float_or(resolved_params.get("tumor_growth_rate"), 0.023),
            "healthy": self._float_or(resolved_params.get("healthy_growth_rate"), 0.015),
        }
        drug_doses = {
            "lenalidomide": self._float_or(resolved_params.get("lenalidomide_dose"), 25.0),
            "bortezomib": self._float_or(resolved_params.get("bortezomib_dose"), 1.3),
            "daratumumab": self._float_or(resolved_params.get("daratumumab_dose"), 16.0),
        }
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
        interaction_strength = self._float_or(resolved_params.get("interaction_strength"), 0.05)
        interaction_matrix = np.eye(len(drug_doses), dtype=float) * interaction_strength
        immune_index = self._float_or(resolved_params.get("immune_compromise_index"), 1.0)
        carrying_t_override = self._float_or(resolved_params.get("carrying_capacity_tumor"), None)
        carrying_h_override = self._float_or(resolved_params.get("carrying_capacity_healthy"), None)

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
        self.artifacts = {key: value for key, value in result_payload.items() if key != "generated_at"}
        self.save(update_fields=["results", "results_summary", "artifacts"])
        return summary

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
        for key in allowed:
            twin_value = twin_payload.get(key)
            if twin_value is None:
                continue
            current = merged.get(key)
            if key not in merged or current in AUTO_SENTINELS or current is None:
                merged[key] = twin_value
        return merged

    @staticmethod
    def _float_or(value, default):
        if value in AUTO_SENTINELS or value is None:
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

from __future__ import annotations

from typing import Any

from .exposure_bridge import exposure_profile_from_payload
from .models import LongitudinalLabResult
from .toxicity_model import summarize_liver_toxicity, summarize_neutropenia_history


DEFAULT_UTILITY_V2_WEIGHTS = {
    "lambda_liver": 0.5,
    "lambda_neutropenia": 0.5,
}

DEFAULT_TOXICITY_COEFFICIENTS = {
    "alpha_ast": 0.35,
    "alpha_alt": 0.35,
    "beta_liver_steatosis": 0.15,
    "gamma_ast_prior_peak": 0.15,
    "gamma_alt_prior_peak": 0.15,
    "alpha_neu": 0.45,
    "beta_prior_neutropenia": 0.25,
}

HEURISTIC_DAILY_DOSE_REFERENCE_MG = 10.0
HEURISTIC_CUMULATIVE_DOSE_REFERENCE_MG = 100.0


def compute_toxicity_dynamics(patient, exposure_profile_payload: dict[str, Any] | None) -> dict[str, Any]:
    profile = _resolve_exposure_profile(exposure_profile_payload)
    if profile is None or not profile.daily_administered_dose_mg:
        return {
            "toxicity_model_status": "unavailable",
            "signal_units": "normalized_risk_signal_0_1",
            "predicted_ast_signal": None,
            "predicted_alt_signal": None,
            "predicted_neu_signal": None,
            "liver_toxicity_signal_0_1": None,
            "neutropenia_signal_0_1": None,
            "toxicity_risk_signal": None,
            "ast_signal_series": [],
            "alt_signal_series": [],
            "liver_signal_series": [],
            "neutropenia_signal_series": [],
            "toxicity_risk_series": [],
            "coefficient_diagnostics": {
                "n_observations_ast": 0,
                "n_observations_alt": 0,
                "n_observations_neu": 0,
                "fitted_or_default": "insufficient_data_default",
                "parameter_status": {
                    "ast": "insufficient_data_default",
                    "alt": "insufficient_data_default",
                    "neu": "insufficient_data_default",
                },
                "coefficient_values": dict(DEFAULT_TOXICITY_COEFFICIENTS),
                "uncertainty_status": "unavailable",
            },
            "limitation": "toxicity signals are prototype risk signals, not validated predictions.",
            "missing_reason": "Exposure metadata unavailable; regenerate run to compute exposure profile.",
        }

    ast_entries = _lab_values(patient, LongitudinalLabResult.ANALYTE_AST)
    alt_entries = _lab_values(patient, LongitudinalLabResult.ANALYTE_ALT)
    neu_entries = _lab_values(patient, LongitudinalLabResult.ANALYTE_NEU)

    diagnostics = estimate_toxicity_coefficients(
        ast_values=ast_entries,
        alt_values=alt_entries,
        neu_values=neu_entries,
    )
    coeffs = diagnostics["coefficient_values"]

    liver_summary = summarize_liver_toxicity(patient)
    neut_summary = summarize_neutropenia_history(patient)
    daily = [float(value) for value in profile.daily_administered_dose_mg]
    cumulative = [float(value) for value in profile.cumulative_dose_mg]
    daily_norm = [_clip01(float(value) / HEURISTIC_DAILY_DOSE_REFERENCE_MG) for value in daily]
    cumulative_norm = [_clip01(float(value) / HEURISTIC_CUMULATIVE_DOSE_REFERENCE_MG) for value in cumulative]

    steatosis_flag = 1.0 if liver_summary.get("hepatic_steatosis") is True else 0.0
    ast_prior_peak = _peak_normalized(ast_entries)
    alt_prior_peak = _peak_normalized(alt_entries)
    prior_neutropenia = 1.0 if neut_summary.get("absolute_neutropenia_event") or _has_low_neutrophil(neu_entries) else 0.0

    ast_baseline = min(0.05 + ast_prior_peak * 0.1, 0.35)
    alt_baseline = min(0.05 + alt_prior_peak * 0.1, 0.35)
    neut_baseline = min(0.05 + prior_neutropenia * 0.1, 0.35)

    ast_signal_series = []
    alt_signal_series = []
    liver_signal_series = []
    neutropenia_signal_series = []
    toxicity_risk_series = []
    for day in profile.time_grid_days:
        ast_signal = _clip01(
            ast_baseline
            + coeffs["alpha_ast"] * daily_norm[day]
            + coeffs["beta_liver_steatosis"] * steatosis_flag
            + coeffs["gamma_ast_prior_peak"] * ast_prior_peak
        )
        alt_signal = _clip01(
            alt_baseline
            + coeffs["alpha_alt"] * daily_norm[day]
            + coeffs["beta_liver_steatosis"] * steatosis_flag
            + coeffs["gamma_alt_prior_peak"] * alt_prior_peak
        )
        liver_signal = _clip01(max(ast_signal, alt_signal))
        neut_signal = _clip01(
            neut_baseline
            + coeffs["alpha_neu"] * cumulative_norm[day]
            + coeffs["beta_prior_neutropenia"] * prior_neutropenia
        )
        risk_signal = _clip01(max(liver_signal, neut_signal))
        ast_signal_series.append({"day": int(day), "value": ast_signal})
        alt_signal_series.append({"day": int(day), "value": alt_signal})
        liver_signal_series.append({"day": int(day), "value": liver_signal})
        neutropenia_signal_series.append({"day": int(day), "value": neut_signal})
        toxicity_risk_series.append({"day": int(day), "value": risk_signal})

    return {
        "toxicity_model_status": "semi_mechanistic_prototype",
        "signal_units": "normalized_risk_signal_0_1",
        "predicted_ast_signal": ast_signal_series[-1]["value"],
        "predicted_alt_signal": alt_signal_series[-1]["value"],
        "predicted_neu_signal": neutropenia_signal_series[-1]["value"],
        "liver_toxicity_signal_0_1": liver_signal_series[-1]["value"],
        "neutropenia_signal_0_1": neutropenia_signal_series[-1]["value"],
        "toxicity_risk_signal": toxicity_risk_series[-1]["value"],
        "ast_signal_series": ast_signal_series,
        "alt_signal_series": alt_signal_series,
        "liver_signal_series": liver_signal_series,
        "neutropenia_signal_series": neutropenia_signal_series,
        "toxicity_risk_series": toxicity_risk_series,
        "coefficient_diagnostics": diagnostics,
        "limitation": "toxicity signals are prototype risk signals, not validated predictions.",
    }


def estimate_toxicity_coefficients(*, ast_values: list[float], alt_values: list[float], neu_values: list[float]) -> dict[str, Any]:
    ast_coeff, ast_status = _coefficient_from_values(ast_values, DEFAULT_TOXICITY_COEFFICIENTS["alpha_ast"])
    alt_coeff, alt_status = _coefficient_from_values(alt_values, DEFAULT_TOXICITY_COEFFICIENTS["alpha_alt"])
    neu_coeff, neu_status = _coefficient_from_values(neu_values, DEFAULT_TOXICITY_COEFFICIENTS["alpha_neu"])
    fitted_statuses = {ast_status, alt_status, neu_status}
    fitted_or_default = "fitted_estimate" if "fitted_estimate" in fitted_statuses else "insufficient_data_default"
    uncertainty_status = "limited_single_patient_estimate" if fitted_or_default == "fitted_estimate" else "insufficient_data_default"
    return {
        "n_observations_ast": len(ast_values),
        "n_observations_alt": len(alt_values),
        "n_observations_neu": len(neu_values),
        "fitted_or_default": fitted_or_default,
        "parameter_status": {
            "ast": ast_status,
            "alt": alt_status,
            "neu": neu_status,
        },
        "coefficient_values": {
            **DEFAULT_TOXICITY_COEFFICIENTS,
            "alpha_ast": ast_coeff,
            "alpha_alt": alt_coeff,
            "alpha_neu": neu_coeff,
        },
        "uncertainty_status": uncertainty_status,
    }


def _coefficient_from_values(values: list[float], default: float) -> tuple[float, str]:
    if len(values) < 3:
        return float(default), "insufficient_data_default"
    normalized = _normalize(values)
    span = max(normalized) - min(normalized)
    estimate = _clip(float(default) + span * 0.4, 0.05, 0.95)
    return estimate, "fitted_estimate"


def _lab_values(patient, analyte: str) -> list[float]:
    return [
        float(item.value)
        for item in LongitudinalLabResult.objects.filter(patient=patient, analyte=analyte).exclude(value__isnull=True).order_by("date", "id")
    ]


def _normalize(values: list[float]) -> list[float]:
    if not values:
        return []
    lo = min(values)
    hi = max(values)
    if abs(hi - lo) < 1.0e-9:
        return [0.5 for _ in values]
    return [(float(value) - lo) / (hi - lo) for value in values]


def _peak_normalized(values: list[float]) -> float:
    normalized = _normalize(values)
    return max(normalized) if normalized else 0.0


def _has_low_neutrophil(values: list[float]) -> bool:
    return any(float(value) < 1.5 for value in values)


def _resolve_exposure_profile(payload: dict[str, Any] | None):
    if not payload:
        return None
    if "daily_administered_dose_mg" in payload:
        return exposure_profile_from_payload(payload)

    profiles = []
    for value in payload.values():
        if isinstance(value, dict) and "daily_administered_dose_mg" in value:
            profile = exposure_profile_from_payload(value)
            if profile is not None:
                profiles.append(profile)
    if not profiles:
        return None

    horizon = max(profile.horizon_days for profile in profiles)
    daily = [0.0 for _ in range(horizon)]
    interruptions = set()
    for profile in profiles:
        for index, value in enumerate(profile.daily_administered_dose_mg):
            daily[index] += float(value)
        interruptions.update(int(day) for day in profile.interruption_days)
    cumulative = []
    running_total = 0.0
    for value in daily:
        running_total += float(value)
        cumulative.append(running_total)
    peak = max(daily) if daily else 0.0
    total = cumulative[-1] if cumulative else 0.0
    return exposure_profile_from_payload(
        {
            "horizon_days": horizon,
            "time_grid_days": list(range(horizon)),
            "drug": "combined_exposure",
            "dose_events": [],
            "daily_administered_dose_mg": daily,
            "cumulative_dose_mg": cumulative,
            "average_daily_dose_mg": total / max(horizon, 1),
            "peak_administered_dose_mg": peak,
            "dosing_frequency_days": None,
            "interruption_days": sorted(interruptions),
            "schedule_type": "combined",
            "schedule_label": "combined_exposure",
            "exposure_resolution_status": "time_resolved_profile_available",
            "total_cumulative_dose_mg": total,
            "interruption_fraction": len(interruptions) / max(horizon, 1),
            "exposure_profile_hash": "combined_exposure",
        }
    )


def _clip01(value: float) -> float:
    return _clip(value, 0.0, 1.0)


def _clip(value: float, lower: float, upper: float) -> float:
    return max(lower, min(float(value), upper))

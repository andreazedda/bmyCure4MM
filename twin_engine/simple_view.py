from __future__ import annotations

from datetime import date
from typing import Any, Iterable

from django.urls import reverse

from clinic.models import Patient, PatientTherapy

from .cockpit import (
    build_assessment_recommendations,
    build_scenario_rows,
    build_validation_panel,
    list_latest_completed_runs_by_label,
)
from .data_explainability import build_data_explainability_context
from .developer_checks import detect_schedule_collapse
from .models import AdverseEvent, CounterfactualRun, LongitudinalLabResult, TherapyInterruption
from .state_model import get_current_twin_state


GLOBAL_BANNER = "Research simulation only. Not a treatment recommendation. Causal effect not identified."


def build_simple_patient_story(patient: Patient, *, include_developer_links: bool = False) -> dict[str, Any]:
    current_state = get_current_twin_state(patient)
    has_runs = patient.counterfactual_runs.filter(status=CounterfactualRun.STATUS_COMPLETED).exists()
    collapse_warnings = detect_schedule_collapse(patient) if has_runs else []
    latest_runs_by_label = list_latest_completed_runs_by_label(patient)
    scenario_rows = build_scenario_rows(patient, latest_runs_by_label, collapse_warnings)
    validation_panel = build_validation_panel(patient, current_state, scenario_rows)
    recommendation = build_assessment_recommendations(patient)
    explainability = build_data_explainability_context(
        patient,
        current_state=current_state,
        scenario_rows=scenario_rows,
        validation_panel=validation_panel,
    )
    data_blocks = _build_data_blocks(
        patient,
        current_state=current_state,
        scenario_rows=scenario_rows,
        include_developer_links=include_developer_links,
    )
    dictionary_lookup = explainability["dictionary_lookup"]
    for block in data_blocks:
        meta = dictionary_lookup.get(block["title"])
        if meta:
            block["data_source"] = meta["data_source"]
            block["model_use"] = meta["model_use"]
            block["missing_data_message"] = meta["missingness_explanation"]
    next_action = _determine_next_action(patient, current_state, scenario_rows, validation_panel, recommendation)
    scenario_story = _build_scenario_story(scenario_rows, validation_panel)
    scenario_explanations = explainability["scenario_explanations"]
    for scenario in scenario_story:
        scenario.update(scenario_explanations.get(scenario["scenario_label"], {}))
    story = {
        "patient_summary": {
            "pseudonym": f"Research Patient {patient.id}",
            "record_count_summary": _record_count_summary(patient, scenario_rows),
            "date_range": _patient_date_range(patient),
            "model_readiness_status": _model_readiness_status(patient, current_state, scenario_rows, validation_panel),
            "next_action": next_action,
        },
        "data_blocks": data_blocks,
        "data_dictionary_blocks": explainability["groups"],
        "classification_legend": explainability["classification_legend"],
        "help_box_steps": explainability["help_box_steps"],
        "model_story": _build_model_story(current_state, scenario_rows, validation_panel),
        "scenario_story": scenario_story,
        "conclusion_cards": _build_conclusion_cards(data_blocks, scenario_rows, validation_panel, next_action),
        "one_minute_summary": _build_one_minute_summary(patient, current_state, scenario_rows, validation_panel, data_blocks),
        "steps": _build_steps(patient, current_state, scenario_rows, validation_panel),
        "input_table": explainability["model_input_map"],
        "lineage_nodes": _build_lineage_nodes(patient, current_state, scenario_rows, validation_panel),
        "allowed_conclusions": [
            "The model can compare simulated trajectories under different schedules.",
            "Some schedules produce different exposure profiles.",
            "Some outputs are more uncertain than others.",
        ],
        "forbidden_conclusions": [
            "This is the best treatment.",
            "This proves clinical superiority.",
            "This estimates a causal effect.",
        ],
        "navigation": {
            "patient_page_url": reverse("clinic:patient_detail", args=[patient.id]),
            "simple_view_url": reverse("twin_engine:simple_research_view", args=[patient.id]),
            "cockpit_url": reverse("twin_engine:research_cockpit", args=[patient.id]),
            "developer_console_url": reverse("twin_engine:developer_console") + f"?patient_id={patient.id}",
            "glossary_url": reverse("twin_engine:research_glossary"),
        },
        "global_banner": GLOBAL_BANNER,
    }
    return story


def _build_data_blocks(
    patient: Patient,
    *,
    current_state,
    scenario_rows: list[dict[str, Any]],
    include_developer_links: bool,
) -> list[dict[str, Any]]:
    disease_counts = _field_counts(
        patient,
        [
            ("FLC ratio", "flc_ratio", LongitudinalLabResult.ANALYTE_FLC_RATIO),
            ("M-protein", "m_protein_g_dl", LongitudinalLabResult.ANALYTE_M_PROTEIN),
            ("kappa FLC", None, LongitudinalLabResult.ANALYTE_KAPPA_FLC),
            ("lambda FLC", None, LongitudinalLabResult.ANALYTE_LAMBDA_FLC),
        ],
    )
    blood_counts = _field_counts(
        patient,
        [
            ("hemoglobin", "hemoglobin_g_dl", LongitudinalLabResult.ANALYTE_HB),
            ("WBC", None, LongitudinalLabResult.ANALYTE_WBC),
            ("NEU", None, LongitudinalLabResult.ANALYTE_NEU),
            ("PLT", None, LongitudinalLabResult.ANALYTE_PLT),
        ],
    )
    liver_counts = _field_counts(
        patient,
        [
            ("AST", None, LongitudinalLabResult.ANALYTE_AST),
            ("ALT", None, LongitudinalLabResult.ANALYTE_ALT),
        ],
    )
    liver_counts["hepatic steatosis"] = patient.adverse_events.filter(event_type=AdverseEvent.TYPE_HEPATIC_STEATOSIS).count()
    liver_counts["hypertransaminasemia events"] = patient.adverse_events.filter(event_type=AdverseEvent.TYPE_HYPERTRANSAMINASEMIA).count()

    therapies = list(patient.therapies.order_by("start_date"))
    interruptions = list(patient.therapy_interruptions.order_by("start_date"))
    adverse_events = list(patient.adverse_events.order_by("date"))
    blocks = [
        _make_block(
            title="Disease markers",
            plain_meaning="These values are used as disease activity signals.",
            fields_used=["FLC ratio", "M-protein", "kappa FLC", "lambda FLC"],
            available_count=sum(disease_counts.values()),
            date_range=_date_range_for_labs(patient, [
                LongitudinalLabResult.ANALYTE_FLC_RATIO,
                LongitudinalLabResult.ANALYTE_M_PROTEIN,
                LongitudinalLabResult.ANALYTE_KAPPA_FLC,
                LongitudinalLabResult.ANALYTE_LAMBDA_FLC,
            ], assessment_fields=["flc_ratio", "m_protein_g_dl"]),
            status=_status_from_field_counts(disease_counts, enough_threshold=2),
            used_by_model=True,
            why_it_matters="The model uses these values to estimate whether the disease-marker signal is stable, increasing, or decreasing.",
            missing_data_message=_missing_message(disease_counts, "Add at least M-protein or free light chain values to support disease-signal interpretation."),
            next_action=_block_action(
                title="Add or verify disease-marker records",
                detail="Confirm structured disease-marker rows before interpreting simulated disease-control signals.",
                href=reverse("clinic:assessment_new", args=[patient.id]),
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#data-availability"),
            limit="These are biomarkers, not a complete disease state.",
        ),
        _make_block(
            title="Blood / marrow markers",
            plain_meaning="These values help describe blood reserve and toxicity-related vulnerability.",
            fields_used=["hemoglobin", "WBC", "NEU", "PLT"],
            available_count=sum(blood_counts.values()),
            date_range=_date_range_for_labs(patient, [
                LongitudinalLabResult.ANALYTE_HB,
                LongitudinalLabResult.ANALYTE_WBC,
                LongitudinalLabResult.ANALYTE_NEU,
                LongitudinalLabResult.ANALYTE_PLT,
            ], assessment_fields=["hemoglobin_g_dl"]),
            status=_status_from_field_counts(blood_counts, enough_threshold=2),
            used_by_model=True,
            why_it_matters="These values help describe blood reserve and toxicity-related vulnerability.",
            missing_data_message=_missing_message(blood_counts, "Structured hemoglobin, white-cell, neutrophil, or platelet rows are incomplete."),
            next_action=_block_action(
                title="Add blood-count records",
                detail="Add structured blood-count fields before treating blood-vulnerability signals as informative.",
                href=reverse("clinic:assessment_new", args=[patient.id]),
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#data-availability"),
            limit="The current model does not fully reconstruct marrow biology.",
        ),
        _make_block(
            title="Liver / toxicity markers",
            plain_meaning="These values constrain how treatment schedules are interpreted.",
            fields_used=["AST", "ALT", "hepatic steatosis", "hypertransaminasemia events"],
            available_count=sum(liver_counts.values()),
            date_range=_merge_ranges(
                _date_range_for_labs(patient, [LongitudinalLabResult.ANALYTE_AST, LongitudinalLabResult.ANALYTE_ALT]),
                _date_range_for_events(adverse_events, {AdverseEvent.TYPE_HEPATIC_STEATOSIS, AdverseEvent.TYPE_HYPERTRANSAMINASEMIA}),
            ),
            status=_status_from_field_counts(liver_counts, enough_threshold=2),
            used_by_model=True,
            why_it_matters="These values constrain how treatment schedules are interpreted.",
            missing_data_message=_missing_message(liver_counts, "Structured AST, ALT, or liver-event context is incomplete."),
            next_action=_block_action(
                title="Review liver and toxicity context",
                detail="Record AST, ALT, and liver-event context before over-interpreting toxicity signals.",
                href=reverse("clinic:assessment_new", args=[patient.id]),
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#toxicity-constraints"),
            limit="The toxicity layer produces prototype risk signals, not validated AST/ALT predictions.",
        ),
        _make_block(
            title="Treatment schedules",
            plain_meaning="Recorded treatment rows are converted into dated dose patterns.",
            fields_used=["drug", "dose", "schedule", "start date", "stop date", "interruptions"],
            available_count=len(therapies),
            date_range=_date_range_for_therapies(therapies),
            status="enough" if therapies else "missing",
            used_by_model=True,
            why_it_matters="The exposure bridge converts treatment schedules into daily dose profiles.",
            missing_data_message="No structured treatment schedule is currently available." if not therapies else "Treatment timing exists, but schedule details may still be partial.",
            next_action=_block_action(
                title="Review treatment schedule details",
                detail="Check dose timing, stop dates, and cycle notes before comparing scenarios.",
                href=reverse("clinic:patient_detail", args=[patient.id]),
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#what-if-scenarios"),
            limit="Dose timing is modeled, but clinical efficacy is not validated.",
        ),
        _make_block(
            title="Interruptions",
            plain_meaning="These rows describe pauses, holds, or changes in exposure over time.",
            fields_used=["drug", "reason", "start date", "stop date"],
            available_count=len(interruptions),
            date_range=_date_range_for_interruptions(interruptions),
            status="enough" if interruptions else "missing",
            used_by_model=True,
            why_it_matters="These rows explain whether the recorded schedule included pauses that affect interpretation.",
            missing_data_message="No structured interruption rows are recorded." if not interruptions else "Interruption context is available.",
            next_action=_block_action(
                title="Check for treatment holds",
                detail="Add interruption rows if the treatment history included pauses or dose holds.",
                href=reverse("clinic:patient_detail", args=[patient.id]),
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:developer_console") + f"?patient_id={patient.id}"),
            limit="An interruption record gives exposure context; it does not prove why an outcome changed.",
        ),
        _make_block(
            title="Adverse events",
            plain_meaning="These rows describe observed safety events that may matter for interpretation.",
            fields_used=["event type", "date", "grade", "action taken"],
            available_count=len(adverse_events),
            date_range=_date_range_for_events(adverse_events),
            status="enough" if adverse_events else "missing",
            used_by_model=True,
            why_it_matters="Observed safety events provide context for toxicity interpretation.",
            missing_data_message="No structured adverse-event rows are recorded." if not adverse_events else "Observed safety context is available.",
            next_action=_block_action(
                title="Review observed safety events",
                detail="Add or verify adverse-event rows before interpreting toxicity summaries as clinically complete.",
                href=reverse("clinic:patient_detail", args=[patient.id]),
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#toxicity-constraints"),
            limit="Observed safety events are descriptive context, not future toxicity predictions.",
        ),
        _make_block(
            title="Twin state",
            plain_meaning="This is the mathematical starting state used before running a simulation.",
            fields_used=["state date", "initialization method", "model version"],
            available_count=patient.twin_states.count(),
            date_range=_date_range_for_states(patient),
            status="enough" if current_state is not None else ("partial" if patient.twin_states.exists() else "missing"),
            used_by_model=True,
            why_it_matters="The simulation needs a mathematical starting state before it can generate trajectories.",
            missing_data_message="No current mathematical starting state is available." if current_state is None else "A current mathematical starting state is available.",
            next_action=_block_action(
                title="Create or review the mathematical starting state",
                detail="Initialize from the strongest available assessment before interpreting scenario outputs.",
                href=reverse("twin_engine:research_cockpit", args=[patient.id]) + "#twin-state",
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#twin-state"),
            limit="A mathematical starting state is a model setup step, not a validated digital copy of the patient.",
        ),
        _make_block(
            title="Simulated scenarios",
            plain_meaning="The scenarios compare alternative treatment schedules under the current mechanistic model.",
            fields_used=["scenario label", "schedule", "exposure pattern", "disease-control signal", "toxicity signal"],
            available_count=len(scenario_rows),
            date_range=_date_range_for_runs(scenario_rows),
            status="enough" if scenario_rows else "missing",
            used_by_model=True,
            why_it_matters="The scenarios compare alternative treatment schedules under the current mechanistic model.",
            missing_data_message="No completed simulated scenarios are currently available." if not scenario_rows else "Completed simulated scenarios are available for plain-language review.",
            next_action=_block_action(
                title="Review scenario differences",
                detail="Compare schedule meaning, exposure pattern, and plain-language interpretation before opening the technical cockpit.",
                href=reverse("twin_engine:simple_research_view", args=[patient.id]) + "#scenario-comparison",
            ),
            developer_link=_developer_link(include_developer_links, reverse("twin_engine:research_cockpit", args=[patient.id]) + "#what-if-scenarios"),
            limit="Scenario differences are model outputs, not identified causal effects.",
        ),
    ]
    return blocks


def _build_model_story(current_state, scenario_rows: list[dict[str, Any]], validation_panel: dict[str, Any]) -> dict[str, Any]:
    backtest_summary = (validation_panel.get("backtest") or {}).get("summary") or {}
    robustness_summary = (validation_panel.get("robustness") or {}).get("summary") or {}
    uncertainty_available = any((row.get("uncertainty") or {}).get("status") == "completed" for row in scenario_rows)
    toxicity_status = "prototype risk signals available" if any(row.get("toxicity_model_status") == "semi_mechanistic_prototype" for row in scenario_rows) else "observed context only"
    calibrated = bool(backtest_summary) or any(row.get("status") == "post_calibration" for row in [])
    return {
        "twin_initialized": "yes" if current_state else "no",
        "calibrated": "yes" if calibrated else "no",
        "simulations_available": "yes" if scenario_rows else "no",
        "uncertainty_available": "yes" if uncertainty_available else "no",
        "backtesting_available": "yes" if backtest_summary.get("status") == "completed" or backtest_summary else "no",
        "toxicity_model_status": toxicity_status,
        "causality_status": "Causal effect not identified. This page does not estimate the effect of changing treatment under causal assumptions.",
        "ranking_stability_available": "yes" if robustness_summary else "no",
    }


def _build_scenario_story(scenario_rows: list[dict[str, Any]], validation_panel: dict[str, Any]) -> list[dict[str, Any]]:
    robustness_summary = (validation_panel.get("robustness") or {}).get("summary") or {}
    robustness_rows = {row.get("scenario_label"): row for row in (robustness_summary.get("rows") or [])}
    stories = []
    for row in scenario_rows:
        robustness_row = robustness_rows.get(row["label"]) or {}
        utility_summary = ((row.get("uncertainty") or {}).get("metric_summaries") or {}).get("research_utility_v2") or {}
        stories.append(
            {
                "scenario_label": row["label"],
                "plain_schedule": _plain_schedule(row),
                "exposure_pattern": _plain_exposure_pattern(row),
                "average_daily_exposure": _format_number(row.get("average_daily_exposure_mg"), suffix=" mg/day"),
                "peak_daily_dose": _format_number(row.get("peak_daily_dose_mg"), suffix=" mg"),
                "tumor_signal_change": _signal_change_label(row.get("tumor_reduction_delta"), positive_label="stronger disease-control signal", negative_label="weaker disease-control signal"),
                "healthy_signal_change": _signal_change_label(_invert_sign(row.get("healthy_loss_delta")), positive_label="less healthy-tissue strain", negative_label="more healthy-tissue strain"),
                "toxicity_signal": _toxicity_label(row),
                "uncertainty_summary": _uncertainty_summary(utility_summary),
                "robustness_label": _ranking_stability_label(robustness_row),
                "plain_interpretation": _plain_interpretation(row, robustness_row, utility_summary),
                "allowed_conclusion": "This simulated schedule changes the model output under the current mechanistic assumptions.",
                "forbidden_conclusion": "Do not read this as proof of the best treatment, clinical superiority, or a causal effect.",
                "technical_metrics": {
                    "tumor_reduction": row.get("tumor_reduction_delta"),
                    "healthy_loss": row.get("healthy_loss_delta"),
                    "durability": row.get("durability_delta"),
                    "utility_v2": row.get("utility_v2"),
                    "p05": utility_summary.get("p05"),
                    "median": utility_summary.get("median"),
                    "p95": utility_summary.get("p95"),
                },
            }
        )
    return stories


def _build_conclusion_cards(
    data_blocks: list[dict[str, Any]],
    scenario_rows: list[dict[str, Any]],
    validation_panel: dict[str, Any],
    next_action: dict[str, str],
) -> list[dict[str, str]]:
    uncertainty_rows = validation_panel.get("uncertainty_rows") or []
    wide_rows = [row["label"] for row in uncertainty_rows if ((row.get("research_utility_v2") or {}).get("uncertainty_classification") == "wide")]
    stable_rows = [row["label"] for row in uncertainty_rows if ((row.get("research_utility_v2") or {}).get("uncertainty_classification") == "narrow")]
    missing_titles = [block["title"] for block in data_blocks if block["status"] == "missing"]
    robustness_summary = (validation_panel.get("robustness") or {}).get("summary") or {}
    contested_rows = [row.get("scenario_label") for row in robustness_summary.get("rows") or [] if row.get("robustness_classification") == "contested"]
    return [
        {
            "title": "What looks stable",
            "detail": ", ".join(stable_rows[:3]) + " stay within a narrow uncertainty range under the current stored checks." if stable_rows else "No scenario is yet marked as stable under the current stored checks.",
        },
        {
            "title": "What looks uncertain",
            "detail": ", ".join(wide_rows[:3]) + " show a wide uncertainty range." if wide_rows else ("Top ranks are contested across stored uncertainty samples." if contested_rows else "Uncertainty remains a research diagnostic and should be read cautiously."),
        },
        {
            "title": "What is missing",
            "detail": ", ".join(missing_titles[:3]) + " still need more structured records." if missing_titles else "The main missing pieces are now about better completeness, not a missing page section.",
        },
        {
            "title": "What is not yet scientifically identified",
            "detail": "Causal effect not identified. The platform compares simulated schedules, not proven treatment effects.",
        },
        {
            "title": "Next model improvement",
            "detail": next_action["detail"],
        },
    ]


def _build_one_minute_summary(
    patient: Patient,
    current_state,
    scenario_rows: list[dict[str, Any]],
    validation_panel: dict[str, Any],
    data_blocks: list[dict[str, Any]],
) -> list[dict[str, str]]:
    available_blocks = [block for block in data_blocks if block["status"] == "enough"]
    partial_blocks = [block for block in data_blocks if block["status"] == "partial"]
    backtest_summary = (validation_panel.get("backtest") or {}).get("summary") or {}
    uncertainty_available = any((row.get("uncertainty") or {}).get("status") == "completed" for row in scenario_rows)
    start_end = _patient_date_range(patient)
    return [
        {
            "title": "Data status",
            "status": "enough" if available_blocks else "partial",
            "detail": f"Structured longitudinal data are available from {start_end}, but some biomarkers are incomplete." if partial_blocks else f"Structured longitudinal data are available from {start_end}.",
        },
        {
            "title": "Model status",
            "status": "enough" if current_state else "missing",
            "detail": "A mathematical starting state exists and simulations are available." if current_state and scenario_rows else ("A mathematical starting state exists, but scenario runs are still missing." if current_state else "No current mathematical starting state is available yet."),
        },
        {
            "title": "Scenario status",
            "status": "enough" if scenario_rows else "missing",
            "detail": f"{len(scenario_rows)} scenarios have been simulated." if scenario_rows else "No simulated scenarios are available yet.",
        },
        {
            "title": "Main limitation",
            "status": "partial" if uncertainty_available or backtest_summary else "missing",
            "detail": "Uncertainty and toxicity outputs remain research diagnostics." if uncertainty_available or backtest_summary else "Observed data can be reviewed, but deeper diagnostics are still missing.",
        },
    ]


def _build_steps(
    patient: Patient,
    current_state,
    scenario_rows: list[dict[str, Any]],
    validation_panel: dict[str, Any],
) -> list[dict[str, str]]:
    backtest_summary = (validation_panel.get("backtest") or {}).get("summary") or {}
    uncertainty_available = any((row.get("uncertainty") or {}).get("status") == "completed" for row in scenario_rows)
    return [
        {
            "step": "1. Data collected",
            "status": "enough" if patient.assessments.exists() else "missing",
            "plain_explanation": "Recorded assessments, treatments, and lab rows are present for this pseudonymized patient." if patient.assessments.exists() else "No assessment history is available yet.",
            "link": "#data-blocks",
        },
        {
            "step": "2. Data structured",
            "status": "enough" if patient.longitudinal_lab_results.exists() else "partial",
            "plain_explanation": "Structured rows are available for several dated fields that the model can read." if patient.longitudinal_lab_results.exists() else "The page relies mostly on assessment snapshots because structured lab rows are limited.",
            "link": "#data-lineage",
        },
        {
            "step": "3. Mathematical starting state created",
            "status": "enough" if current_state else "missing",
            "plain_explanation": "A twin state means a mathematical starting state for simulation, not a direct copy of the patient." if current_state else "No mathematical starting state has been created yet.",
            "link": "#model-status",
        },
        {
            "step": "4. Model calibrated",
            "status": "enough" if backtest_summary else "partial",
            "plain_explanation": "Observed-versus-predicted checks are available." if backtest_summary else "A technical fit check is still limited or missing.",
            "link": "#one-minute-summary",
        },
        {
            "step": "5. Scenarios simulated",
            "status": "enough" if scenario_rows else "missing",
            "plain_explanation": "Alternative schedules have been simulated under the current mechanistic model." if scenario_rows else "No scenario runs are available yet.",
            "link": "#scenario-comparison",
        },
        {
            "step": "6. Outputs compared",
            "status": "enough" if scenario_rows else "missing",
            "plain_explanation": "The page compares schedule meaning, exposure pattern, and plain-language interpretation." if scenario_rows else "There are no scenario outputs to compare yet.",
            "link": "#what-the-model-uses",
        },
        {
            "step": "7. Limitations checked",
            "status": "enough" if uncertainty_available or backtest_summary else "partial",
            "plain_explanation": "The page calls out uncertainty ranges, missing data, and non-clinical limits." if uncertainty_available or backtest_summary else "The platform still needs more limitation diagnostics before interpretation becomes stronger.",
            "link": "#conclusion-limits",
        },
    ]


def _build_input_table(data_blocks: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows = []
    for block in data_blocks:
        rows.append(
            {
                "input_group": block["title"],
                "used_fields": ", ".join(block["fields_used"]),
                "model_role": block["why_it_matters"],
                "available": "Yes" if block["used_by_model"] and block["available_count"] else "No",
                "missing": block["missing_data_message"],
                "interpretation_risk": block["limit"],
            }
        )
    return rows


def _build_lineage_nodes(
    patient: Patient,
    current_state,
    scenario_rows: list[dict[str, Any]],
    validation_panel: dict[str, Any],
) -> list[dict[str, str]]:
    uncertainty_available = any((row.get("uncertainty") or {}).get("status") == "completed" for row in scenario_rows)
    robustness_available = bool(((validation_panel.get("robustness") or {}).get("summary") or {}).get("rows"))
    return [
        {"label": "Clinical documents", "status": "partial", "detail": "Source records exist, but this page shows only structured research-safe summaries.", "href": reverse("clinic:patient_detail", args=[patient.id])},
        {"label": "Structured patient records", "status": "enough" if patient.assessments.exists() else "missing", "detail": "Assessments, therapies, and structured lab rows are the visible input layer.", "href": "#data-blocks"},
        {"label": "Model inputs", "status": "enough" if patient.assessments.exists() else "missing", "detail": "Only selected structured fields feed the current model.", "href": "#what-the-model-uses"},
        {"label": "Twin state", "status": "enough" if current_state else "missing", "detail": "Twin state means mathematical starting state.", "href": "#model-status"},
        {"label": "Simulated scenarios", "status": "enough" if scenario_rows else "missing", "detail": "Alternative schedules are simulated from the current mathematical starting state.", "href": "#scenario-comparison"},
        {"label": "Uncertainty / ranking stability checks", "status": "enough" if uncertainty_available or robustness_available else "partial", "detail": "These checks show whether rankings stay stable under stored perturbation tests.", "href": "#scenario-comparison"},
        {"label": "Interpretation cards", "status": "enough", "detail": "The final cards summarize what is stable, uncertain, missing, and out of scope.", "href": "#conclusion-limits"},
    ]


def _make_block(**kwargs: Any) -> dict[str, Any]:
    return kwargs


def _record_count_summary(patient: Patient, scenario_rows: list[dict[str, Any]]) -> str:
    return (
        f"{patient.assessments.count()} assessments, "
        f"{patient.longitudinal_lab_results.count()} structured lab rows, "
        f"{patient.therapies.count()} treatments, "
        f"{patient.twin_states.count()} mathematical starting states, "
        f"{len(scenario_rows)} simulated scenarios"
    )


def _model_readiness_status(
    patient: Patient,
    current_state,
    scenario_rows: list[dict[str, Any]],
    validation_panel: dict[str, Any],
) -> str:
    if not patient.assessments.exists():
        return "missing"
    if current_state is None or not scenario_rows:
        return "partial"
    backtest_summary = (validation_panel.get("backtest") or {}).get("summary") or {}
    return "enough" if backtest_summary else "partial"


def _determine_next_action(
    patient: Patient,
    current_state,
    scenario_rows: list[dict[str, Any]],
    validation_panel: dict[str, Any],
    recommendation: dict[str, Any],
) -> dict[str, str]:
    backtest_summary = (validation_panel.get("backtest") or {}).get("summary") or {}
    robustness_summary = (validation_panel.get("robustness") or {}).get("summary") or {}
    uncertainty_rows = validation_panel.get("uncertainty_rows") or []
    has_wide_uncertainty = any(((row.get("research_utility_v2") or {}).get("uncertainty_classification") == "wide") for row in uncertainty_rows)
    if not patient.assessments.exists():
        return _block_action(
            title="Add the first structured assessment",
            detail="No structured assessment is available yet, so the page cannot explain model inputs or run a mathematical starting state.",
            href=reverse("clinic:assessment_new", args=[patient.id]),
        )
    if not recommendation.get("has_usable_assessment"):
        return _block_action(
            title="Fill in the minimum model-input fields",
            detail="At least disease-marker and hemoglobin fields are needed before the mathematical starting state becomes reproducible.",
            href=reverse("clinic:assessment_new", args=[patient.id]),
        )
    if current_state is None:
        return _block_action(
            title="Create the mathematical starting state",
            detail="Initialize the model from the strongest available assessment before interpreting any scenario output.",
            href=reverse("twin_engine:research_cockpit", args=[patient.id]) + "#twin-state",
        )
    if not backtest_summary:
        return _block_action(
            title="Run or review the fit check",
            detail="Observed-versus-predicted checks help explain where the current model is stronger or weaker.",
            href=reverse("twin_engine:research_cockpit", args=[patient.id]) + "#validation-and-uncertainty",
        )
    if not scenario_rows:
        return _block_action(
            title="Run scenario comparisons",
            detail="No simulated scenario is available yet, so the page cannot compare alternative schedules.",
            href=reverse("twin_engine:research_cockpit", args=[patient.id]) + "#what-if-scenarios",
        )
    if has_wide_uncertainty or robustness_summary:
        return _block_action(
            title="Check ranking stability before choosing a story",
            detail="Inspect which scenario differences remain stable and which remain fragile under the stored uncertainty checks.",
            href=reverse("twin_engine:research_cockpit", args=[patient.id]) + "#validation-and-uncertainty",
        )
    return _block_action(
        title="Open the scientific cockpit for technical detail",
        detail="The simple page is ready; move to the technical cockpit only if you need formulas, fit diagnostics, or developer checks.",
        href=reverse("twin_engine:research_cockpit", args=[patient.id]),
    )


def _block_action(*, title: str, detail: str, href: str) -> dict[str, str]:
    return {"title": title, "detail": detail, "href": href}


def _developer_link(include_developer_links: bool, href: str) -> str | None:
    return href if include_developer_links else None


def _patient_date_range(patient: Patient) -> str:
    ranges = [
        _date_range_from_dates([assessment.date for assessment in patient.assessments.all()]),
        _date_range_from_dates([lab.date for lab in patient.longitudinal_lab_results.all()]),
        _date_range_from_dates([therapy.start_date for therapy in patient.therapies.all() if therapy.start_date]),
    ]
    return _merge_ranges(*ranges)


def _field_counts(
    patient: Patient,
    field_specs: list[tuple[str, str | None, str | None]],
) -> dict[str, int]:
    assessments = list(patient.assessments.all())
    counts: dict[str, int] = {}
    for label, assessment_field, analyte in field_specs:
        assessment_count = 0
        if assessment_field:
            assessment_count = sum(1 for assessment in assessments if getattr(assessment, assessment_field, None) is not None)
        lab_count = 0
        if analyte:
            lab_count = patient.longitudinal_lab_results.filter(analyte=analyte, value__isnull=False).count()
        counts[label] = assessment_count + lab_count
    return counts


def _status_from_field_counts(field_counts: dict[str, int], *, enough_threshold: int) -> str:
    available_fields = sum(1 for count in field_counts.values() if count > 0)
    if available_fields >= enough_threshold:
        return "enough"
    if available_fields > 0:
        return "partial"
    return "missing"


def _missing_message(field_counts: dict[str, int], fallback: str) -> str:
    missing_fields = [label for label, count in field_counts.items() if count == 0]
    if not missing_fields:
        return "No major gap is flagged in this block right now."
    return fallback + " Missing now: " + ", ".join(missing_fields) + "."


def _date_range_for_labs(
    patient: Patient,
    analytes: list[str],
    *,
    assessment_fields: list[str] | None = None,
) -> str:
    lab_dates = list(patient.longitudinal_lab_results.filter(analyte__in=analytes).values_list("date", flat=True))
    assessment_dates: list[date] = []
    for field in assessment_fields or []:
        assessment_dates.extend(
            [assessment.date for assessment in patient.assessments.all() if getattr(assessment, field, None) is not None]
        )
    return _date_range_from_dates([*lab_dates, *assessment_dates])


def _date_range_for_events(events: Iterable[Any], allowed_types: set[str] | None = None) -> str:
    filtered = [event.date for event in events if allowed_types is None or getattr(event, "event_type", None) in allowed_types]
    return _date_range_from_dates(filtered)


def _date_range_for_therapies(therapies: list[PatientTherapy]) -> str:
    dates = [therapy.start_date for therapy in therapies if therapy.start_date]
    dates.extend([therapy.end_date for therapy in therapies if therapy.end_date])
    return _date_range_from_dates(dates)


def _date_range_for_interruptions(interruptions: list[TherapyInterruption]) -> str:
    dates = [row.start_date for row in interruptions if row.start_date]
    dates.extend([row.end_date for row in interruptions if row.end_date])
    return _date_range_from_dates(dates)


def _date_range_for_states(patient: Patient) -> str:
    return _date_range_from_dates(list(patient.twin_states.values_list("state_date", flat=True)))


def _date_range_for_runs(scenario_rows: list[dict[str, Any]]) -> str:
    return _date_range_from_dates([row["run"].created_at.date() for row in scenario_rows if row.get("run")])


def _date_range_from_dates(values: Iterable[date | None]) -> str:
    dates = sorted({value for value in values if value is not None})
    if not dates:
        return "No dated records"
    if len(dates) == 1:
        return dates[0].isoformat()
    return f"{dates[0].isoformat()} to {dates[-1].isoformat()}"


def _merge_ranges(*ranges: str) -> str:
    dates: list[date] = []
    for text in ranges:
        if not text or text == "No dated records":
            continue
        parts = [part.strip() for part in text.split("to")]
        for part in parts:
            try:
                dates.append(date.fromisoformat(part))
            except ValueError:
                continue
    return _date_range_from_dates(dates)


def _plain_schedule(row: dict[str, Any]) -> str:
    drug = row.get("primary_drug") or "Primary drug"
    schedule_type = (row.get("schedule_type") or "recorded schedule").replace("_", " ")
    peak = row.get("peak_daily_dose_mg")
    average = row.get("average_daily_exposure_mg")
    if peak and average and peak != average:
        return f"{drug} with a {schedule_type} pattern, peaking near {peak:.1f} mg and averaging {average:.1f} mg/day."
    if peak:
        return f"{drug} with a {schedule_type} pattern around {peak:.1f} mg."
    return f"{drug} with a {schedule_type} pattern."


def _plain_exposure_pattern(row: dict[str, Any]) -> str:
    temporal_profile = row.get("temporal_profile") or "unavailable"
    average = row.get("average_daily_exposure_mg")
    peak = row.get("peak_daily_dose_mg")
    if average is None and peak is None:
        return f"Exposure pattern: {temporal_profile}."
    return f"Exposure pattern: {temporal_profile}; average {average or 0:.1f} mg/day and peak {peak or 0:.1f} mg."


def _signal_change_label(value: float | None, *, positive_label: str, negative_label: str) -> str:
    if value is None:
        return "Not available"
    if value > 0.05:
        return positive_label
    if value < -0.05:
        return negative_label
    return "little visible change"


def _invert_sign(value: float | None) -> float | None:
    if value is None:
        return None
    return value * -1.0


def _toxicity_label(row: dict[str, Any]) -> str:
    liver = row.get("liver_signal")
    neutropenia = row.get("neutropenia_signal")
    highest = max([value for value in [liver, neutropenia] if value is not None], default=None)
    if highest is None:
        return "No stored toxicity summary"
    if highest >= 0.6:
        return "higher research toxicity signal"
    if highest >= 0.3:
        return "moderate research toxicity signal"
    return "lower research toxicity signal"


def _uncertainty_summary(summary: dict[str, Any]) -> str:
    if not summary:
        return "No stored uncertainty range"
    classification = summary.get("uncertainty_classification") or "unavailable"
    median = summary.get("median")
    if median is None:
        return f"Stored uncertainty range is {classification}."
    return f"Stored uncertainty range is {classification} around a median of {median:.3f}."


def _ranking_stability_label(robustness_row: dict[str, Any]) -> str:
    classification = robustness_row.get("robustness_classification")
    if not classification:
        return "No ranking-stability check stored"
    return classification.replace("_", " ")


def _plain_interpretation(
    row: dict[str, Any],
    robustness_row: dict[str, Any],
    utility_summary: dict[str, Any],
) -> str:
    disease = _signal_change_label(row.get("tumor_reduction_delta"), positive_label="stronger disease-control signal", negative_label="weaker disease-control signal")
    healthy = _signal_change_label(_invert_sign(row.get("healthy_loss_delta")), positive_label="less healthy-tissue strain", negative_label="more healthy-tissue strain")
    stability = _ranking_stability_label(robustness_row)
    uncertainty = utility_summary.get("uncertainty_classification")
    parts = [f"This scenario shows {disease}"]
    if healthy != "Not available":
        parts.append(f"with {healthy}")
    if uncertainty:
        parts.append(f"and a {uncertainty} uncertainty range")
    if stability and stability != "No ranking-stability check stored":
        parts.append(f"while ranking stability is {stability}")
    return " ".join(parts) + "."


def _format_number(value: float | None, *, suffix: str = "") -> str:
    if value is None:
        return "Not available"
    return f"{value:.2f}{suffix}"

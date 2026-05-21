from __future__ import annotations

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

from django.apps import apps

Regimen = apps.get_model("clinic", "Regimen")

from . import explain, forms, models
from .presets import PRESETS
from .models_help import HelpArticle
from .pharmaco import registry as pharmaco_registry
from .permissions import is_editor


def scenario_list(request):
    scenarios = models.Scenario.objects.filter(active=True).prefetch_related("recommended_regimens")
    twin_assessment_id = (request.GET.get("twin_assessment_id") or "").strip()
    twin_label = ""
    twin_error = ""
    if twin_assessment_id and twin_assessment_id.isdigit() and getattr(request.user, "is_authenticated", False):
        try:
            from clinic.models import Assessment
            from .access import accessible_assessments

            assessment_pk = int(twin_assessment_id)
            a = accessible_assessments(request.user, base_qs=Assessment.objects.select_related("patient")).filter(pk=assessment_pk).first()
            if a and getattr(a, "patient", None):
                first = (a.patient.first_name or "").strip()
                last = (a.patient.last_name or "").strip()
                full_name = (f"{first} {last}").strip() or last or first
                twin_label = f"{a.patient.mrn} · {full_name} · {a.date}"
            else:
                twin_error = "Selected assessment is not accessible (permissions)."
        except Exception:
            twin_error = "Could not resolve selected assessment."

    quickstart_scenario = scenarios.order_by("pk").first()

    context = {
        "scenarios": scenarios,
        "twin_assessment_id": twin_assessment_id,
        "twin_label": twin_label,
        "twin_error": twin_error,
        "quickstart_scenario": quickstart_scenario,
    }
    return render(request, "simulator/scenario_list.html", context)


@login_required
def getting_started(request):
    """Getting started page with tutorials and practice scenarios."""
    return render(request, "simulator/getting_started.html")


@login_required
def interactive_tutorial(request):
    """Interactive step-by-step tutorial for creating a patient and running first simulation."""
    from clinic.models import Patient, Assessment
    
    # Check if demo patients exist
    demo_patients = Patient.objects.filter(mrn__startswith='DEMO').prefetch_related('assessments')
    has_demo_data = demo_patients.exists()
    
    context = {
        'has_demo_data': has_demo_data,
        'demo_patients': demo_patients,
        'demo_patients_count': demo_patients.count(),
    }
    return render(request, "simulator/interactive_tutorial.html", context)


@login_required
def visibility_diagnostics(request):
    """Explain what assessments a user can/can't use (and why)."""
    from clinic.models import Assessment
    from .access import DEMO_MRN_PREFIX, accessible_assessments

    all_qs = Assessment.objects.select_related("patient")
    total_assessments = all_qs.count()

    is_privileged = request.user.is_staff or is_editor(request.user)
    if is_privileged:
        accessible_count = total_assessments
        inaccessible_count = 0
        inaccessible_sample = []
    else:
        access_q = Q(patient__owner=request.user) | Q(patient__mrn__startswith=DEMO_MRN_PREFIX)
        accessible_count = accessible_assessments(request.user, base_qs=all_qs).count()
        inaccessible_count = all_qs.exclude(access_q).count()
        inaccessible_sample = list(all_qs.exclude(access_q).order_by("-date")[:50])

    context = {
        "is_privileged": is_privileged,
        "total_assessments": total_assessments,
        "accessible_count": accessible_count,
        "inaccessible_count": inaccessible_count,
        "inaccessible_sample": inaccessible_sample,
        "demo_prefix": DEMO_MRN_PREFIX,
    }
    return render(request, "simulator/visibility_diagnostics.html", context)


@login_required
def scenario_detail(request, pk: int):
    scenario = get_object_or_404(
        models.Scenario.objects.prefetch_related("recommended_regimens", "attempts__selected_regimen", "attempts__user"),
        pk=pk,
        active=True,
    )
    # Hide legacy/empty entries (e.g. saved without regimen/response/notes and without simulation outputs)
    # to reduce beginner confusion.
    attempts_all = list(scenario.attempts.all())
    attempts = [
        a
        for a in attempts_all
        if a.selected_regimen_id
        or (a.predicted_response or "").strip()
        or (a.notes or "").strip()
        or bool(a.results_summary)
    ]
    hidden_attempts_count = max(0, len(attempts_all) - len(attempts))
    form = forms.SimulationAttemptForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        attempt = form.save(commit=False)
        attempt.scenario = scenario
        attempt.user = request.user
        regimen = attempt.selected_regimen
        if regimen and scenario.recommended_regimens.filter(pk=regimen.pk).exists():
            attempt.is_guideline_aligned = True
        attempt.save()
        if attempt.is_guideline_aligned:
            messages.success(
                request,
                "Plan recorded. Great choice—aligned with the guideline set for this scenario.",
            )
        else:
            messages.warning(
                request,
                "Plan recorded. Review the guideline notes below to compare approaches.",
            )
        return redirect(reverse("simulator:scenario_detail", args=[scenario.pk]))

    recommended_regimens = scenario.recommended_regimens.order_by("name")
    regimen_names = ", ".join(regimen.name for regimen in recommended_regimens) if recommended_regimens else "No guideline regimen linked yet."
    editor = is_editor(request.user)
    available_regimens = (
        Regimen.objects.exclude(pk__in=recommended_regimens.values_list("pk", flat=True)).order_by("name")
        if editor
        else Regimen.objects.none()
    )
    latest_simulation = (
        scenario.attempts.exclude(results_summary={}).order_by("-submitted").first()
    )
    latest_summary = latest_simulation.results_summary if latest_simulation else None
    latest_results = latest_simulation.results if latest_simulation else {}
    game_mode = request.GET.get("game") == "1"
    from .game import compute_game_metrics
    game = compute_game_metrics(latest_summary) if game_mode else None
    latest_warnings: list[str] = []
    if latest_summary:
        healthy_loss = latest_summary.get("healthy_loss")
        if healthy_loss is not None:
            if healthy_loss > 0.3:
                latest_warnings.append("Healthy cell loss exceeded 30% in the most recent simulation—dose reduction recommended.")
            elif healthy_loss > 0.2:
                latest_warnings.append("Healthy cell loss above 20% warrants close monitoring for toxicity.")
        tumor_reduction = latest_summary.get("tumor_reduction")
        if tumor_reduction is not None and tumor_reduction < 0:
            latest_warnings.append("Latest simulation predicted tumor regrowth (negative reduction). Adjust parameters and re-run.")
    def _profile_with_ranges(drug: str):
        profile = pharmaco_registry.get_drug_profile(drug)
        if not profile:
            return None
        profile_copy = dict(profile)
        dose = profile_copy.get("dose_range", {})
        span = f"{dose.get('min')}–{dose.get('max')} {profile_copy.get('unit', '')}".strip()
        profile_copy["range_en"] = f"Allowed: {span}"
        profile_copy["range_it"] = f"Consentito: {span}"
        return profile_copy

    drug_profiles = {
        "lenalidomide": _profile_with_ranges("lenalidomide"),
        "bortezomib": _profile_with_ranges("bortezomib"),
        "daratumumab": _profile_with_ranges("daratumumab"),
        "carfilzomib": _profile_with_ranges("carfilzomib"),
    }
    guide_slugs = ["quickstart", "optimization_lab"]
    guides: dict[str, dict[str, dict[str, str]]] = {}
    articles = {article.slug: article for article in HelpArticle.objects.filter(slug__in=guide_slugs)}
    for slug in guide_slugs:
        article = articles.get(slug)
        if not article:
            continue
        guides[slug] = {
            "en": article.as_lang("en"),
            "it": article.as_lang("it"),
        }
    preset_descriptions = {
        key: {
            "en": preset.get("description_en", ""),
            "it": preset.get("description_it", ""),
            "story_en": preset.get("story_en", {}),
            "story_it": preset.get("story_it", {}),
        }
        for key, preset in PRESETS.items()
    }
    help_index = [
        {"slug": article.slug, "title_en": article.title_en, "title_it": article.title_it}
        for article in HelpArticle.objects.order_by("slug")
    ]

    twin_assessment_id = (request.GET.get("twin_assessment_id") or "").strip()
    twin_label = ""
    twin_patient_pk = None
    twin_error = ""
    sim_initial = {}
    if twin_assessment_id.isdigit():
        sim_initial["twin_assessment_id"] = int(twin_assessment_id)
    prefill_twin = bool(sim_initial)

    if prefill_twin:
        try:
            from clinic.models import Assessment
            from .access import accessible_assessments

            assessment_pk = int(twin_assessment_id)
            a = (
                accessible_assessments(request.user, base_qs=Assessment.objects.select_related("patient"))
                .filter(pk=assessment_pk)
                .first()
            )
            if a and getattr(a, "patient", None):
                first = (a.patient.first_name or "").strip()
                last = (a.patient.last_name or "").strip()
                full_name = (f"{first} {last}").strip() or last or first
                twin_label = f"{a.patient.mrn} · {full_name} · {a.date}"
                twin_patient_pk = a.patient.pk
            else:
                twin_error = "Selected assessment is not accessible (permissions)."
        except Exception:
            twin_error = "Could not resolve selected assessment."

    simulation_runs_qs = scenario.attempts.exclude(results_summary={}).select_related("user").order_by("-submitted")
    decision_logs_qs = (
        scenario.attempts.filter(
            Q(selected_regimen__isnull=False)
            | ~Q(predicted_response="")
            | ~Q(notes="")
        )
        .select_related("user", "selected_regimen")
        .order_by("-submitted")
    )
    simulation_runs_count = simulation_runs_qs.count()
    decision_logs_count = decision_logs_qs.count()
    simulation_runs = list(simulation_runs_qs[:20])
    decision_logs = list(decision_logs_qs[:20])

    context = {
        "scenario": scenario,
        "attempts": attempts,
        "hidden_attempts_count": hidden_attempts_count,
        "form": form,
        "recommended_regimens": recommended_regimens,
        "regimen_names": regimen_names,
        "is_editor": editor,
        "available_regimens": available_regimens,
        "simulation_parameter_form": forms.SimulationParameterForm(user=request.user, initial=sim_initial),
        "prefill_twin": prefill_twin,
        "latest_simulation": latest_simulation,
        "latest_simulation_summary": latest_summary,
        "latest_simulation_results": latest_results,
        "latest_simulation_warnings": latest_warnings,
        "game_mode": game_mode,
        "game": game,
        "drug_profiles": drug_profiles,
        "sim_form_help_it": forms.SIMULATION_FORM_HELP_TEXT_IT,
        "sim_form_help_en": forms.SIMULATION_FORM_HELP_TEXT_EN,
        "kpi": explain.KPI,
        "guide_articles": guides,
        "help_index": help_index,
        "preset_descriptions": preset_descriptions,
        "twin_assessment_id": twin_assessment_id,
        "twin_label": twin_label,
        "twin_patient_pk": twin_patient_pk,
        "twin_error": twin_error,
        "simulation_runs": simulation_runs,
        "simulation_runs_count": simulation_runs_count,
        "decision_logs": decision_logs,
        "decision_logs_count": decision_logs_count,
    }
    return render(request, "simulator/scenario_detail.html", context)


@login_required
def algorithm_transparency(request):
    """Display exploratory simulation rule documentation for audit and review."""
    from .decision_algorithm import get_algorithm
    
    algorithm = get_algorithm()
    patient_context = _build_algorithm_patient_context(request)
    
    context = {
        "version": algorithm.get("version"),
        "last_updated": algorithm.get("last_updated"),
        "page_purpose": "Use this page to inspect how the simulator turns encoded thresholds and rule logic into exploratory workflow labels.",
        "interpretation_status": {
            "title": "Algorithm logic transparency",
            "summary": "This page documents heuristic and literature-informed rules used for exploratory simulation support.",
            "limit": "This page is not a prescribing tool and does not establish patient-specific comparative benefit.",
        },
        "rule_family_rows": _build_algorithm_rule_family_rows(),
        "threshold_sections": _build_algorithm_threshold_sections(algorithm.get("thresholds", {})),
        "decision_rule_rows": _build_algorithm_decision_rule_rows(algorithm.get("decision_rules", [])),
        "risk_context_rows": _build_algorithm_risk_context_rows(algorithm.get("risk_stratification", {})),
        "cytogenetic_flag_rows": _build_algorithm_cytogenetic_flag_rows(algorithm.get("high_risk_cytogenetics", {})),
        "data_sources": algorithm.get("data_sources", []),
        "what_can_be_concluded": [
            "The page explains how rule-based labels are produced.",
            "The page helps audit which inputs affect simulation-support labels.",
            "The page exposes threshold assumptions and limitations.",
        ],
        "what_cannot_be_concluded": [
            "It cannot establish individual regimen advantage.",
            "It cannot infer patient-specific clinical benefit.",
            "It cannot prove outcome under an alternative intervention.",
            "It cannot replace counterfactual simulation, calibration, or clinical validation.",
        ],
        "next_actions": _build_algorithm_next_actions(patient_context),
    }
    return render(request, "simulator/algorithm_transparency.html", context)


def _build_algorithm_rule_family_rows():
    return [
        {
            "rule_family": "Metric threshold bands",
            "inputs_used": "tumor_reduction, healthy_loss, time_to_recurrence",
            "classification": "THRESHOLD-BASED",
            "purpose": "Turns continuous simulation outputs into stable rule-based context labels used by the interface.",
            "limitation": "Threshold bands do not prove benefit, superiority, or individual outcome.",
        },
        {
            "rule_family": "Decision-rule outputs",
            "inputs_used": "threshold labels, time_horizon, combined output states",
            "classification": "SIMULATION SUPPORT",
            "purpose": "Surfaces exploratory rule outputs, constraint flags, and workflow cues when encoded triggers fire.",
            "limitation": "These outputs document simulator behavior; they do not prescribe care.",
        },
        {
            "rule_family": "Risk/context staging",
            "inputs_used": "beta-2 microglobulin, albumin, LDH, cytogenetic context",
            "classification": "LITERATURE-INFORMED",
            "purpose": "Documents how published staging references are translated into risk/context flags and reference anchors.",
            "limitation": "Population staging references do not establish an individual outcome.",
        },
        {
            "rule_family": "Narrative label mapping",
            "inputs_used": "priority levels, threshold outcomes, evidence lists",
            "classification": "HEURISTIC",
            "purpose": "Maps encoded results into readable simulation-support labels and audit text for users.",
            "limitation": "Readable labels are platform-defined explanations, not clinical instructions.",
        },
        {
            "rule_family": "Unmodeled or absent factors",
            "inputs_used": "comorbidities, patient preferences, local practice, undocumented covariates",
            "classification": "UNKNOWN",
            "purpose": "Highlights decision-relevant factors that are not encoded in this page's rule families.",
            "limitation": "Missing inputs cannot be inferred by the rule set or by this page.",
        },
    ]


def _build_algorithm_threshold_sections(thresholds):
    category_specs = {
        "efficacy": {
            "title": "Efficacy thresholds",
            "summary": "Thresholds that label tumor-reduction outputs for simulation support review.",
            "variable": "tumor_reduction",
            "basis_label": "HEURISTIC",
            "workflow_change": "Changes the efficacy simulation-support label used in summaries and rule evaluation.",
            "not_prove": "Does not prove clinical response quality or regimen advantage.",
        },
        "toxicity": {
            "title": "Toxicity thresholds",
            "summary": "Thresholds that label simulated healthy-cell loss for constraint review.",
            "variable": "healthy_loss",
            "basis_label": "HEURISTIC",
            "workflow_change": "Changes the toxicity or constraint flag and can activate toxicity-focused rule outputs.",
            "not_prove": "Does not prove tolerability for an individual patient.",
        },
        "durability": {
            "title": "Durability thresholds",
            "summary": "Thresholds that label recurrence timing for durability context inside the workflow.",
            "variable": "time_to_recurrence",
            "basis_label": "LITERATURE-INFORMED",
            "workflow_change": "Changes the durability or risk/context flag and can activate longer-horizon workflow cues.",
            "not_prove": "Does not prove real-world recurrence timing for an individual patient.",
        },
    }

    sections = []
    for category, spec in category_specs.items():
        category_thresholds = thresholds.get(category, {})
        items = []
        for level, threshold in category_thresholds.items():
            items.append(
                {
                    "level": level.replace("_", " ").title(),
                    "description": threshold.get("description_en", ""),
                    "condition": threshold.get("condition", ""),
                    "variable": spec["variable"],
                    "classification": "THRESHOLD-BASED",
                    "basis_label": spec["basis_label"],
                    "workflow_change": spec["workflow_change"],
                    "not_prove": spec["not_prove"],
                }
            )
        sections.append(
            {
                "title": spec["title"],
                "summary": spec["summary"],
                "items": items,
            }
        )
    return sections


def _build_algorithm_decision_rule_rows(decision_rules):
    rows = []
    for rule in decision_rules:
        rows.append(
            {
                "id": rule.get("id", "RULE"),
                "name": rule.get("name_en", "Unnamed rule"),
                "trigger_condition": rule.get("trigger_condition", ""),
                "priority": str(rule.get("priority", "medium")).upper(),
                "inputs_used": _build_algorithm_inputs_used(rule.get("trigger_condition", "")),
                "classification": "SIMULATION SUPPORT",
                "output_label": _soften_algorithm_action(rule.get("id", ""), rule.get("action_en", "")),
                "purpose": _build_algorithm_rule_purpose(rule.get("id", "")),
                "limitation": "Documents how the platform reacts to encoded triggers; it does not establish patient-specific clinical benefit.",
                "evidence": rule.get("evidence", []),
            }
        )
    return rows


def _build_algorithm_risk_context_rows(risk_stratification):
    rows = []
    for key, stage in risk_stratification.items():
        stage_label = key.replace("_", " ")
        median_os = stage.get("median_OS_months")
        os_label = f"{median_os} months" if median_os is not None else "Not reached"
        rows.append(
            {
                "name": stage_label,
                "criteria": stage.get("criteria", ""),
                "classification": "LITERATURE-INFORMED",
                "purpose": "Creates a risk/context flag and displays published survival anchors for audit.",
                "limitation": "Population staging anchors do not establish an individual regimen advantage or outcome.",
                "reference_values": f"Median PFS {stage.get('median_PFS_months')} months; Median OS {os_label}; 5-year survival {stage.get('five_year_survival_percent')}%",
            }
        )
    return rows


def _build_algorithm_cytogenetic_flag_rows(high_risk_cytogenetics):
    rows = []
    for _, cyto in high_risk_cytogenetics.items():
        rows.append(
            {
                "name": cyto.get("name", "Unknown cytogenetic flag"),
                "gene": cyto.get("gene", ""),
                "impact": cyto.get("impact_en", ""),
                "classification": "LITERATURE-INFORMED",
                "purpose": "Adds a risk/context flag documenting how high-risk cytogenetic findings are surfaced to the simulator workflow.",
                "limitation": "The flag does not prove that any individual regimen will outperform another.",
            }
        )
    return rows


def _build_algorithm_next_actions(patient_context):
    secondary = []
    if patient_context is not None:
        patient_id = patient_context["patient_id"]
        secondary = [
            {
                "label": "Open Simple Research View",
                "href": reverse("twin_engine:simple_research_view", args=[patient_id]),
                "disabled": False,
            },
            {
                "label": "Open Scientific Cockpit",
                "href": reverse("twin_engine:research_cockpit", args=[patient_id]),
                "disabled": False,
            },
            {
                "label": "Back to Patient",
                "href": patient_context["href"],
                "disabled": False,
            },
        ]
    else:
        secondary = [
            {
                "label": "Open Simple Research View",
                "href": "#",
                "disabled": True,
            },
            {
                "label": "Open Scientific Cockpit",
                "href": "#",
                "disabled": True,
            },
        ]

    return {
        "primary": {
            "label": "Back to simulator scenarios",
            "href": reverse("simulator:scenario_list"),
        },
        "secondary": secondary,
        "patient_context_note": (
            None
            if patient_context is not None
            else "Patient-specific research links appear when the page is opened with a patient context."
        ),
    }


def _build_algorithm_patient_context(request):
    patient_id = (request.GET.get("patient_id") or "").strip()
    if not patient_id.isdigit():
        return None

    from clinic import models as clinic_models
    from clinic.views import can_edit_patient

    patient = clinic_models.Patient.objects.filter(pk=int(patient_id)).first()
    if patient is None or not can_edit_patient(request.user, patient):
        return None

    return {
        "patient_id": patient.pk,
        "href": reverse("clinic:patient_detail", args=[patient.pk]),
    }


def _build_algorithm_inputs_used(trigger_condition):
    variables = []
    for variable in ["tumor_reduction", "healthy_loss", "time_to_recurrence", "time_horizon"]:
        if variable in trigger_condition:
            variables.append(variable)
    return ", ".join(variables) if variables else "Encoded workflow state"


def _build_algorithm_rule_purpose(rule_id):
    purposes = {
        "RULE_001": "Flags high-toxicity simulation outputs for constraint review.",
        "RULE_002": "Flags borderline-toxicity outputs for closer monitoring context.",
        "RULE_003": "Flags low-efficacy outputs for escalation-context review inside the simulator.",
        "RULE_004": "Flags simulated failure states and surfaces alternative exploratory comparison context.",
        "RULE_005": "Flags early-recurrence patterns for longer-horizon exploratory review.",
        "RULE_006": "Flags a favorable balance zone while keeping alternative exploratory review available.",
    }
    return purposes.get(rule_id, "Documents how an encoded rule changes simulation-support labels.")


def _soften_algorithm_action(rule_id, action_text):
    overrides = {
        "RULE_001": "Adds a constraint flag for lower-intensity or shorter-duration exploratory review.",
        "RULE_002": "Adds a monitoring/context flag for lower-intensity exploratory review when toxicity signals are present.",
        "RULE_003": "Adds an escalation-context flag for higher-intensity or longer-horizon exploratory review.",
        "RULE_004": "Adds a failure-context flag and surfaces alternative regimen families for exploratory comparison.",
        "RULE_005": "Adds a durability-context flag for longer-horizon exploratory review.",
        "RULE_006": "Adds a balance flag showing the current simulation sits in the platform's favorable zone.",
    }
    return overrides.get(rule_id, f"Exploratory rule output: {action_text}" if action_text else "Exploratory rule output.")

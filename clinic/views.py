from __future__ import annotations

import json
from collections import Counter
from datetime import date
import logging

from django.contrib.auth.decorators import login_required, user_passes_test
from django.db.models import Count, Prefetch, F
from django.http import HttpRequest, HttpResponse, HttpResponseForbidden, JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_http_methods

import django_filters

from . import forms, models
from .models_symptoms import SymptomAssessment
from .forms_symptoms import SymptomAssessmentForm, QuickSymptomForm

embed_debug_logger = logging.getLogger("embed_debug")


def _interpret_latest_simulation(summary: dict[str, object] | None, parameters: dict[str, object] | None) -> dict[str, object] | None:
    """Comprehensive decision support: interpretation + actionable recommendations."""
    if not summary:
        return None

    def f(name: str) -> float | None:
        v = summary.get(name)
        try:
            return float(v) if v is not None else None
        except Exception:
            return None

    tumor_reduction = f("tumor_reduction")
    healthy_loss = f("healthy_loss")
    time_to_recurrence = f("time_to_recurrence")

    def band_label(value: float | None, *, good_ge: float | None = None, bad_lt: float | None = None) -> str:
        if value is None:
            return "unknown"
        if good_ge is not None and value >= good_ge:
            return "good"
        if bad_lt is not None and value < bad_lt:
            return "bad"
        return "caution"

    # Heuristics (UX-only): meant for beginner orientation, not clinical advice.
    tumor_label = band_label(tumor_reduction, good_ge=0.50, bad_lt=0.0)
    # Lower healthy loss is better.
    if healthy_loss is None:
        healthy_label = "unknown"
    elif healthy_loss < 0.20:
        healthy_label = "good"
    elif healthy_loss < 0.30:
        healthy_label = "caution"
    else:
        healthy_label = "bad"

    # Longer time to recurrence is better (if present).
    if time_to_recurrence is None:
        recurrence_label = "unknown"
    elif time_to_recurrence >= 180:
        recurrence_label = "good"
    elif time_to_recurrence >= 90:
        recurrence_label = "caution"
    else:
        recurrence_label = "bad"

    # Overall: worst-of with a slight preference for toxicity (healthy_loss) warnings.
    labels = [tumor_label, healthy_label, recurrence_label]
    if "bad" in labels:
        overall = "bad"
    elif "caution" in labels:
        overall = "caution"
    elif "good" in labels:
        overall = "good"
    else:
        overall = "unknown"

    # Extract current parameters (for recommendations)
    param = parameters or {}
    time_horizon = param.get("time_horizon")
    try:
        time_horizon = int(time_horizon) if time_horizon is not None else 168
    except Exception:
        time_horizon = 168

    # Generate actionable recommendations based on results.
    recommendations: list[dict[str, str]] = []

    # Scenario 1: High toxicity (healthy_loss ≥ 0.30)
    if healthy_loss is not None and healthy_loss >= 0.30:
        recommendations.append({
            "issue_en": "High toxicity (healthy cell loss ≥30%)",
            "issue_it": "Alta tossicità (perdita di cellule sane ≥30%)",
            "action_en": "Reduce drug doses by 20–30% or shorten time horizon",
            "action_it": "Riduci le dosi dei farmaci del 20–30% o accorcia l'orizzonte temporale",
            "rationale_en": "Too much damage to healthy plasma cells. Lower doses preserve immune function.",
            "rationale_it": "Troppo danno alle plasmacellule sane. Dosi più basse preservano la funzione immunitaria.",
            "icon": "⚠️",
            "priority": "high",
        })
    elif healthy_loss is not None and healthy_loss >= 0.20:
        recommendations.append({
            "issue_en": "Moderate toxicity (healthy cell loss 20–30%)",
            "issue_it": "Tossicità moderata (perdita di cellule sane 20–30%)",
            "action_en": "Consider reducing doses by 10–15% if patient shows clinical toxicity signs",
            "action_it": "Considera di ridurre le dosi del 10–15% se il paziente mostra segni clinici di tossicità",
            "rationale_en": "Borderline toxicity. Monitor closely; reduce if side effects appear.",
            "rationale_it": "Tossicità al limite. Monitora attentamente; riduci se compaiono effetti collaterali.",
            "icon": "⚠️",
            "priority": "medium",
        })

    # Scenario 2: Poor efficacy (tumor_reduction < 0.30)
    if tumor_reduction is not None and tumor_reduction < 0.30:
        if healthy_loss is None or healthy_loss < 0.30:
            # Low efficacy but tolerable toxicity → try increasing doses or horizon
            recommendations.append({
                "issue_en": "Low tumor reduction (<30%)",
                "issue_it": "Bassa riduzione tumorale (<30%)",
                "action_en": "Increase drug doses by 15–25% or extend time horizon to 224–280 days",
                "action_it": "Aumenta le dosi dei farmaci del 15–25% o estendi l'orizzonte a 224–280 giorni",
                "rationale_en": "More aggressive therapy may improve response if toxicity remains acceptable.",
                "rationale_it": "Una terapia più aggressiva può migliorare la risposta se la tossicità resta accettabile.",
                "icon": "📈",
                "priority": "high",
            })

    # Scenario 3: Negative tumor reduction (tumor growth)
    if tumor_reduction is not None and tumor_reduction < 0:
        recommendations.append({
            "issue_en": "Tumor growth (negative reduction)",
            "issue_it": "Crescita tumorale (riduzione negativa)",
            "action_en": "Switch to a different regimen or significantly increase doses",
            "action_it": "Cambia regime terapeutico o aumenta significativamente le dosi",
            "rationale_en": "Current regimen is ineffective. Consider alternative drug combinations.",
            "rationale_it": "Il regime attuale è inefficace. Considera combinazioni alternative di farmaci.",
            "icon": "🚨",
            "priority": "critical",
        })

    # Scenario 4: Short time to recurrence (if measured)
    if time_to_recurrence is not None and time_to_recurrence < 90:
        if time_horizon < 200:
            recommendations.append({
                "issue_en": "Early recurrence predicted (<90 days)",
                "issue_it": "Recidiva precoce prevista (<90 giorni)",
                "action_en": "Extend time horizon to 224–280 days to simulate longer treatment",
                "action_it": "Estendi l'orizzonte a 224–280 giorni per simulare un trattamento più lungo",
                "rationale_en": "Longer therapy duration may delay recurrence and improve durability.",
                "rationale_it": "Una durata più lunga può ritardare la recidiva e migliorare la durabilità.",
                "icon": "⏱️",
                "priority": "medium",
            })

    # Scenario 5: Good balance
    if overall == "good":
        recommendations.append({
            "issue_en": "Favorable balance (good efficacy, acceptable toxicity)",
            "issue_it": "Equilibrio favorevole (buona efficacia, tossicità accettabile)",
            "action_en": "Fine-tune by testing ±10% dose variations or compare alternative regimens",
            "action_it": "Ottimizza testando variazioni di dose ±10% o confronta regimi alternativi",
            "rationale_en": "Current settings look promising. Minor adjustments may further optimize.",
            "rationale_it": "Le impostazioni attuali sono promettenti. Piccole modifiche possono ottimizzare ulteriormente.",
            "icon": "✅",
            "priority": "low",
        })

    # Sort by priority
    priority_order = {"critical": 0, "high": 1, "medium": 2, "low": 3}
    recommendations.sort(key=lambda r: priority_order.get(r.get("priority", "low"), 99))

    return {
        "overall": overall,
        "tumor_reduction_label": tumor_label,
        "healthy_loss_label": healthy_label,
        "time_to_recurrence_label": recurrence_label,
        "recommendations": recommendations,
        "has_recommendations": len(recommendations) > 0,
    }

is_staff = user_passes_test(lambda u: u.is_staff)

DEMO_MRN_PREFIX = "DEMO"


def can_edit_patient(user, patient: models.Patient) -> bool:
    if not getattr(user, "is_authenticated", False):
        return False
    if getattr(user, "is_staff", False):
        return True
    if patient.owner_id and patient.owner_id == user.id:
        return True
    if patient.owner_id is None and (patient.mrn or "").upper().startswith(DEMO_MRN_PREFIX):
        return True
    return False


@login_required
def dashboard(request: HttpRequest) -> HttpResponse:
    from django.db.models import Q, Max, OuterRef, Subquery
    from django.utils import timezone
    from datetime import timedelta
    
    patient_count = models.Patient.objects.count()
    
    # Create queryset for response counting (don't slice yet)
    recent_assessments_qs = (
        models.Assessment.objects.select_related("patient")
        .all()
        .order_by("-date")
    )
    
    # Response distribution from recent assessments (count before slicing)
    response_counts = {}
    for code, label in models.Assessment.RESPONSE_CHOICES:
        response_counts[code] = recent_assessments_qs.filter(response=code).count()
    
    # Now slice for display (latest 10)
    recent_assessments = recent_assessments_qs[:10]
    
    # R-ISS distribution from latest assessments (SQLite-compatible approach)
    # Get the latest assessment date for each patient
    latest_assessment_dates = models.Assessment.objects.filter(
        patient=OuterRef('patient')
    ).order_by('-date').values('date')[:1]
    
    # Get assessments that match the latest date for each patient
    latest_assessments = models.Assessment.objects.annotate(
        latest_date=Subquery(latest_assessment_dates)
    ).filter(date=F('latest_date'))
    
    riss_counts = {
        "I": latest_assessments.filter(r_iss="I").count(),
        "II": latest_assessments.filter(r_iss="II").count(),
        "III": latest_assessments.filter(r_iss="III").count(),
    }
    
    # Active therapies (no end date)
    active_therapies = models.PatientTherapy.objects.filter(
        end_date__isnull=True
    ).select_related('patient', 'regimen').count()
    
    # Patients needing attention (progressive disease or high FLC ratio)
    thirty_days_ago = timezone.now().date() - timedelta(days=30)
    patients_needing_attention = models.Patient.objects.filter(
        Q(assessments__response="PD", assessments__date__gte=thirty_days_ago) |
        Q(assessments__flc_ratio__gt=10, assessments__date__gte=thirty_days_ago)
    ).distinct().count()
    
    # High-risk cytogenetics count
    high_risk_codes = {"del(17p)", "t(4;14)", "t(14;16)", "1q+"}
    high_risk_patients = models.Patient.objects.filter(
        cytogenetics__abnormality__code__in=high_risk_codes
    ).distinct().count()
    
    context = {
        "patient_count": patient_count,
        "recent_assessments": recent_assessments,
        "riss_counts": riss_counts,
        "response_counts": response_counts,
        "active_therapies": active_therapies,
        "patients_needing_attention": patients_needing_attention,
        "high_risk_patients": high_risk_patients,
    }
    return render(request, "clinic/dashboard.html", context)


class PatientFilter(django_filters.FilterSet):
    r_iss = django_filters.ChoiceFilter(
        field_name="assessments__r_iss",
        label="R-ISS",
        choices=models.Assessment.R_ISS_CHOICES,
    )
    high_risk = django_filters.BooleanFilter(
        label="High-risk cytogenetics",
        method="filter_high_risk",
    )

    class Meta:
        model = models.Patient
        fields: list[str] = []

    def filter_high_risk(self, queryset, name, value):
        if value:
            high_risk_codes = {"del(17p)", "t(4;14)", "t(14;16)", "1q+"}
            return queryset.filter(
                cytogenetics__abnormality__code__in=high_risk_codes
            )
        return queryset


@login_required
def patient_list(request: HttpRequest) -> HttpResponse:
    query = request.GET.get("q", "").strip()
    qs = models.Patient.objects.all().prefetch_related("assessments")
    if query:
        qs = qs.filter(last_name__icontains=query)
    filterset = PatientFilter(request.GET or None, queryset=qs)
    if "r_iss" in filterset.form.fields:
        filterset.form.fields["r_iss"].widget.attrs.update({"class": "form-select"})
    if "high_risk" in filterset.form.fields:
        filterset.form.fields["high_risk"].widget.attrs.update({"class": "form-check-input"})
    patients = (
        filterset.qs.annotate(assessment_count=Count("assessments")).order_by("last_name")
    )
    context = {
        "patients": patients,
        "filterset": filterset,
        "query": query,
    }
    return render(request, "clinic/patient_list.html", context)


@login_required
def patient_detail(request: HttpRequest, pk: int) -> HttpResponse:
    patient = get_object_or_404(
        models.Patient.objects.prefetch_related(
            Prefetch("assessments", queryset=models.Assessment.objects.order_by("-date")),
            Prefetch(
                "cytogenetics",
                queryset=models.PatientCytogenetics.objects.select_related("abnormality").order_by("-detected_on"),
            ),
            Prefetch(
                "therapies",
                queryset=models.PatientTherapy.objects.select_related("regimen").order_by("-start_date"),
            ),
        ),
        pk=pk,
    )
    assessments = patient.assessments.all()
    latest_assessment = assessments.first()
    latest_assessment_id = latest_assessment.pk if latest_assessment else None
    latest_assessment_date_iso = (
        latest_assessment.date.isoformat() if latest_assessment and latest_assessment.date else ""
    )

    editable = can_edit_patient(request.user, patient)

    therapy_form = forms.PatientTherapyForm()
    if request.method == "POST":
        action = request.POST.get("action")
        if action == "add_therapy":
            if not editable:
                return HttpResponseForbidden("Not allowed")
            therapy_form = forms.PatientTherapyForm(request.POST)
            if therapy_form.is_valid():
                therapy = therapy_form.save(commit=False)
                therapy.patient = patient
                therapy.save()
                return redirect(reverse("clinic:patient_detail", args=[patient.id]))
    chart_points = _build_patient_detail_chart_points(patient, assessments)
    twin_chart_has_data = any(
        point.get("ldh") is not None or point.get("beta2m") is not None or point.get("riss") is not None
        for point in chart_points
    )

    # IMPORTANT: serialize as real JSON for safe embedding into JS (null instead of Python None).
    chart_points_json = json.dumps(chart_points)

    therapy_spans = [
        {
            "name": t.regimen.name,
            "start": t.start_date.isoformat() if t.start_date else None,
            "end": t.end_date.isoformat() if t.end_date else None,
        }
        for t in patient.therapies.all()
    ]

    therapy_spans_json = json.dumps(therapy_spans)

    # Observed effects: compare last assessment before therapy start vs last assessment during/after.
    # Use the already-prefetched assessments list to avoid N+1 queries.
    therapy_effects: dict[int, dict[str, object]] = {}
    try:
        from simulator.twin import build_patient_twin  # local import to avoid hard dependency at import time
    except Exception:  # pragma: no cover
        build_patient_twin = None

    # assessments is already ordered by -date (from prefetch); sort ascending for boundary lookups
    all_assessments_asc = sorted(assessments, key=lambda a: a.date)

    def maybe_float(v):
        try:
            return float(v) if v is not None else None
        except Exception:
            return None

    for t in patient.therapies.all():
        # baseline: latest assessment on or before therapy start
        baseline = None
        for a in reversed(all_assessments_asc):
            if a.date <= t.start_date:
                baseline = a
                break
        # follow: latest assessment on or before therapy end (or latest overall)
        follow_until = t.end_date
        follow = None
        for a in reversed(all_assessments_asc):
            if follow_until is None or a.date <= follow_until:
                follow = a
                break

        baseline_m = maybe_float(getattr(baseline, "m_protein_g_dl", None)) if baseline else None
        follow_m = maybe_float(getattr(follow, "m_protein_g_dl", None)) if follow else None
        delta_m = (follow_m - baseline_m) if (baseline_m is not None and follow_m is not None) else None

        baseline_risk = None
        follow_risk = None
        if build_patient_twin and baseline:
            baseline_risk = maybe_float(build_patient_twin(baseline).get("risk_score"))
        if build_patient_twin and follow:
            follow_risk = maybe_float(build_patient_twin(follow).get("risk_score"))
        delta_risk = (
            (follow_risk - baseline_risk)
            if (baseline_risk is not None and follow_risk is not None)
            else None
        )

        therapy_effects[t.pk] = {
            "baseline_date": baseline.date.isoformat() if baseline else None,
            "follow_date": follow.date.isoformat() if follow else None,
            "delta_m": delta_m,
            "baseline_risk": baseline_risk,
            "follow_risk": follow_risk,
            "delta_risk": delta_risk,
        }

    therapy_rows = [
        {"therapy": t, "eff": therapy_effects.get(t.pk)}
        for t in patient.therapies.all()
    ]

    regimen_count = models.Regimen.objects.count()
    regimen_add_url = reverse("clinic:regimen_new")
    regimen_list_url = reverse("clinic:regimen_list")

    # Optional: provide a 1-click path into Simulator for beginners.
    quickstart_scenario_pk = None
    try:
        from simulator.models import Scenario

        quickstart_scenario_pk = Scenario.objects.filter(active=True).order_by("pk").values_list("pk", flat=True).first()
    except Exception:
        quickstart_scenario_pk = None

    # Optional: surface the most recent Simulator results for this patient (based on latest Assessment snapshot).
    latest_simulation_attempt = None
    latest_simulation_summary = None
    latest_simulation_artifacts = None
    latest_simulation_scenario_url = None
    latest_simulation_interpretation = None
    if latest_assessment_id:
        try:
            from django.db.models import Q
            from simulator.models import SimulationAttempt

            aid = latest_assessment_id
            latest_simulation_attempt = (
                SimulationAttempt.objects.select_related("scenario")
                .filter(
                    Q(parameters__twin_assessment_id=aid)
                    | Q(parameters__twin_assessment_id=str(aid))
                    | Q(parameters__assessment_id=aid)
                    | Q(parameters__assessment_id=str(aid))
                )
                .order_by("-submitted")
                .first()
            )

            # Fallback: if no run exists for the *latest* snapshot, show the most recent run
            # for any snapshot belonging to this patient.
            if not latest_simulation_attempt:
                patient_assessment_ids = [a.pk for a in assessments if a.pk is not None]
                patient_assessment_ids_str = [str(aid) for aid in patient_assessment_ids]
                if patient_assessment_ids:
                    latest_simulation_attempt = (
                        SimulationAttempt.objects.select_related("scenario")
                        .filter(
                            Q(parameters__twin_assessment_id__in=patient_assessment_ids)
                            | Q(parameters__twin_assessment_id__in=patient_assessment_ids_str)
                            | Q(parameters__assessment_id__in=patient_assessment_ids)
                            | Q(parameters__assessment_id__in=patient_assessment_ids_str)
                        )
                        .order_by("-submitted")
                        .first()
                    )

            if latest_simulation_attempt:
                latest_simulation_summary = latest_simulation_attempt.results_summary or None
                latest_simulation_artifacts = latest_simulation_attempt.artifacts or None
                latest_simulation_interpretation = _interpret_latest_simulation(
                    latest_simulation_summary if isinstance(latest_simulation_summary, dict) else None,
                    latest_simulation_attempt.parameters if isinstance(latest_simulation_attempt.parameters, dict) else None,
                )
                latest_simulation_scenario_url = (
                    reverse("simulator:scenario_detail", args=[latest_simulation_attempt.scenario_id])
                    + f"?twin_assessment_id={aid}#simulationResults"
                )

                plot_url = None
                if isinstance(latest_simulation_artifacts, dict):
                    plot_url = latest_simulation_artifacts.get("plot")
                embed_debug_logger.info(
                    "patient_detail patient_id=%s assessment_id=%s attempt_id=%s scenario_id=%s plot_url=%r page_path=%s",
                    patient.id,
                    aid,
                    latest_simulation_attempt.id,
                    latest_simulation_attempt.scenario_id,
                    plot_url,
                    request.get_full_path(),
                )
        except Exception:
            latest_simulation_attempt = None

    research_twin_state = None
    research_latest_residual = None
    research_last_calibration_status = "Not initialized"
    research_missing_required_data: list[str] = ["assessment"] if latest_assessment is None else []
    research_counterfactual_count = 0
    research_completed_counterfactual_count = 0
    research_provenance_count = 0
    research_simple_url = reverse("twin_engine:simple_research_view", args=[patient.id])
    research_cockpit_url = reverse("twin_engine:research_cockpit", args=[patient.id])
    research_developer_url = reverse("twin_engine:developer_console") + f"?patient_id={patient.id}"
    research_glossary_url = reverse("twin_engine:research_glossary")
    research_twin_url = research_cockpit_url
    interpretation_lab_count = len(assessments)
    interpretation_lab_date_range = _format_date_range([assessment.date for assessment in assessments])
    interpretation_therapy_count = len(patient.therapies.all())
    interpretation_therapy_date_range = _format_date_range(
        [therapy.start_date for therapy in patient.therapies.all()] + [therapy.end_date for therapy in patient.therapies.all()]
    )
    interpretation_simulation_count = 0
    interpretation_adverse_event_count = 0
    interpretation_adverse_event_categories: list[str] = []
    try:
        from twin_engine.models import AdverseEvent, CounterfactualRun, ObservationResidual, SimulationRunMetadata
        from twin_engine.state_model import get_current_twin_state
        from twin_engine.validators import validate_assessment_minimum_fields
        from simulator.models import SimulationAttempt
        from django.db.models import Q

        research_twin_state = get_current_twin_state(patient)
        if latest_assessment is not None:
            try:
                research_missing_required_data = validate_assessment_minimum_fields(latest_assessment)["missing"]
            except Exception as exc:
                research_missing_required_data = [str(exc)]
        if research_twin_state is not None:
            research_latest_residual = (
                ObservationResidual.objects.filter(twin_state=research_twin_state)
                .select_related("assessment")
                .order_by("-created_at")
                .first()
            )
            research_last_calibration_status = research_twin_state.get_method_display()
            research_provenance_count = SimulationRunMetadata.objects.filter(twin_state=research_twin_state).count()
        research_counterfactual_count = patient.counterfactual_runs.count()
        research_completed_counterfactual_count = patient.counterfactual_runs.filter(
            status=CounterfactualRun.STATUS_COMPLETED
        ).count()

        patient_assessment_ids = [assessment.pk for assessment in assessments if assessment.pk is not None]
        patient_assessment_ids_str = [str(assessment_id) for assessment_id in patient_assessment_ids]
        if patient_assessment_ids:
            interpretation_simulation_count = SimulationAttempt.objects.filter(
                Q(parameters__twin_assessment_id__in=patient_assessment_ids)
                | Q(parameters__twin_assessment_id__in=patient_assessment_ids_str)
                | Q(parameters__assessment_id__in=patient_assessment_ids)
                | Q(parameters__assessment_id__in=patient_assessment_ids_str)
            ).count()

        adverse_events = list(patient.adverse_events.order_by("-date"))
        interpretation_adverse_event_count = len(adverse_events)
        adverse_event_labels = dict(AdverseEvent.EVENT_TYPE_CHOICES)
        interpretation_adverse_event_categories = [
            adverse_event_labels.get(event_type, event_type)
            for event_type, _count in Counter(event.event_type for event in adverse_events).most_common(3)
        ]
    except Exception:
        research_twin_state = None

    if research_completed_counterfactual_count == 0:
        interpretation_current_conclusion = "No treatment-comparison conclusion available."
    else:
        interpretation_current_conclusion = (
            f"{research_completed_counterfactual_count} exploratory what-if scenario"
            f"{'s' if research_completed_counterfactual_count != 1 else ''} available for review in Simple Research View or Scientific Cockpit."
        )

    interpretation_data_available = []
    if interpretation_lab_count:
        interpretation_data_available.append(
            f"Lab snapshots: {interpretation_lab_count} ({interpretation_lab_date_range})."
        )
    else:
        interpretation_data_available.append("Lab snapshots: none recorded.")
    if interpretation_therapy_count:
        interpretation_data_available.append(
            f"Therapies: {interpretation_therapy_count} ({interpretation_therapy_date_range})."
        )
    else:
        interpretation_data_available.append("Therapies: none recorded.")
    if interpretation_adverse_event_count:
        interpretation_data_available.append(
            f"Adverse events: {interpretation_adverse_event_count} ({', '.join(interpretation_adverse_event_categories)})."
        )
    else:
        interpretation_data_available.append("Adverse events: no structured categories recorded.")

    if interpretation_simulation_count or research_completed_counterfactual_count:
        interpretation_run_status = (
            f"Mechanistic simulation records: {interpretation_simulation_count}. "
            f"Exploratory what-if runs completed: {research_completed_counterfactual_count}."
        )
    else:
        interpretation_run_status = "No completed mechanistic simulation or exploratory what-if comparison is currently recorded for this patient."

    if research_missing_required_data:
        research_recommended_next_step = {
            "label": "Complete minimum research input fields",
            "detail": "The latest assessment is missing fields required for a reproducible twin initialization.",
            "href": reverse("clinic:assessment_new", args=[patient.id]) if editable else research_twin_url,
        }
    elif research_twin_state is None:
        research_recommended_next_step = {
            "label": "Initialize research twin",
            "detail": "Create a PatientTwinState from a dated assessment before calibration or what-if simulation.",
            "href": research_twin_url + "#twin-state",
        }
    elif research_latest_residual is None:
        research_recommended_next_step = {
            "label": "Run or review calibration",
            "detail": "Calibration residuals are needed to understand how well the model reproduces observed biomarkers.",
            "href": research_twin_url + "#calibration-quality",
        }
    elif research_completed_counterfactual_count == 0:
        research_recommended_next_step = {
            "label": "Run predefined what-if scenarios",
            "detail": "Completed mechanistic runs are needed before trajectory and utility comparisons are meaningful.",
            "href": research_twin_url + "#what-if-scenarios",
        }
    elif research_provenance_count == 0:
        research_recommended_next_step = {
            "label": "Regenerate provenance metadata",
            "detail": "Traceability is incomplete until simulation metadata records exist for the current twin state.",
            "href": research_twin_url + "#provenance",
        }
    else:
        research_recommended_next_step = {
            "label": "Review full research cockpit",
            "detail": "Data, twin state, calibration, scenarios, and provenance are present; review limitations before interpreting outputs.",
            "href": research_twin_url,
        }

    research_interpretation_status = {
        "headline": "Research interpretation status",
        "can_answer": "No. This page cannot determine whether the patient could have been treated better.",
        "evidence_required": [
            "A defined utility function, for example U(a)=benefit(a)-toxicity(a)-uncertainty(a).",
            "Comparable exploratory what-if runs with documented assumptions, calibration, and uncertainty diagnostics.",
            "Structured disease, therapy, and adverse-event data across time for the same patient context.",
        ],
        "data_available": interpretation_data_available,
        "run_status": interpretation_run_status,
        "current_conclusion": interpretation_current_conclusion,
        "can_conclude": (
            "The platform can compare mechanistic simulation outputs and exploratory what-if branches, "
            "but only as research interpretation and not clinical recommendation."
        ),
        "cannot_conclude": (
            "The platform cannot prove what would have happened clinically, cannot establish causal effect, "
            "and cannot prove that one therapy would have produced a better real-world outcome on this page."
        ),
        "utility_formula": "U(a)=benefit(a)-toxicity(a)-uncertainty(a)",
        "utility_explanation": "Better treatment is undefined until a utility function is specified.",
        "plain_language": "The platform can compare simulated strategies, but cannot prove what would have happened clinically.",
        "next_action": research_recommended_next_step,
    }

    context = {
        "patient": patient,
        "assessments": assessments,
        "latest_assessment_id": latest_assessment_id,
        "latest_assessment_date_iso": latest_assessment_date_iso,
        "cytogenetics": patient.cytogenetics.all(),
        "therapies": patient.therapies.all(),
        "assessment_form": forms.AssessmentForm(),
        "therapy_form": therapy_form,
        "chart_points": chart_points,
        "twin_chart_has_data": twin_chart_has_data,
        "therapy_spans": therapy_spans,
        "chart_points_json": chart_points_json,
        "therapy_spans_json": therapy_spans_json,
        "can_edit_patient": editable,
        "regimen_count": regimen_count,
        "regimen_add_url": regimen_add_url,
        "regimen_list_url": regimen_list_url,
        "therapy_rows": therapy_rows,
        "quickstart_scenario_pk": quickstart_scenario_pk,
        "latest_simulation_attempt": latest_simulation_attempt,
        "latest_simulation_summary": latest_simulation_summary,
        "latest_simulation_artifacts": latest_simulation_artifacts,
        "latest_simulation_scenario_url": latest_simulation_scenario_url,
        "latest_simulation_interpretation": latest_simulation_interpretation,
        "research_twin_state": research_twin_state,
        "research_latest_residual": research_latest_residual,
        "research_last_calibration_status": research_last_calibration_status,
        "research_missing_required_data": research_missing_required_data,
        "research_counterfactual_count": research_counterfactual_count,
        "research_completed_counterfactual_count": research_completed_counterfactual_count,
        "research_provenance_count": research_provenance_count,
        "research_simple_url": research_simple_url,
        "research_cockpit_url": research_cockpit_url,
        "research_developer_url": research_developer_url,
        "research_glossary_url": research_glossary_url,
        "research_twin_url": research_twin_url,
        "research_recommended_next_step": research_recommended_next_step,
        "research_interpretation_status": research_interpretation_status,
    }
    return render(request, "clinic/patient_detail.html", context)


def _format_date_range(values: list[date | None]) -> str:
    available = sorted(value for value in values if value is not None)
    if not available:
        return "no dated records"
    if available[0] == available[-1]:
        return available[0].isoformat()
    return f"{available[0].isoformat()} to {available[-1].isoformat()}"


def _build_patient_detail_chart_points(patient: models.Patient, assessments) -> list[dict[str, object]]:
    points_by_date: dict[str, dict[str, object]] = {}
    for assessment in reversed(assessments):
        date_key = assessment.date.isoformat()
        points_by_date[date_key] = {
            "date": date_key,
            "m": float(assessment.m_protein_g_dl) if assessment.m_protein_g_dl is not None else None,
            "flc": float(assessment.flc_ratio) if assessment.flc_ratio is not None else None,
            "ldh": float(assessment.ldH_u_l) if assessment.ldH_u_l is not None else None,
            "beta2m": float(assessment.beta2m_mg_l) if assessment.beta2m_mg_l is not None else None,
            "riss": (1 if assessment.r_iss == "I" else 2 if assessment.r_iss == "II" else 3 if assessment.r_iss == "III" else None),
        }
    try:
        from twin_engine.models import LongitudinalLabResult
    except Exception:
        return list(points_by_date.values())

    analyte_to_key = {
        LongitudinalLabResult.ANALYTE_M_PROTEIN: "m",
        LongitudinalLabResult.ANALYTE_FLC_RATIO: "flc",
        LongitudinalLabResult.ANALYTE_LDH: "ldh",
        LongitudinalLabResult.ANALYTE_BETA2M: "beta2m",
    }
    labs = LongitudinalLabResult.objects.filter(patient=patient, analyte__in=analyte_to_key).order_by("date", "id")
    for lab in labs:
        if lab.value is None:
            continue
        date_key = lab.date.isoformat()
        points_by_date.setdefault(
            date_key,
            {"date": date_key, "m": None, "flc": None, "ldh": None, "beta2m": None, "riss": None},
        )
        points_by_date[date_key][analyte_to_key[lab.analyte]] = float(lab.value)
    return [points_by_date[key] for key in sorted(points_by_date)]


@login_required
def regimen_list(request: HttpRequest) -> HttpResponse:
    regimens = models.Regimen.objects.all().order_by("name")
    context = {
        "regimens": regimens,
    }
    return render(request, "clinic/regimen_list.html", context)


@login_required
def regimen_new(request: HttpRequest) -> HttpResponse:
    if request.method == "POST":
        form = forms.RegimenForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect(reverse("clinic:regimen_list"))
    else:
        form = forms.RegimenForm()
    return render(request, "clinic/regimen_form.html", {"form": form})


@login_required
@is_staff
def patient_new(request: HttpRequest) -> HttpResponse:
    """Create a new patient."""
    if request.method == "POST":
        form = forms.PatientForm(request.POST)
        if form.is_valid():
            patient = form.save()
            return redirect(reverse("clinic:patient_detail", args=[patient.id]))
    else:
        form = forms.PatientForm()
    
    context = {
        "form": form,
        "title": "Add New Patient",
    }
    return render(request, "clinic/patient_form.html", context)


@login_required
@is_staff
def patient_edit(request: HttpRequest, pk: int) -> HttpResponse:
    """Edit an existing patient."""
    patient = get_object_or_404(models.Patient, pk=pk)
    
    if request.method == "POST":
        form = forms.PatientForm(request.POST, instance=patient)
        if form.is_valid():
            patient = form.save()
            return redirect(reverse("clinic:patient_detail", args=[patient.id]))
    else:
        form = forms.PatientForm(instance=patient)
    
    context = {
        "form": form,
        "patient": patient,
        "title": f"Edit Patient: {patient}",
    }
    return render(request, "clinic/patient_form.html", context)


@login_required
def assessment_new(request: HttpRequest, patient_id: int) -> HttpResponse:
    patient = get_object_or_404(models.Patient, pk=patient_id)
    if not can_edit_patient(request.user, patient):
        return HttpResponseForbidden("Not allowed")
    if request.method == "POST":
        form = forms.AssessmentForm(request.POST)
        if form.is_valid():
            assessment = form.save(commit=False)
            assessment.patient = patient
            assessment.save()
            return redirect(reverse("clinic:patient_detail", args=[patient.id]))
    else:
        form = forms.AssessmentForm()
    context = {
        "patient": patient,
        "form": form,
    }
    return render(request, "clinic/assessment_form.html", context)


# ══════════════════════════════════════════════════════════════════════════════
# SYMPTOM ASSESSMENT VIEWS
# ══════════════════════════════════════════════════════════════════════════════

@login_required
def symptom_assessment_new(request: HttpRequest, patient_id: int) -> HttpResponse:
    """Create a new symptom assessment for a patient."""
    patient = get_object_or_404(models.Patient, pk=patient_id)
    if not can_edit_patient(request.user, patient):
        return HttpResponseForbidden("Not allowed")
    
    if request.method == "POST":
        form = SymptomAssessmentForm(request.POST)
        if form.is_valid():
            assessment = form.save(commit=False)
            assessment.patient = patient
            assessment.save()
            return redirect(reverse("clinic:patient_detail", args=[patient.id]))
    else:
        form = SymptomAssessmentForm()
    
    context = {
        "patient": patient,
        "form": form,
    }
    return render(request, "clinic/symptom_assessment_form.html", context)


@login_required
def symptom_assessment_list(request: HttpRequest, patient_id: int) -> HttpResponse:
    """List all symptom assessments for a patient."""
    patient = get_object_or_404(models.Patient, pk=patient_id)
    assessments = SymptomAssessment.objects.filter(patient=patient).order_by("-assessment_date")
    
    context = {
        "patient": patient,
        "assessments": assessments,
    }
    return render(request, "clinic/symptom_assessment_list.html", context)


# ══════════════════════════════════════════════════════════════════════════════
# PROGNOSIS & TIMELINE VIEWS
# ══════════════════════════════════════════════════════════════════════════════

@login_required
def prognosis_timeline(request: HttpRequest, patient_id: int) -> HttpResponse:
    """Display exploratory prognosis context and timeline for a patient."""
    from simulator.prognosis import estimate_prognosis
    
    patient = get_object_or_404(models.Patient, pk=patient_id)
    
    # Get latest assessment for R-ISS and response
    latest_assessment = patient.assessments.first()
    
    # Get cytogenetics
    cytogenetics = list(
        patient.cytogenetics.values_list("abnormality__code", flat=True)
    )
    
    # Get latest symptom assessment for ECOG
    latest_symptoms = SymptomAssessment.objects.filter(patient=patient).first()
    ecog = latest_symptoms.ecog_status if latest_symptoms else None
    
    # Count therapy lines
    therapy_count = patient.therapies.count()
    line_of_therapy = max(therapy_count, 1)
    structured_r_iss = latest_assessment.r_iss if latest_assessment and latest_assessment.r_iss else None
    response_value = latest_assessment.response if latest_assessment and latest_assessment.response else None
    response_label = latest_assessment.get_response_display() if response_value else None
    
    # Build patient parameters
    patient_params = {
        "r_iss": structured_r_iss or "II",
        "cytogenetics": cytogenetics,
        "age": patient.age,
        "ecog": ecog,
        "response": response_value,
        "line_of_therapy": line_of_therapy,
    }
    
    # Calculate prognosis estimate
    estimate = estimate_prognosis(**patient_params)
    confidence_score = float(estimate.confidence or 0.0)
    confidence_percentage = round(confidence_score * 100)
    if confidence_score <= 0.10:
        confidence_warning = "Model confidence is extremely limited. Treat this page as exploratory context only."
        confidence_alert_class = "alert-danger"
    elif confidence_score < 0.50:
        confidence_warning = "Model confidence is limited. Do not use this page for treatment comparison."
        confidence_alert_class = "alert-warning"
    else:
        confidence_warning = "Model confidence is still heuristic and does not establish individual prognosis."
        confidence_alert_class = "alert-info"
    
    # Calculate intermediate timepoints (3m, 6m) by interpolation
    import math
    def survival_prob(median: float, months: float) -> float:
        if median <= 0:
            return 0.0
        lambda_rate = math.log(2) / median
        return math.exp(-lambda_rate * months) * 100
    
    pfs_3m = survival_prob(estimate.median_pfs_months, 3)
    pfs_6m = survival_prob(estimate.median_pfs_months, 6)
    timeline_rows = [
        {
            "horizon": horizon,
            "pfs_probability": f"{probability}%",
            "interpretation": _build_prognosis_timeline_interpretation(probability),
            "caveat": "Population/heuristic estimate; not an individual prediction.",
        }
        for horizon, probability in [
            ("3 months", round(pfs_3m)),
            ("6 months", round(pfs_6m)),
            ("12 months", round(estimate.pfs_12m * 100)),
            ("24 months", round(estimate.pfs_24m * 100)),
            ("36 months", round(estimate.pfs_36m * 100)),
        ]
    ]

    adjustment_rows = _build_prognosis_adjustment_rows(estimate.modifiers_applied)
    risk_group_label = _format_risk_category_label(estimate.risk_category)
    patient_detail_url = reverse("clinic:patient_detail", args=[patient.id])
    research_simple_url = reverse("twin_engine:simple_research_view", args=[patient.id])
    research_cockpit_url = reverse("twin_engine:research_cockpit", args=[patient.id])
    simulation_url = reverse("simulator:scenario_list")
    if latest_assessment is not None:
        simulation_url = f"{simulation_url}?twin_assessment_id={latest_assessment.pk}"

    data_source_rows = [
        {
            "variable": "R-ISS",
            "value": structured_r_iss or "II (default fallback)",
            "classification": "RAW STRUCTURED" if structured_r_iss else "UNKNOWN",
            "role": (
                "Structured stage used to select the literature-derived baseline survival values."
                if structured_r_iss
                else "Fallback stage used because no structured R-ISS value was recorded."
            ),
        },
        {
            "variable": "Age",
            "value": f"{patient.age} years",
            "classification": "DERIVED",
            "role": "Derived from birth date and used as a heuristic age modifier.",
        },
        {
            "variable": "Line of therapy",
            "value": str(line_of_therapy) if therapy_count else "1 (default from no recorded therapies)",
            "classification": "DERIVED" if therapy_count else "HEURISTIC",
            "role": (
                "Derived from recorded therapy count and used as relapse context."
                if therapy_count
                else "Default context when no explicit therapy line is recorded."
            ),
        },
        {
            "variable": "Response / progression context",
            "value": response_label or "Not recorded",
            "classification": "RAW STRUCTURED" if response_label else "UNKNOWN",
            "role": (
                "Structured response context that modifies the heuristic estimate."
                if response_label
                else "No structured response context was available for this estimate."
            ),
        },
        {
            "variable": "Adjustments applied",
            "value": "; ".join(estimate.modifiers_applied) if estimate.modifiers_applied else "No explicit adjustment beyond baseline stage.",
            "classification": "HEURISTIC",
            "role": "Rule-based modifiers that shift or contextualize the baseline estimate.",
        },
        {
            "variable": "Reference",
            "value": estimate.reference,
            "classification": "LITERATURE-BASED",
            "role": "Literature source used for the baseline numerical survival values in this estimate.",
        },
        {
            "variable": "Median PFS",
            "value": f"{estimate.median_pfs_months:.1f} months",
            "classification": "HEURISTIC",
            "role": "Adjusted exploratory progression-free survival estimate using literature baseline values plus rule-based modifiers.",
        },
        {
            "variable": "Median OS",
            "value": f"{estimate.median_os_months:.1f} months",
            "classification": "HEURISTIC",
            "role": "Adjusted exploratory overall survival estimate using literature baseline values plus rule-based modifiers.",
        },
    ]

    estimate_cards = [
        {
            "label": "Risk group",
            "value": risk_group_label,
            "detail": "Derived from the literature-based baseline stage plus applied heuristic modifiers.",
        },
        {
            "label": "Median PFS",
            "value": f"{estimate.median_pfs_months:.1f} months",
            "detail": "Exploratory progression-free survival estimate.",
        },
        {
            "label": "Median OS",
            "value": f"{estimate.median_os_months:.1f} months",
            "detail": "Exploratory overall survival estimate.",
        },
        {
            "label": "Confidence / completeness score",
            "value": f"{confidence_percentage}%",
            "detail": "Input completeness and heuristic reliability score, not clinical certainty.",
        },
    ]
    
    context = {
        "patient": patient,
        "page_purpose": "This page organizes exploratory prognosis context from structured patient factors, literature-derived baseline survival values, and heuristic modifiers.",
        "interpretation_status": {
            "title": "Prognosis / timeline interpretation status",
            "summary": "This page shows heuristic/literature-informed prognosis estimates, not a patient-specific clinical prediction.",
            "limit": "This page cannot determine whether an alternative treatment would have changed this patient’s outcome.",
        },
        "confidence_panel": {
            "title": "Model confidence / uncertainty",
            "score_label": "Model completeness/confidence score",
            "percentage": confidence_percentage,
            "explanation": "This score reflects input completeness and heuristic reliability, not clinical certainty.",
            "warning": confidence_warning,
            "alert_class": confidence_alert_class,
        },
        "data_source_rows": data_source_rows,
        "estimate_cards": estimate_cards,
        "timeline_rows": timeline_rows,
        "adjustment_rows": adjustment_rows,
        "reference_info": {
            "title": "Reference",
            "label": "Literature source used for baseline numerical values",
            "value": estimate.reference,
            "detail": "The current code uses literature-derived baseline survival values and then applies rule-based heuristic modifiers to generate the displayed estimate.",
        },
        "what_can_be_concluded": [
            "The page can provide exploratory prognosis context.",
            "The page can show which patient factors modify the heuristic estimate.",
            "The page can highlight uncertainty and missing evidence.",
        ],
        "what_cannot_be_concluded": [
            "It cannot prove individual survival.",
            "It cannot determine treatment superiority.",
            "It cannot infer what would have happened under an alternative therapy.",
            "It cannot replace counterfactual simulation or clinical validation.",
        ],
        "next_actions": {
            "primary": {
                "label": "Open Simple Research View",
                "href": research_simple_url,
            },
            "secondary": [
                {
                    "label": "Back to Patient",
                    "href": patient_detail_url,
                },
                {
                    "label": "Open Scientific Cockpit",
                    "href": research_cockpit_url,
                },
                {
                    "label": "Start exploratory simulation",
                    "href": simulation_url,
                },
            ],
        },
    }
    return render(request, "clinic/prognosis_timeline.html", context)


def _build_prognosis_timeline_interpretation(probability_pct: int) -> str:
    if probability_pct >= 80:
        return "Higher heuristic progression-free probability at this horizon."
    if probability_pct >= 60:
        return "Moderate heuristic progression-free probability at this horizon."
    if probability_pct >= 40:
        return "Meaningful uncertainty about remaining progression-free by this horizon."
    return "Lower heuristic progression-free probability at this horizon."


def _build_prognosis_adjustment_rows(modifiers_applied: list[str]) -> list[dict[str, str]]:
    if not modifiers_applied:
        return [
            {
                "trigger": "No explicit adjustment recorded beyond the baseline stage lookup.",
                "source_type": "HEURISTIC",
                "effect_direction": "No extra direction recorded beyond the baseline literature lookup.",
                "estimate_effect": "Context only",
            }
        ]

    rows: list[dict[str, str]] = []
    for modifier in modifiers_applied:
        modifier_lower = modifier.lower()
        source_type = "HEURISTIC"
        effect_direction = "Rule-based modifier applied to the estimate."
        estimate_effect = "Changes PFS/OS estimate"

        if modifier_lower.startswith("cytogenetics:"):
            source_type = "RAW STRUCTURED"
            effect_direction = "Higher-risk cytogenetics usually shift the estimate downward."
        elif modifier_lower.startswith("age"):
            source_type = "DERIVED"
            effect_direction = "Age band changes the estimate through a stored heuristic modifier."
        elif modifier_lower.startswith("ecog"):
            source_type = "RAW STRUCTURED"
            effect_direction = "Worse functional status shifts the estimate downward."
        elif modifier_lower.startswith("response:"):
            source_type = "RAW STRUCTURED"
            effect_direction = _build_prognosis_response_effect_direction(modifier_lower)
        elif "mrd" in modifier_lower:
            source_type = "RAW STRUCTURED"
            effect_direction = "MRD status changes the estimate according to the stored heuristic modifier."
        elif modifier_lower.startswith("line"):
            source_type = "DERIVED"
            effect_direction = "Later therapy lines shift the estimate downward."

        rows.append(
            {
                "trigger": modifier,
                "source_type": source_type,
                "effect_direction": effect_direction,
                "estimate_effect": estimate_effect,
            }
        )
    return rows


def _build_prognosis_response_effect_direction(modifier_lower: str) -> str:
    if "stringent complete response" in modifier_lower or "complete response" in modifier_lower or "very good partial response" in modifier_lower:
        return "Deeper recorded response shifts the estimate upward relative to partial response."
    if "stable disease" in modifier_lower or "progressive disease" in modifier_lower:
        return "Less favorable recorded response shifts the estimate downward relative to partial response."
    return "Recorded response modifies the estimate relative to the stored baseline response assumption."


def _format_risk_category_label(risk_category: str) -> str:
    labels = {
        "standard": "Standard",
        "intermediate": "Intermediate",
        "high": "High",
        "very_high": "Very High",
    }
    return labels.get(risk_category, risk_category.replace("_", " ").title())


@login_required
@require_http_methods(["GET"])
def prognosis_api(request: HttpRequest, patient_id: int) -> JsonResponse:
    """JSON API for prognosis data."""
    from simulator.prognosis import estimate_prognosis
    
    patient = get_object_or_404(models.Patient, pk=patient_id)
    if not can_edit_patient(request.user, patient):
        return JsonResponse({"error": "forbidden"}, status=403)
    
    # Get parameters from request or patient data
    latest_assessment = patient.assessments.first()
    cytogenetics = list(patient.cytogenetics.values_list("abnormality__code", flat=True))
    latest_symptoms = SymptomAssessment.objects.filter(patient=patient).first()
    
    params = {
        "r_iss": request.GET.get("r_iss") or (latest_assessment.r_iss if latest_assessment else "II"),
        "cytogenetics": request.GET.getlist("cytogenetics") or cytogenetics,
        "age": int(request.GET.get("age", patient.age)),
        "ecog": int(request.GET.get("ecog")) if request.GET.get("ecog") else (latest_symptoms.ecog_status if latest_symptoms else None),
        "response": request.GET.get("response") or (latest_assessment.response if latest_assessment else None),
        "mrd_status": request.GET.get("mrd_status"),
        "line_of_therapy": int(request.GET.get("line", 1)),
    }
    
    estimate = estimate_prognosis(**params)
    
    return JsonResponse({
        "patient_id": patient_id,
        "parameters": params,
        "estimate": estimate.to_dict(),
    })


# ══════════════════════════════════════════════════════════════════════════════
# REGIMEN SUGGESTER VIEWS
# ══════════════════════════════════════════════════════════════════════════════

@login_required
def regimen_suggestions(request: HttpRequest, patient_id: int) -> HttpResponse:
    """Display treatment regimen suggestions for a patient."""
    from simulator.regimen_suggester import suggest_regimens
    
    patient = get_object_or_404(models.Patient, pk=patient_id)
    
    # Get patient characteristics
    latest_assessment = patient.assessments.first()
    latest_symptoms = SymptomAssessment.objects.filter(patient=patient).first()
    
    # Determine transplant eligibility based on age
    age = patient.age
    transplant_eligible = age < 70  # Simplified heuristic
    
    # Get cytogenetics risk
    high_risk_cytos = ["del(17p)", "t(4;14)", "t(14;16)", "1q21"]
    patient_cytos = list(patient.cytogenetics.values_list("abnormality__code", flat=True))
    has_high_risk = any(
        any(hr.lower() in pc.lower() for hr in high_risk_cytos) 
        for pc in patient_cytos
    )
    
    # Get prior therapies
    prior_therapies = list(
        patient.therapies.values_list("regimen__components", flat=True)
    )
    prior_agents = []
    for pt in prior_therapies:
        if pt:
            prior_agents.extend([a.strip() for a in pt.split(",")])
    
    # Get neuropathy grade from symptoms
    neuropathy = 0
    if latest_symptoms:
        neuropathy = latest_symptoms.max_neuropathy_grade
    
    # Count therapy lines
    line_of_therapy = max(patient.therapies.count(), 1)
    
    # Get suggestions
    suggestions = suggest_regimens(
        age=age,
        transplant_eligible=transplant_eligible,
        ecog=latest_symptoms.ecog_status if latest_symptoms else 0,
        r_iss=latest_assessment.r_iss if latest_assessment else None,
        high_risk_cytogenetics=has_high_risk,
        line_of_therapy=line_of_therapy,
        prior_therapies=prior_agents,
        neuropathy_grade=neuropathy,
    )
    
    context = {
        "patient": patient,
        "suggestions": suggestions,
        "patient_info": {
            "age": age,
            "transplant_eligible": transplant_eligible,
            "high_risk_cytogenetics": has_high_risk,
            "line_of_therapy": line_of_therapy,
            "neuropathy_grade": neuropathy,
        },
    }
    return render(request, "clinic/regimen_suggestions.html", context)


@login_required
@require_http_methods(["GET"])
def regimen_suggestions_api(request: HttpRequest, patient_id: int) -> JsonResponse:
    """JSON API for regimen suggestions."""
    from simulator.regimen_suggester import suggest_regimens
    
    patient = get_object_or_404(models.Patient, pk=patient_id)
    if not can_edit_patient(request.user, patient):
        return JsonResponse({"error": "forbidden"}, status=403)
    
    # Parse parameters from request
    age = int(request.GET.get("age", patient.age))
    transplant_eligible = request.GET.get("transplant_eligible", "").lower() == "true"
    ecog = int(request.GET.get("ecog", 0))
    r_iss = request.GET.get("r_iss")
    high_risk = request.GET.get("high_risk", "").lower() == "true"
    line = int(request.GET.get("line", 1))
    neuropathy = int(request.GET.get("neuropathy", 0))
    
    suggestions = suggest_regimens(
        age=age,
        transplant_eligible=transplant_eligible,
        ecog=ecog,
        r_iss=r_iss,
        high_risk_cytogenetics=high_risk,
        line_of_therapy=line,
        neuropathy_grade=neuropathy,
    )
    
    return JsonResponse({
        "patient_id": patient_id,
        "suggestions": suggestions,
    })


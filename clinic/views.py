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

from mmportal.governance import (
    EpistemicLabel,
    build_model_relative_diagnostics,
    governance_metadata,
)

from . import forms, models
from .models_symptoms import SymptomAssessment
from .forms_symptoms import SymptomAssessmentForm, QuickSymptomForm

embed_debug_logger = logging.getLogger("embed_debug")


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
    latest_simulation_diagnostics = None
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
                latest_simulation_diagnostics = build_model_relative_diagnostics(
                    latest_simulation_summary if isinstance(latest_simulation_summary, dict) else None,
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
        research_next_step = {
            "label": "Complete minimum research input fields",
            "detail": "The latest assessment is missing fields required for a reproducible twin initialization.",
            "href": reverse("clinic:assessment_new", args=[patient.id]) if editable else research_twin_url,
        }
    elif research_twin_state is None:
        research_next_step = {
            "label": "Initialize research twin",
            "detail": "Create a PatientTwinState from a dated assessment before calibration or what-if simulation.",
            "href": research_twin_url + "#twin-state",
        }
    elif research_latest_residual is None:
        research_next_step = {
            "label": "Run or review calibration",
            "detail": "Calibration residuals are needed to understand how well the model reproduces observed biomarkers.",
            "href": research_twin_url + "#calibration-quality",
        }
    elif research_completed_counterfactual_count == 0:
        research_next_step = {
            "label": "Run predefined what-if scenarios",
            "detail": "Completed mechanistic runs are needed before trajectory and utility comparisons are meaningful.",
            "href": research_twin_url + "#what-if-scenarios",
        }
    elif research_provenance_count == 0:
        research_next_step = {
            "label": "Regenerate provenance metadata",
            "detail": "Traceability is incomplete until simulation metadata records exist for the current twin state.",
            "href": research_twin_url + "#provenance",
        }
    else:
        research_next_step = {
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
        "utility_explanation": "A cross-scenario ordering is undefined until a research utility function is specified.",
        "plain_language": "The platform can compare simulated strategies, but cannot prove what would have happened clinically.",
        "next_action": research_next_step,
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
        "latest_simulation_diagnostics": latest_simulation_diagnostics,
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
        "research_next_step": research_next_step,
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
    """Fail closed until individual prognosis has governed validation."""
    patient = get_object_or_404(models.Patient, pk=patient_id)

    context = {
        "patient": patient,
        "page_purpose": (
            "This page records the E1 scientific gate for individual prognosis "
            "and links back to observed longitudinal data."
        ),
        "interpretation_status": {
            "title": "PATIENT_SPECIFIC_PREDICTION_NOT_VALIDATED",
            "summary": (
                "No patient-specific PFS, OS, survival probability, or benefit "
                "estimate is emitted at the current evidence level."
            ),
            "limit": (
                "Population references and heuristic modifiers have not passed "
                "the governed validation required for individual prediction."
            ),
        },
        "confidence_panel": {
            "title": "Scientific validation gate",
            "score_label": "Validated individual-prediction confidence",
            "percentage": 0,
            "explanation": (
                "Zero represents unavailable governed validation, not a "
                "probability about this patient."
            ),
            "warning": (
                "INSUFFICIENT_VALIDATION: individual prognosis is suppressed "
                "and cannot be recovered by adding a disclaimer."
            ),
            "alert_class": "alert-danger",
        },
        "data_source_rows": [
            {
                "variable": "Observed assessment timeline",
                "value": f"{patient.assessments.count()} dated assessment(s)",
                "classification": "OBSERVED",
                "role": "Available for descriptive longitudinal review only.",
            },
            {
                "variable": "Recorded therapy timeline",
                "value": f"{patient.therapies.count()} dated therapy record(s)",
                "classification": "OBSERVED",
                "role": "Available for chronology; no counterfactual effect is inferred.",
            },
            {
                "variable": "Individual prognosis model",
                "value": "Not validated",
                "classification": "UNKNOWN",
                "role": "Suppressed by claims-policy-v1.",
            },
        ],
        "estimate_cards": [
            {
                "label": "Scientific status",
                "value": "INSUFFICIENT_VALIDATION",
                "detail": "No individual numerical prognosis is available.",
            }
        ],
        "timeline_rows": [],
        "adjustment_rows": [],
        "reference_info": {
            "title": "Validation requirement",
            "label": "Required before numerical individual prognosis",
            "value": "Separately governed source verification and external validation",
            "detail": (
                "The dormant heuristic research module is not evidence of "
                "patient-level applicability."
            ),
        },
        "what_can_be_concluded": [
            "Dated observed assessments and therapies can be reviewed.",
            "Missing validation is an explicit scientific result.",
        ],
        "what_cannot_be_concluded": [
            "Individual PFS, OS, or survival probability.",
            "Patient benefit under an observed or alternative therapy.",
            "A comparative or causal treatment effect.",
        ],
        "next_actions": {
            "primary": {
                "label": "Back to observed patient timeline",
                "href": reverse("clinic:patient_detail", args=[patient.id]),
            },
            "secondary": [
                {
                    "label": "Open Simple Research View",
                    "href": reverse("twin_engine:simple_research_view", args=[patient.id]),
                },
                {
                    "label": "Open Scientific Cockpit",
                    "href": reverse("twin_engine:research_cockpit", args=[patient.id]),
                },
            ],
        },
    }
    return render(request, "clinic/prognosis_timeline.html", context)


@login_required
@require_http_methods(["GET"])
def prognosis_api(request: HttpRequest, patient_id: int) -> JsonResponse:
    """Fail-closed JSON status for the unavailable individual prognosis."""
    patient = get_object_or_404(models.Patient, pk=patient_id)
    if not can_edit_patient(request.user, patient):
        return JsonResponse({"error": "forbidden"}, status=403)

    return JsonResponse({
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.UNKNOWN,
            output_kind="individual_prognosis_gate",
        ),
        "patient_id": patient_id,
        "status": "PATIENT_SPECIFIC_PREDICTION_NOT_VALIDATED",
        "estimate": None,
        "limitations": [
            "No patient-specific PFS, OS, survival probability, or benefit estimate is emitted.",
            "Population heuristics have not passed governed external validation for individual use.",
        ],
    })



# ══════════════════════════════════════════════════════════════════════════════
# EXPLORATORY REGIMEN-CONTEXT VIEWS
# ══════════════════════════════════════════════════════════════════════════════

@login_required
def regimen_suggestions(request: HttpRequest, patient_id: int) -> HttpResponse:
    """Display exploratory regimen context for a patient."""
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
    prior_agents = list(dict.fromkeys(agent for agent in prior_agents if agent))
    
    # Get neuropathy grade from symptoms
    neuropathy = 0
    if latest_symptoms:
        neuropathy = latest_symptoms.max_neuropathy_grade
    
    # Count therapy lines
    therapy_count = patient.therapies.count()
    line_of_therapy = max(therapy_count, 1)
    ecog = latest_symptoms.ecog_status if latest_symptoms else 0
    disease_setting = (
        "Frontline / newly diagnosed context"
        if line_of_therapy == 1
        else "Relapsed / previously treated context"
    )
    creatinine_value = None
    if latest_assessment and latest_assessment.creatinine_mg_dl is not None:
        creatinine_value = float(latest_assessment.creatinine_mg_dl)
    
    # Get suggestions
    suggestions = suggest_regimens(
        age=age,
        transplant_eligible=transplant_eligible,
        ecog=ecog,
        r_iss=latest_assessment.r_iss if latest_assessment else None,
        high_risk_cytogenetics=has_high_risk,
        line_of_therapy=line_of_therapy,
        prior_therapies=prior_agents,
        neuropathy_grade=neuropathy,
    )

    missing_inputs: list[str] = []
    if latest_assessment is None:
        missing_inputs.append("structured disease assessment")
    if not patient_cytos:
        missing_inputs.append("cytogenetic detail")
    if not prior_agents:
        missing_inputs.append("prior therapy components")
    if latest_symptoms is None:
        missing_inputs.append("symptom / toxicity assessment")
    if creatinine_value is None and not (latest_symptoms and latest_symptoms.crab_renal):
        missing_inputs.append("renal constraint summary")

    patient_detail_url = reverse("clinic:patient_detail", args=[patient.id])
    research_simple_url = reverse("twin_engine:simple_research_view", args=[patient.id])
    research_cockpit_url = reverse("twin_engine:research_cockpit", args=[patient.id])
    simulation_url = reverse("simulator:scenario_list")
    if latest_assessment is not None:
        simulation_url = f"{simulation_url}?twin_assessment_id={latest_assessment.pk}"

    context = {
        "patient": patient,
        "page_purpose": "Use this page to inspect how the current rule set groups regimen hypotheses and constraints from the available structured inputs.",
        "interpretation_status": {
            "title": "Exploratory regimen context",
            "summary": "This page organizes regimen-related context using heuristic rules whose source-level evidence requires verification. It is not a clinical treatment recommendation.",
            "limit": "This page cannot support a patient-specific comparative clinical conclusion.",
        },
        "data_source_rows": _build_regimen_data_source_rows(
            disease_setting=disease_setting,
            therapy_count=therapy_count,
            line_of_therapy=line_of_therapy,
            prior_agents=prior_agents,
            transplant_eligible=transplant_eligible,
            patient_cytos=patient_cytos,
            has_high_risk=has_high_risk,
            latest_assessment=latest_assessment,
            latest_symptoms=latest_symptoms,
            neuropathy=neuropathy,
            creatinine_value=creatinine_value,
        ),
        "regimen_sections": _build_regimen_sections(suggestions, missing_inputs=missing_inputs),
        "constraint_cards": _build_regimen_constraint_cards(suggestions.get("avoid", []), missing_inputs=missing_inputs),
        "what_can_be_concluded": [
            "This page can organize regimen-related hypotheses.",
            "This page can expose constraints and prior-therapy context.",
            "This page can help decide what to inspect next in the Simple Research View or Scientific Cockpit.",
        ],
        "what_cannot_be_concluded": [
            "It cannot establish clinical superiority.",
            "It cannot infer patient-specific benefit.",
            "It cannot determine what would have happened under another therapy.",
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
                {
                    "label": "View Algorithm Transparency",
                    "href": reverse("simulator:algorithm_transparency"),
                },
            ],
        },
    }
    return render(request, "clinic/regimen_suggestions.html", context)


def _build_regimen_data_source_rows(
    *,
    disease_setting: str,
    therapy_count: int,
    line_of_therapy: int,
    prior_agents: list[str],
    transplant_eligible: bool,
    patient_cytos: list[str],
    has_high_risk: bool,
    latest_assessment: models.Assessment | None,
    latest_symptoms: SymptomAssessment | None,
    neuropathy: int,
    creatinine_value: float | None,
) -> list[dict[str, str]]:
    cytogenetic_value = "No structured cytogenetic record surfaced on this page"
    cytogenetic_classification = "UNKNOWN"
    if patient_cytos:
        cytogenetic_value = ", ".join(patient_cytos)
        if has_high_risk:
            cytogenetic_value = f"{cytogenetic_value} (high-risk pattern present)"
        cytogenetic_classification = "RAW STRUCTURED"

    renal_parts: list[str] = []
    renal_classification = "UNKNOWN"
    if creatinine_value is not None:
        renal_parts.append(f"Creatinine {creatinine_value:.1f} mg/dL")
        renal_classification = "RAW STRUCTURED"
    if latest_symptoms and latest_symptoms.crab_renal:
        renal_parts.append("Structured CRAB renal flag present")
        renal_classification = "RAW STRUCTURED"
    if renal_parts:
        renal_parts.append("No hepatic variable is surfaced on this page")
        renal_value = "; ".join(renal_parts)
    else:
        renal_value = "Renal/hepatic constraints are not captured in the current page inputs"

    therapy_line_classification = "DERIVED" if therapy_count else "HEURISTIC"
    therapy_line_value = (
        f"Line {line_of_therapy}"
        if therapy_count
        else "Line 1 (default because no structured therapy history was recorded)"
    )

    prior_therapy_value = ", ".join(prior_agents) if prior_agents else "No structured prior therapy components recorded"
    prior_therapy_classification = "RAW STRUCTURED" if prior_agents else "UNKNOWN"

    neuropathy_value = "No structured symptom / toxicity assessment surfaced on this page"
    neuropathy_classification = "UNKNOWN"
    if latest_symptoms:
        neuropathy_value = f"Neuropathy grade {neuropathy}; ECOG {latest_symptoms.ecog_status}"
        neuropathy_classification = "RAW STRUCTURED"

    r_iss_value = latest_assessment.r_iss if latest_assessment and latest_assessment.r_iss else "Not recorded"
    r_iss_classification = "RAW STRUCTURED" if latest_assessment and latest_assessment.r_iss else "UNKNOWN"

    return [
        {
            "factor": "Disease setting",
            "value": disease_setting,
            "classification": therapy_line_classification,
            "role": "Derived from therapy line and used to route the page toward frontline or relapsed regimen buckets.",
        },
        {
            "factor": "Therapy line",
            "value": therapy_line_value,
            "classification": therapy_line_classification,
            "role": "Current line used by the rule set to switch between frontline, relapse, and later-line contexts.",
        },
        {
            "factor": "Prior therapies",
            "value": prior_therapy_value,
            "classification": prior_therapy_classification,
            "role": "Structured regimen components collapsed into agent exposure history for heuristic grouping.",
        },
        {
            "factor": "Transplant status",
            "value": "Eligible (age-based heuristic)" if transplant_eligible else "Ineligible (age-based heuristic)",
            "classification": "HEURISTIC",
            "role": "Age-based heuristic used by this page to separate transplant-oriented and non-transplant regimen contexts.",
        },
        {
            "factor": "Disease stage (R-ISS)",
            "value": r_iss_value,
            "classification": r_iss_classification,
            "role": "Structured staging context available to the rule set, but not the sole driver of regimen grouping.",
        },
        {
            "factor": "High-risk cytogenetics",
            "value": cytogenetic_value,
            "classification": cytogenetic_classification,
            "role": "Structured cytogenetic findings can shift the rule set toward more intensive or higher-risk exploratory contexts.",
        },
        {
            "factor": "Neuropathy / toxicity constraints",
            "value": neuropathy_value,
            "classification": neuropathy_classification,
            "role": "Structured symptom inputs are used to flag or deprioritize regimens that may worsen toxicity.",
        },
        {
            "factor": "Renal / hepatic constraints",
            "value": renal_value,
            "classification": renal_classification,
            "role": "Renal signals can affect lenalidomide or bortezomib handling; hepatic context is not surfaced in this page's current inputs.",
        },
        {
            "factor": "Regimen grouping logic",
            "value": "Bucketed by line, prior exposure, transplant heuristic, toxicity flags, and later-line investigational rules.",
            "classification": "HEURISTIC",
            "role": "Rule-based grouping determines whether a regimen appears as higher-priority, alternative, or constraint-flagged context.",
        },
    ]


def _build_regimen_sections(
    suggestions: dict[str, object], *, missing_inputs: list[str]
) -> list[dict[str, object]]:
    section_specs = [
        (
            "preferred",
            "Higher-priority exploratory context",
            "Rule-based bucket for regimens surfaced most prominently from the currently available inputs.",
        ),
        (
            "alternative",
            "Alternative exploratory context",
            "Rule-based bucket for other regimen contexts that remain plausible under the currently available inputs.",
        ),
        (
            "consider_in_clinical_trial",
            "Alternative exploratory context: later-line / investigational",
            "Later-line or investigational contexts surfaced by the current rules when therapy history is more advanced.",
        ),
    ]

    sections: list[dict[str, object]] = []
    for key, title, description in section_specs:
        regimens = suggestions.get(key) or []
        if not regimens:
            continue
        sections.append(
            {
                "title": title,
                "description": description,
                "cards": _build_regimen_cards(regimens, missing_inputs=missing_inputs),
            }
        )
    return sections


def _build_regimen_cards(
    regimens: list[dict[str, object]], *, missing_inputs: list[str]
) -> list[dict[str, object]]:
    cards: list[dict[str, object]] = []
    for regimen in regimens:
        trial_list = list(regimen.get("key_trials") or [])
        classification_labels = ["HEURISTIC"]
        classification_labels.append("LITERATURE-INFORMED" if trial_list else "NEEDS_EVIDENCE")
        cards.append(
            {
                "name": regimen.get("name", "Unnamed regimen"),
                "full_name": regimen.get("full_name", ""),
                "components": regimen.get("components", []),
                "why_it_appears": (
                    "This catalog entry was surfaced by the current heuristic "
                    "bucket rules for the provided inputs. The rule does not "
                    "estimate patient benefit."
                ),
                "classification_labels": classification_labels,
                "logic_basis": (
                    "Heuristic grouping with named source candidates."
                    if trial_list
                    else "Heuristic grouping with limited literature metadata surfaced on this page."
                ),
                "uncertainty_caveat": _build_regimen_uncertainty_caveat(missing_inputs),
                "population_signal": regimen.get("expected_response_rate") or "Population response signal not surfaced",
                "evidence_level": regimen.get("evidence_level") or "Not stated",
                "key_trials": trial_list,
                "considerations": [],
                "contraindications": [],
            }
        )
    return cards


def _build_regimen_constraint_cards(
    avoid_items: list[dict[str, object]], *, missing_inputs: list[str]
) -> list[dict[str, str]]:
    cards: list[dict[str, str]] = []
    for item in avoid_items:
        cards.append(
            {
                "agent": str(item.get("agent") or "Constraint flag"),
                "reason": (
                    "The provided inputs activated a configured research "
                    "constraint for this catalog entry."
                ),
                "classification": "HEURISTIC",
                "logic_basis": "Constraint flag from the current toxicity or comorbidity rules.",
                "uncertainty_caveat": _build_regimen_uncertainty_caveat(missing_inputs),
            }
        )
    return cards


def _build_regimen_uncertainty_caveat(missing_inputs: list[str]) -> str:
    base = "This card reflects heuristic grouping; its source-level evidence requires verification and it does not infer patient-specific benefit."
    if not missing_inputs:
        return base
    return f"{base} Missing structured inputs for this page: {', '.join(missing_inputs)}."


@login_required
@require_http_methods(["GET"])
def regimen_suggestions_api(request: HttpRequest, patient_id: int) -> JsonResponse:
    """JSON API for non-prescriptive regimen catalog context."""
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
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.HEURISTIC,
            output_kind="regimen_context",
        ),
        "patient_id": patient_id,
        "exploratory_reference_sets": _build_regimen_sections(
            suggestions,
            missing_inputs=[],
        ),
        "constraint_flags": _build_regimen_constraint_cards(
            suggestions.get("avoid", []),
            missing_inputs=[],
        ),
        "limitations": [
            "Catalog buckets are heuristic and their source-level evidence requires verification.",
            "No patient-specific benefit, safety, or comparative clinical conclusion is identified.",
        ],
    })

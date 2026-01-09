from __future__ import annotations

import json
import logging

from django.contrib.auth.decorators import login_required, user_passes_test
from django.db.models import Count, Prefetch, F
from django.http import HttpRequest, HttpResponse, HttpResponseForbidden
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

import django_filters

from . import forms, models

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
    chart_points = [
        {
            "date": assessment.date.isoformat(),
            "m": float(assessment.m_protein_g_dl) if assessment.m_protein_g_dl is not None else None,
            "flc": float(assessment.flc_ratio) if assessment.flc_ratio is not None else None,
            "ldh": float(assessment.ldH_u_l) if assessment.ldH_u_l is not None else None,
            "beta2m": float(assessment.beta2m_mg_l) if assessment.beta2m_mg_l is not None else None,
            "riss": (1 if assessment.r_iss == "I" else 2 if assessment.r_iss == "II" else 3 if assessment.r_iss == "III" else None),
        }
        for assessment in reversed(assessments)
    ]

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
    therapy_effects: dict[int, dict[str, object]] = {}
    try:
        from simulator.twin import build_patient_twin  # local import to avoid hard dependency at import time
    except Exception:  # pragma: no cover
        build_patient_twin = None

    for t in patient.therapies.all():
        baseline = (
            models.Assessment.objects.filter(patient=patient, date__lte=t.start_date)
            .order_by("-date")
            .first()
        )
        follow_until = t.end_date or None
        follow_qs = models.Assessment.objects.filter(patient=patient)
        if follow_until:
            follow_qs = follow_qs.filter(date__lte=follow_until)
        follow = follow_qs.order_by("-date").first()

        def maybe_float(v):
            try:
                return float(v) if v is not None else None
            except Exception:
                return None

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
    }
    return render(request, "clinic/patient_detail.html", context)


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

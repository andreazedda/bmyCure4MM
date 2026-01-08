from __future__ import annotations

import json

from django.contrib.auth.decorators import login_required, user_passes_test
from django.db.models import Count, Prefetch, F
from django.http import HttpRequest, HttpResponse, HttpResponseForbidden
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

import django_filters

from . import forms, models

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

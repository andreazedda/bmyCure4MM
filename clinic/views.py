from __future__ import annotations

from django.contrib.auth.decorators import login_required, user_passes_test
from django.db.models import Count, Prefetch, F
from django.http import HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

import django_filters

from . import forms, models

is_staff = user_passes_test(lambda u: u.is_staff)


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
    chart_points = [
        {
            "date": assessment.date.isoformat(),
            "m": float(assessment.m_protein_g_dl) if assessment.m_protein_g_dl is not None else None,
            "flc": float(assessment.flc_ratio) if assessment.flc_ratio is not None else None,
        }
        for assessment in reversed(assessments)
    ]

    context = {
        "patient": patient,
        "assessments": assessments,
        "cytogenetics": patient.cytogenetics.all(),
        "therapies": patient.therapies.all(),
        "assessment_form": forms.AssessmentForm(),
        "therapy_form": forms.PatientTherapyForm(),
        "chart_points": chart_points,
    }
    return render(request, "clinic/patient_detail.html", context)


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
@is_staff
def assessment_new(request: HttpRequest, patient_id: int) -> HttpResponse:
    patient = get_object_or_404(models.Patient, pk=patient_id)
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

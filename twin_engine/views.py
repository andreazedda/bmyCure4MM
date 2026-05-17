from __future__ import annotations

import json
from datetime import timedelta
from pathlib import Path

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied, ValidationError
from django.http import HttpRequest, HttpResponse, HttpResponseBadRequest
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.utils import timezone

from clinic.models import Assessment, Patient

from .calibration import calibrate_patient_parameters
from .causal import build_causal_assumption_set, classify_estimand
from .cockpit import (
    build_assessment_recommendations,
    build_research_cockpit_context,
    build_research_glossary,
    run_developer_checks,
    summarize_checks,
    write_local_feedback,
)
from .counterfactual import run_counterfactual
from .models import CausalAssumptionSet, CounterfactualRun, ObservationResidual, SimulationRunMetadata
from .provenance import CURRENT_MODEL_VERSION, hash_json, record_simulation_metadata
from .report_builder import write_json_artifact
from .simulation_bridge import run_patient_simulation
from .simple_view import build_simple_patient_story
from .state_model import get_current_twin_state, initialize_from_assessment, serialize_state
from .therapy_schedule import build_therapy_schedule
from .toxicity_model import compute_toxicity_constraints
from .validators import validate_assessment_minimum_fields, validate_patient_access


@login_required
def patient_twin_detail(request: HttpRequest, patient_id: int) -> HttpResponse:
    return simple_research_view(request, patient_id)


@login_required
def simple_research_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    patient = _get_research_patient(request, patient_id)
    story = build_simple_patient_story(patient, include_developer_links=request.user.is_staff)
    context = {
        "patient": patient,
        "story": story,
    }
    return render(request, "twin_engine/simple_research_view.html", context)


@login_required
def research_cockpit_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    patient = _get_research_patient(request, patient_id)
    context = build_research_cockpit_context(patient, include_developer_checks=request.user.is_staff)
    return render(request, "twin_engine/research_cockpit.html", context)


@login_required
def research_glossary_view(request: HttpRequest) -> HttpResponse:
    return render(request, "twin_engine/research_glossary.html", {"terms": build_research_glossary()})


@login_required
def initialize_twin_from_assessment_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    patient = _get_research_patient(request, patient_id)

    if request.method != "POST":
        context = build_research_cockpit_context(patient, include_developer_checks=request.user.is_staff)
        context["initialize_guidance"] = build_assessment_recommendations(patient)
        return render(request, "twin_engine/initialize_twin.html", context)

    assessment_id = request.POST.get("assessment_id")
    if assessment_id:
        assessment = get_object_or_404(Assessment, pk=assessment_id, patient=patient)
    else:
        guidance = build_assessment_recommendations(patient)
        recommended = guidance.get("recommended")
        assessment = recommended["assessment"] if recommended else patient.assessments.order_by("-date").first()
        if assessment is None:
            messages.error(request, "No assessments available to initialize the research twin.")
            return redirect(reverse("twin_engine:research_cockpit", args=[patient.id]))

    try:
        initialize_from_assessment(assessment, user=request.user)
    except ValidationError as exc:
        messages.error(request, "; ".join(exc.messages))
        return redirect(reverse("twin_engine:initialize_twin_from_assessment", args=[patient.id]))
    messages.success(request, "Research twin initialized from the selected assessment.")
    return redirect(reverse("twin_engine:research_cockpit", args=[patient.id]))


@login_required
def calibrate_twin_from_history_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    if request.method != "POST":
        return HttpResponseBadRequest("POST required")
    patient = _get_research_patient(request, patient_id)

    current_state = get_current_twin_state(patient)
    assessments = list(patient.assessments.order_by("date"))
    if not assessments:
        messages.error(request, "Calibration requires at least one assessment.")
        return redirect(reverse("twin_engine:research_cockpit", args=[patient.id]))
    if current_state is None:
        current_state = initialize_from_assessment(assessments[0], user=request.user)

    start_date = request.POST.get("start_date") or None
    end_date = request.POST.get("end_date") or None
    if start_date:
        assessments = [item for item in assessments if item.date.isoformat() >= start_date]
    if end_date:
        assessments = [item for item in assessments if item.date.isoformat() <= end_date]

    try:
        result = calibrate_patient_parameters(
            patient,
            current_state,
            assessments,
            patient.therapies.all(),
        )
    except ValidationError as exc:
        messages.error(request, "; ".join(exc.messages))
    else:
        messages.success(
            request,
            f"Calibration completed with objective {result['optimizer']['objective']:.4f}.",
        )
    return redirect(reverse("twin_engine:research_cockpit", args=[patient.id]))


@login_required
def run_patient_simulation_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    if request.method != "POST":
        return HttpResponseBadRequest("POST required")
    patient = _get_research_patient(request, patient_id)
    base_state = _get_base_state(patient, request.POST.get("base_state_id"))
    horizon_days = int(request.POST.get("horizon_days") or 90)
    end_date = base_state.state_date + timedelta(days=horizon_days)
    therapy_schedule = build_therapy_schedule(patient, base_state.state_date, end_date)
    simulation_result = run_patient_simulation(base_state, therapy_schedule=therapy_schedule, horizon_days=horizon_days)
    toxicity_constraints = compute_toxicity_constraints(patient)

    run_token = int(timezone.now().timestamp())
    patient_reference = hash_json({"patient_id": patient.id})[:12]
    report_payload = {
        "label": "research simulation",
        "counterfactual_type": "mechanistic model counterfactual",
        "causal_interpretation": "unvalidated exploratory hypothesis",
        "patient_reference": patient_reference,
        "base_state": serialize_state(base_state),
        "summary": simulation_result["summary"],
        "predicted_biomarkers": simulation_result["summary"].get("predicted_biomarkers"),
        "toxicity_constraints": toxicity_constraints,
        "solver_inputs": simulation_result["solver_inputs"]["raw_parameters"],
        "missing_doses": simulation_result["solver_inputs"]["missing_doses"],
        "generated_at": timezone.now().isoformat(),
        "user_id": request.user.id,
    }
    artifact_url, _ = write_json_artifact(
        "patient_simulation",
        report_payload,
        patient_id=patient.id,
        run_id=run_token,
        folder_name="research_reports",
    )
    metadata = record_simulation_metadata(
        twin_state=base_state,
        model_version=CURRENT_MODEL_VERSION,
        solver_name="MathematicalModel",
        input_payload={
            "patient_id": patient.id,
            "base_state_id": base_state.id,
            "horizon_days": horizon_days,
            "therapy_schedule": therapy_schedule,
        },
        solver_parameters=simulation_result["solver_inputs"]["raw_parameters"],
        output_payload={
            "artifact_url": artifact_url,
            "summary": simulation_result["summary"],
        },
    )

    context = _build_research_context(patient)
    context.update(
        {
            "simulation_result": simulation_result,
            "simulation_artifact_url": artifact_url,
            "simulation_metadata": metadata,
            "toxicity_constraints": toxicity_constraints,
        }
    )
    return render(request, "twin_engine/patient_simulation_result.html", context)


@login_required
def run_counterfactual_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    if request.method != "POST":
        return HttpResponseBadRequest("POST required")
    patient = _get_research_patient(request, patient_id)
    base_state = _get_base_state(patient, request.POST.get("base_state_id"))
    horizon_days = int(request.POST.get("horizon_days") or 90)

    intervention_json = request.POST.get("intervention_json") or "{}"
    try:
        intervention_definition = json.loads(intervention_json)
    except json.JSONDecodeError:
        messages.error(request, "Intervention JSON is invalid.")
        return redirect(reverse("twin_engine:research_cockpit", args=[patient.id]))

    try:
        counterfactual_run = run_counterfactual(
            patient,
            base_state,
            intervention_definition,
            horizon_days,
            user=request.user,
        )
    except Exception as exc:
        messages.error(request, f"Counterfactual run failed: {exc}")
        return redirect(reverse("twin_engine:research_cockpit", args=[patient.id]))

    messages.success(request, "Counterfactual research run completed.")
    return redirect(reverse("twin_engine:counterfactual_report", args=[patient.id, counterfactual_run.id]))


@login_required
def counterfactual_report_view(request: HttpRequest, patient_id: int, run_id: int) -> HttpResponse:
    patient = _get_research_patient(request, patient_id)
    counterfactual_run = get_object_or_404(
        CounterfactualRun.objects.select_related("base_twin_state", "alternative_regimen"),
        pk=run_id,
        patient=patient,
    )
    report_payload = _load_json_artifact(counterfactual_run.report_artifact) or counterfactual_run.simulation_summary
    context = {
        "patient": patient,
        "counterfactual_run": counterfactual_run,
        "report_payload": report_payload,
        "toxicity_constraints": (report_payload or {}).get("toxicity_constraints") or (counterfactual_run.simulation_summary or {}).get("toxicity_constraints") or {},
    }
    return render(request, "twin_engine/counterfactual_report.html", context)


@login_required
def causal_assumption_set_view(request: HttpRequest, patient_id: int) -> HttpResponse:
    patient = _get_research_patient(request, patient_id)

    if request.method == "POST":
        try:
            graph_definition = json.loads(request.POST.get("graph_definition") or "{}")
            intervention = json.loads(request.POST.get("intervention") or "{}")
            outcome = json.loads(request.POST.get("outcome") or "{}")
            assumptions = json.loads(request.POST.get("assumptions") or "{}")
            adjustment_set = json.loads(request.POST.get("adjustment_set") or "[]")
            identification_status = request.POST.get("identification_status") or CausalAssumptionSet.IDENT_MECHANISTIC_ONLY
            build_causal_assumption_set(
                patient=patient,
                created_by=request.user,
                name=request.POST.get("name") or f"Patient {patient.id} research assumptions",
                graph_definition=graph_definition,
                variables=list(graph_definition.get("nodes", [])),
                intervention=intervention,
                outcome=outcome,
                adjustment_set=adjustment_set,
                assumptions=assumptions,
                identification_status=identification_status,
                notes=request.POST.get("notes", ""),
            )
        except (ValidationError, json.JSONDecodeError) as exc:
            messages.error(request, f"Causal assumption set could not be saved: {exc}")
        else:
            messages.success(request, "Causal assumption set saved.")
        return redirect(reverse("twin_engine:causal_assumption_set", args=[patient.id]))

    sets = patient.causal_assumption_sets.order_by("-created_at")
    assumption_rows = [
        {
            "assumption_set": assumption_set,
            "classification": classify_estimand(
                graph_definition=assumption_set.graph_definition,
                intervention=assumption_set.intervention,
                outcome=assumption_set.outcome,
                adjustment_set=assumption_set.adjustment_set,
                identification_status=assumption_set.identification_status,
            ),
        }
        for assumption_set in sets
    ]
    context = {
        "patient": patient,
        "assumption_rows": assumption_rows,
    }
    return render(request, "twin_engine/causal_assumptions.html", context)


@login_required
def developer_console_view(request: HttpRequest) -> HttpResponse:
    if not request.user.is_staff:
        raise PermissionDenied("Developer console requires staff access.")
    if request.method == "POST":
        feedback_payload = {
            "patient_id": request.POST.get("patient_id") or "",
            "category": request.POST.get("category") or "general",
            "message": request.POST.get("message") or "",
        }
        if feedback_payload["message"].strip():
            write_local_feedback(request.user, feedback_payload)
            messages.success(request, "Local-only developer feedback captured.")
        else:
            messages.warning(request, "Feedback message was empty.")
        return redirect(reverse("twin_engine:developer_console"))
    check_groups = run_developer_checks()
    focus_patient = None
    focus_patient_id = request.GET.get("patient_id") or ""
    if focus_patient_id.isdigit():
        focus_patient = Patient.objects.filter(pk=int(focus_patient_id)).first()

    navigation_links = {
        "patient_page_url": reverse("clinic:patient_list"),
        "simple_view_url": "",
        "cockpit_url": "",
        "developer_console_url": reverse("twin_engine:developer_console"),
        "glossary_url": reverse("twin_engine:research_glossary"),
    }
    if focus_patient is not None:
        navigation_links.update(
            {
                "patient_page_url": reverse("clinic:patient_detail", args=[focus_patient.id]),
                "simple_view_url": reverse("twin_engine:simple_research_view", args=[focus_patient.id]),
                "cockpit_url": reverse("twin_engine:research_cockpit", args=[focus_patient.id]),
                "developer_console_url": reverse("twin_engine:developer_console") + f"?patient_id={focus_patient.id}",
            }
        )
    context = {
        "check_summary": summarize_checks(check_groups),
        "patients_with_twins": Patient.objects.filter(twin_states__isnull=False).distinct().order_by("id")[:50],
        "focus_patient": focus_patient,
        "navigation_links": navigation_links,
    }
    return render(request, "twin_engine/developer_console.html", context)


def _get_research_patient(request: HttpRequest, patient_id: int) -> Patient:
    patient = get_object_or_404(Patient.objects.prefetch_related("assessments", "therapies"), pk=patient_id)
    if not validate_patient_access(request.user, patient):
        raise PermissionDenied("Not allowed")
    return patient


def _get_base_state(patient: Patient, base_state_id: str | None):
    if base_state_id:
        return get_object_or_404(patient.twin_states, pk=base_state_id)
    state = get_current_twin_state(patient)
    if state is None:
        guidance = build_assessment_recommendations(patient)
        recommended = guidance.get("recommended")
        if recommended is None:
            raise ValidationError("No assessment available to initialize a research twin state.")
        state = initialize_from_assessment(recommended["assessment"])
    return state


def _build_research_context(patient: Patient) -> dict[str, object]:
    latest_assessment = patient.assessments.order_by("-date").first()
    current_state = get_current_twin_state(patient)
    latest_residual = None
    last_calibration_status = "Not initialized"
    provenance_records = SimulationRunMetadata.objects.none()
    if current_state is not None:
        latest_residual = ObservationResidual.objects.filter(twin_state=current_state).select_related("assessment").order_by("-created_at").first()
        last_calibration_status = current_state.get_method_display()
        provenance_records = SimulationRunMetadata.objects.filter(twin_state=current_state).order_by("-created_at")[:10]

    if latest_assessment is not None:
        try:
            field_status = validate_assessment_minimum_fields(latest_assessment)
            missing_required_data = field_status["missing"]
        except ValidationError as exc:
            missing_required_data = list(exc.messages)
    else:
        missing_required_data = ["assessment"]

    schedule_validation = None
    if current_state is not None:
        readiness_schedule = build_therapy_schedule(
            patient,
            current_state.state_date,
            current_state.state_date + timedelta(days=90),
        )
        schedule_validation = readiness_schedule.get("validation")
    toxicity_constraints = compute_toxicity_constraints(patient)

    return {
        "patient": patient,
        "current_twin_state": current_state,
        "latest_residual": latest_residual,
        "last_calibration_status": last_calibration_status,
        "missing_required_data": missing_required_data,
        "schedule_validation": schedule_validation,
        "counterfactual_runs": patient.counterfactual_runs.select_related("base_twin_state", "alternative_regimen").order_by("-created_at")[:10],
        "causal_assumption_sets": patient.causal_assumption_sets.order_by("-created_at")[:10],
        "provenance_records": provenance_records,
        "latest_assessment": latest_assessment,
        "toxicity_constraints": toxicity_constraints,
    }


def _load_json_artifact(artifact_url: str):
    if not artifact_url:
        return None
    media_url = settings.MEDIA_URL.rstrip("/")
    relative = artifact_url[len(media_url) + 1 :] if artifact_url.startswith(media_url + "/") else artifact_url.lstrip("/")
    path = Path(settings.MEDIA_ROOT) / relative
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))

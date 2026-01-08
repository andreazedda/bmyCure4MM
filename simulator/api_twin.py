from __future__ import annotations

from django.contrib.auth.decorators import login_required
from django.http import HttpRequest, JsonResponse
from django.urls import reverse
from django.views.decorators.http import require_GET

from clinic.models import Assessment
from simulator.twin import build_patient_twin
from simulator.permissions import is_editor

from .access import get_accessible_assessment_by_id


@login_required
@require_GET
def twin_preview(request: HttpRequest) -> JsonResponse:
    is_privileged = request.user.is_staff or is_editor(request.user)

    raw_id = request.GET.get("id")
    if not raw_id:
        return JsonResponse({"error": "Missing query parameter: id"}, status=400)

    try:
        assessment_id = int(raw_id)
    except (TypeError, ValueError):
        return JsonResponse({"error": "Invalid id"}, status=400)

    if is_privileged:
        assessment = Assessment.objects.select_related("patient").filter(pk=assessment_id).first()
    else:
        assessment = get_accessible_assessment_by_id(request.user, assessment_id)
    if not assessment:
        return JsonResponse({"error": "Assessment not found"}, status=404)

    payload = build_patient_twin(assessment)
    patient = getattr(assessment, "patient", None)
    patient_id = getattr(patient, "pk", None)
    response = {
        "assessment": {
            "id": assessment.pk,
            "date": str(getattr(assessment, "date", "")),
            "patient": {
                "id": patient_id,
                "mrn": getattr(patient, "mrn", None),
                "name": f"{getattr(patient, 'first_name', '')} {getattr(patient, 'last_name', '')}".strip() or None,
                "detail_url": reverse("clinic:patient_detail", args=[patient_id]) if patient_id else None,
            },
            "inputs": {
                "r_iss": getattr(assessment, "r_iss", None),
                "ldh_u_l": float(assessment.ldH_u_l) if getattr(assessment, "ldH_u_l", None) is not None else None,
                "beta2m_mg_l": float(assessment.beta2m_mg_l) if getattr(assessment, "beta2m_mg_l", None) is not None else None,
                "flc_ratio": float(assessment.flc_ratio) if getattr(assessment, "flc_ratio", None) is not None else None,
            },
        },
        "twin": payload,
    }
    return JsonResponse(response)

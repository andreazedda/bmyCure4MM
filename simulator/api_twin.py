from __future__ import annotations

from django.contrib.auth.decorators import login_required
from django.http import HttpRequest, JsonResponse
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
    response = {
        "assessment": {
            "id": assessment.pk,
            "date": str(getattr(assessment, "date", "")),
            "patient": {
                "mrn": getattr(getattr(assessment, "patient", None), "mrn", None),
            },
        },
        "twin": payload,
    }
    return JsonResponse(response)

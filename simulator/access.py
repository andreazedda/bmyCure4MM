from __future__ import annotations

from django.db.models import Q, QuerySet

from clinic.models import Assessment

from .permissions import is_editor


DEMO_MRN_PREFIX = "DEMO"


def accessible_assessments(user, *, base_qs: QuerySet | None = None) -> QuerySet:
    """Return assessments accessible to the given user.

    Policy:
    - staff/editor: all assessments
    - non-staff: assessments for owned patients OR demo patients (MRN startswith DEMO)
    """
    qs = base_qs if base_qs is not None else Assessment.objects.all()

    if not user or not getattr(user, "is_authenticated", False):
        return qs.none()

    if getattr(user, "is_staff", False) or is_editor(user):
        return qs

    return qs.filter(Q(patient__owner=user) | Q(patient__mrn__startswith=DEMO_MRN_PREFIX))


def get_accessible_assessment_by_id(user, assessment_id: int):
    qs = accessible_assessments(user, base_qs=Assessment.objects.select_related("patient"))
    return qs.filter(pk=assessment_id).first()

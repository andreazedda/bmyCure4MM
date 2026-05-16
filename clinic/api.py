from __future__ import annotations

from rest_framework import routers, serializers, viewsets
from django.db.models import Q

from . import models


class PatientSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Patient
        fields = "__all__"


class AssessmentSerializer(serializers.ModelSerializer):
    class Meta:
        model = models.Assessment
        fields = "__all__"


class PatientViewSet(viewsets.ModelViewSet):
    queryset = models.Patient.objects.none()
    serializer_class = PatientSerializer

    def get_queryset(self):
        user = self.request.user
        queryset = models.Patient.objects.all()
        if not getattr(user, "is_authenticated", False):
            return queryset.none()
        if getattr(user, "is_staff", False):
            return queryset
        return queryset.filter(Q(owner=user) | Q(owner__isnull=True, mrn__startswith="DEMO"))


class AssessmentViewSet(viewsets.ModelViewSet):
    queryset = models.Assessment.objects.none()
    serializer_class = AssessmentSerializer

    def get_queryset(self):
        user = self.request.user
        queryset = models.Assessment.objects.select_related("patient").all()
        if not getattr(user, "is_authenticated", False):
            return queryset.none()
        if getattr(user, "is_staff", False):
            return queryset
        return queryset.filter(Q(patient__owner=user) | Q(patient__owner__isnull=True, patient__mrn__startswith="DEMO"))


router = routers.DefaultRouter()
router.register(r"patients", PatientViewSet)
router.register(r"assessments", AssessmentViewSet)

urlpatterns = router.urls


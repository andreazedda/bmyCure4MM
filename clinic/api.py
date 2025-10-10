from __future__ import annotations

from rest_framework import routers, serializers, viewsets

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
    queryset = models.Patient.objects.all()
    serializer_class = PatientSerializer


class AssessmentViewSet(viewsets.ModelViewSet):
    queryset = models.Assessment.objects.select_related("patient").all()
    serializer_class = AssessmentSerializer


router = routers.DefaultRouter()
router.register(r"patients", PatientViewSet)
router.register(r"assessments", AssessmentViewSet)

urlpatterns = router.urls


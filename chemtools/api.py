from __future__ import annotations

from rest_framework import routers

from . import viewsets


router = routers.DefaultRouter()
router.register(r"jobs", viewsets.ChemJobViewSet, basename="chemjob")

urlpatterns = router.urls

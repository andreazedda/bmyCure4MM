from __future__ import annotations

from django.contrib import admin

from . import models


@admin.register(models.Scenario)
class ScenarioAdmin(admin.ModelAdmin):
    list_display = ("title", "clinical_stage", "risk_stratification", "active", "updated")
    list_filter = ("clinical_stage", "active")
    search_fields = ("title", "summary", "risk_stratification")
    filter_horizontal = ("recommended_regimens",)


@admin.register(models.SimulationAttempt)
class SimulationAttemptAdmin(admin.ModelAdmin):
    list_display = ("scenario", "user", "selected_regimen", "predicted_response", "confidence", "is_guideline_aligned", "submitted")
    list_filter = ("is_guideline_aligned", "predicted_response")
    search_fields = ("scenario__title", "user__username", "notes")
    raw_id_fields = ("scenario", "user", "selected_regimen")


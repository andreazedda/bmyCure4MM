from django.contrib import admin

from .models import (
    CausalAssumptionSet,
    CounterfactualRun,
    ObservationResidual,
    PatientTwinState,
    SimulationRunMetadata,
)


@admin.register(PatientTwinState)
class PatientTwinStateAdmin(admin.ModelAdmin):
    list_display = ("patient", "state_date", "is_current", "method", "model_version")
    list_filter = ("is_current", "method", "model_version")
    search_fields = ("patient__mrn", "patient__last_name")


@admin.register(ObservationResidual)
class ObservationResidualAdmin(admin.ModelAdmin):
    list_display = ("patient", "twin_state", "assessment", "rmse", "mae", "created_at")
    list_filter = ("created_at",)
    search_fields = ("patient__mrn", "patient__last_name")


@admin.register(CounterfactualRun)
class CounterfactualRunAdmin(admin.ModelAdmin):
    list_display = ("patient", "base_twin_state", "status", "alternative_regimen", "created_at")
    list_filter = ("status", "created_at")
    search_fields = ("patient__mrn", "patient__last_name")


@admin.register(CausalAssumptionSet)
class CausalAssumptionSetAdmin(admin.ModelAdmin):
    list_display = ("name", "patient", "identification_status", "created_at")
    list_filter = ("identification_status", "created_at")
    search_fields = ("name", "patient__mrn", "patient__last_name")


@admin.register(SimulationRunMetadata)
class SimulationRunMetadataAdmin(admin.ModelAdmin):
    list_display = ("model_version", "solver_name", "twin_state", "counterfactual_run", "created_at")
    list_filter = ("model_version", "solver_name", "created_at")

from django.contrib import admin

from .models import (
    CausalAssumptionSet,
    CounterfactualRun,
    ObservationResidual,
    PatientTwinState,
    ResearchRunInvalidation,
    ResearchRunManifest,
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
    list_display = ("model_version", "solver_name", "manifest", "twin_state", "counterfactual_run", "created_at")
    list_filter = ("model_version", "solver_name", "created_at")


class ImmutableResearchAdmin(admin.ModelAdmin):
    def has_add_permission(self, request):
        return False

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return False


@admin.register(ResearchRunManifest)
class ResearchRunManifestAdmin(ImmutableResearchAdmin):
    list_display = ("run_id", "model_id", "model_version", "dataset_version", "status", "created_at")
    list_filter = ("model_id", "model_version", "status", "intended_use_level", "epistemic_label")
    search_fields = ("run_id", "dataset_id", "dataset_sha256", "record_subset_sha256")


@admin.register(ResearchRunInvalidation)
class ResearchRunInvalidationAdmin(ImmutableResearchAdmin):
    list_display = ("manifest", "change_kind", "change_sha256", "created_at")
    list_filter = ("change_kind", "created_at")
    search_fields = ("manifest__run_id", "change_sha256", "reason")

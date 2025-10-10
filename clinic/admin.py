from __future__ import annotations

from django.contrib import admin

from . import models


class AssessmentInline(admin.TabularInline):
    model = models.Assessment
    extra = 0
    fields = ("date", "m_protein_g_dl", "flc_ratio", "response")


@admin.register(models.Patient)
class PatientAdmin(admin.ModelAdmin):
    list_display = ("mrn", "last_name", "first_name", "diagnosis_date")
    search_fields = ("mrn", "last_name", "first_name")
    list_filter = ("sex",)
    inlines = [AssessmentInline]


@admin.register(models.CytogeneticAbnormality)
class CytogeneticAbnormalityAdmin(admin.ModelAdmin):
    list_display = ("code", "description")
    search_fields = ("code", "description")


@admin.register(models.PatientCytogenetics)
class PatientCytogeneticsAdmin(admin.ModelAdmin):
    list_display = ("patient", "abnormality", "detected_on", "method")
    search_fields = ("patient__mrn", "patient__last_name", "abnormality__code")
    list_filter = ("method", "detected_on")


@admin.register(models.Assessment)
class AssessmentAdmin(admin.ModelAdmin):
    list_display = ("patient", "date", "m_protein_g_dl", "flc_ratio", "response", "r_iss")
    search_fields = ("patient__mrn", "patient__last_name")
    list_filter = ("response", "r_iss")
    date_hierarchy = "date"


@admin.register(models.Regimen)
class RegimenAdmin(admin.ModelAdmin):
    list_display = ("name", "line", "intent")
    search_fields = ("name", "components")


@admin.register(models.PatientTherapy)
class PatientTherapyAdmin(admin.ModelAdmin):
    list_display = ("patient", "regimen", "start_date", "end_date", "outcome")
    search_fields = ("patient__mrn", "patient__last_name", "regimen__name")
    list_filter = ("regimen__line",)

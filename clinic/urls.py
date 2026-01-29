from django.urls import path

from . import views

urlpatterns = [
    path("", views.dashboard, name="dashboard"),
    path("patients/", views.patient_list, name="patient_list"),
    path("regimens/", views.regimen_list, name="regimen_list"),
    path("regimens/new/", views.regimen_new, name="regimen_new"),
    path("patients/new/", views.patient_new, name="patient_new"),
    path("patients/<int:pk>/", views.patient_detail, name="patient_detail"),
    path("patients/<int:pk>/edit/", views.patient_edit, name="patient_edit"),
    path(
        "patients/<int:patient_id>/assessments/new/",
        views.assessment_new,
        name="assessment_new",
    ),
    # Symptom Assessments
    path(
        "patients/<int:patient_id>/symptoms/new/",
        views.symptom_assessment_new,
        name="symptom_assessment_new",
    ),
    path(
        "patients/<int:patient_id>/symptoms/",
        views.symptom_assessment_list,
        name="symptom_assessment_list",
    ),
    # Prognosis & Timeline
    path(
        "patients/<int:patient_id>/prognosis/",
        views.prognosis_timeline,
        name="prognosis_timeline",
    ),
    path(
        "api/patients/<int:patient_id>/prognosis/",
        views.prognosis_api,
        name="prognosis_api",
    ),
    # Regimen Suggestions
    path(
        "patients/<int:patient_id>/regimens/",
        views.regimen_suggestions,
        name="regimen_suggestions",
    ),
    path(
        "api/patients/<int:patient_id>/regimens/",
        views.regimen_suggestions_api,
        name="regimen_suggestions_api",
    ),
]


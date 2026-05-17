from django.urls import path

from . import views


app_name = "twin_engine"


urlpatterns = [
    path("developer/", views.developer_console_view, name="developer_console"),
    path("glossary/", views.research_glossary_view, name="research_glossary"),
    path("patient/<int:patient_id>/simple/", views.simple_research_view, name="simple_research_view"),
    path("patient/<int:patient_id>/cockpit/", views.research_cockpit_view, name="research_cockpit"),
    path("patients/<int:patient_id>/twin/", views.patient_twin_detail, name="patient_twin_detail"),
    path(
        "patients/<int:patient_id>/twin/initialize/",
        views.initialize_twin_from_assessment_view,
        name="initialize_twin_from_assessment",
    ),
    path(
        "patients/<int:patient_id>/twin/calibrate/",
        views.calibrate_twin_from_history_view,
        name="calibrate_twin_from_history",
    ),
    path(
        "patients/<int:patient_id>/simulation/",
        views.run_patient_simulation_view,
        name="run_patient_simulation",
    ),
    path(
        "patients/<int:patient_id>/counterfactual/run/",
        views.run_counterfactual_view,
        name="run_counterfactual",
    ),
    path(
        "patients/<int:patient_id>/counterfactual/<int:run_id>/",
        views.counterfactual_report_view,
        name="counterfactual_report",
    ),
    path(
        "patients/<int:patient_id>/causal/",
        views.causal_assumption_set_view,
        name="causal_assumption_set",
    ),
]

from django.urls import path

from . import api_twin, views, views_cohort, views_manage

app_name = "simulator"

urlpatterns = [
    path("", views.scenario_list, name="scenario_list"),
    path("getting-started/", views.getting_started, name="getting_started"),
    path("tutorial/", views.interactive_tutorial, name="interactive_tutorial"),
    path("diagnostics/visibility/", views.visibility_diagnostics, name="visibility_diagnostics"),
    path("manage/", views_manage.ScenarioManageListView.as_view(), name="scenario_manage_list"),
    path("scenarios/new/", views_manage.ScenarioCreateView.as_view(), name="scenario_create"),
    path("scenarios/<int:pk>/edit/", views_manage.ScenarioUpdateView.as_view(), name="scenario_edit"),
    path("scenarios/<int:pk>/duplicate/", views_manage.duplicate_scenario, name="scenario_duplicate"),
    path(
        "scenarios/<int:pk>/regimens/add/",
        views_manage.add_regimen_to_scenario,
        name="scenario_regimen_add",
    ),
    path(
        "scenarios/<int:pk>/regimens/remove/",
        views_manage.remove_regimen_from_scenario,
        name="scenario_regimen_remove",
    ),
    path(
        "scenarios/<int:pk>/simulate/",
        views_manage.simulate_scenario,
        name="scenario_simulate",
    ),
    path("api/twin/preview/", api_twin.twin_preview, name="twin_preview"),
    path("regimens/", views_manage.RegimenListView.as_view(), name="regimen_list"),
    path("regimens/new/", views_manage.RegimenCreateView.as_view(), name="regimen_create"),
    path("regimens/<int:pk>/edit/", views_manage.RegimenUpdateView.as_view(), name="regimen_edit"),
    path("regimens/<int:pk>/delete/", views_manage.RegimenDeleteView.as_view(), name="regimen_delete"),
    path("attempts/<int:pk>/export/", views_manage.export_attempt, name="attempt_export"),
    path("attempts/<int:pk>/plot/embed/", views_manage.attempt_plot_embed, name="attempt_plot_embed"),
    # Optimization experiments
    path("experiments/", views_manage.ExperimentsView.as_view(), name="experiments"),
    path("experiments/run/<int:pk>/", views_manage.run_experiment, name="run_experiment"),
    # Virtual cohort (in-silico trials)
    path("cohort/", views_cohort.CohortView.as_view(), name="cohort"),
    path("cohort/run/", views_cohort.cohort_run, name="cohort_run"),
    path("cohort/export/<str:cohort_id>/", views_cohort.cohort_export, name="cohort_export"),
    path("<int:pk>/", views.scenario_detail, name="scenario_detail"),
]

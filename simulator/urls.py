from django.urls import path

from . import views
from . import views_manage

urlpatterns = [
    path("", views.scenario_list, name="scenario_list"),
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
    path("regimens/", views_manage.RegimenListView.as_view(), name="regimen_list"),
    path("regimens/new/", views_manage.RegimenCreateView.as_view(), name="regimen_create"),
    path("regimens/<int:pk>/edit/", views_manage.RegimenUpdateView.as_view(), name="regimen_edit"),
    path("regimens/<int:pk>/delete/", views_manage.RegimenDeleteView.as_view(), name="regimen_delete"),
    path("<int:pk>/", views.scenario_detail, name="scenario_detail"),
]

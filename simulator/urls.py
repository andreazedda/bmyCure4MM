from django.urls import path

from . import views

urlpatterns = [
    path("", views.scenario_list, name="scenario_list"),
    path("<int:pk>/", views.scenario_detail, name="scenario_detail"),
]


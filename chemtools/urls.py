from django.urls import path

from . import views

urlpatterns = [
    path("", views.tools_home, name="tools_home"),
    path("params/", views.drug_params, name="drug_params"),
    path("binding/", views.binding_viz, name="binding_viz"),
    path("similarity/", views.similarity, name="similarity"),
    path("job/<int:pk>/", views.job_detail, name="job_detail"),
    path("job/<int:pk>/retry/", views.retry_job, name="job_retry"),
    path("job/<int:pk>/status.json", views.job_status, name="job_status"),
    path("job/<int:pk>/enriched-data.json", views.job_enriched_data, name="job_enriched_data"),
]

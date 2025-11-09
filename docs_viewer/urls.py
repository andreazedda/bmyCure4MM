"""
URL configuration for documentation viewer.
"""
from django.urls import path
from . import views

app_name = 'docs_viewer'

urlpatterns = [
    path('', views.docs_index, name='index'),
    path('view/<path:path>', views.docs_view, name='view'),
    path('raw/<path:path>', views.docs_raw, name='raw'),
    path('search/', views.docs_search, name='search'),
]

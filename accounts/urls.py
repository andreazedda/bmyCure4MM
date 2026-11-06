from django.urls import path

from . import views

app_name = "accounts"
urlpatterns = [
    path("access/request/", views.request_access, name="request"),
    path("access/request/submitted/", views.submitted, name="submitted"),
    path("access/request/status/", views.request_status, name="status"),
    path("access/verify/", views.capability, {"purpose": "verify"}, name="verify"),
    path("access/activate/", views.capability, {"purpose": "activate"}, name="activation"),
    path("owner/activate/", views.capability, {"purpose": "owner"}, name="owner_activation"),
    path("accounts/login/", views.sign_in, name="login"),
    path("accounts/logout/", views.sign_out, name="logout"),
    path("accounts/security/", views.owner_security, name="security"),
    path("password/reset/", views.reset_request, name="reset"),
    path("password/reset/confirm/", views.capability, {"purpose": "reset"}, name="reset_confirm"),
    path("feedback/", views.feedback, name="feedback"),
    path("manage/access-requests/", views.manage_requests, name="manage_requests"),
    path("manage/users/", views.manage_users, name="manage_users"),
    path("manage/feedback/", views.manage_feedback, name="manage_feedback"),
    path("manage/controls/", views.switches, name="switches"),
    path("manage/notifications/retry/", views.retry_notification, name="retry_notification"),
]

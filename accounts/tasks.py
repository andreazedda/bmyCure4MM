import logging

from celery import shared_task
from django.conf import settings

logger = logging.getLogger(__name__)


def request_delivery():
    if not settings.ACCOUNT_EMAIL_TASKS_ENABLED:
        return
    try:
        deliver_notifications.apply_async(retry=False)
    except Exception:
        logger.warning("account_notification_dispatch_queued")


@shared_task(
    ignore_result=True, autoretry_for=(Exception,), retry_kwargs={"max_retries": 5}, retry_backoff=10
)
def deliver_notifications():
    from .models import Notification
    from .notifications import dispatch_notification, housekeeping

    housekeeping()
    for _ in range(5):
        if not dispatch_notification():
            break
    if Notification.objects.filter(state__in=["queued", "failed"], attempts__lt=5).exists():
        deliver_notifications.apply_async(countdown=60, retry=False)

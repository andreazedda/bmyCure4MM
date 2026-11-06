from django.core.management.base import BaseCommand

from accounts.notifications import dispatch_notification, housekeeping


class Command(BaseCommand):
    help = "Deliver up to five committed registration or ticket notifications."

    def handle(self, **options):
        housekeeping()
        for _ in range(5):
            if not dispatch_notification():
                break

import json
import re
import secrets

from accounts import services
from accounts.models import AccessPolicy, Notification

TEST_SETTINGS = {
    "PUBLIC_BASE_URL": "https://member.example.test",
    "PLATFORM_OWNER_EMAIL": "owner@example.test",
    "PLATFORM_OWNER_USERNAME": "synthetic-owner",
    "PLATFORM_OWNER_DISPLAY_NAME": "Synthetic Owner",
}


def token_for(identity, purpose):
    item = Notification.objects.filter(identity=identity, kind=purpose).latest("created_at")
    body = json.loads(services.cipher().decrypt(item.encrypted_payload.encode()))["body"]
    return re.search(r"#token=([A-Za-z0-9_-]+)", body).group(1)


def open_policy():
    policy, _ = AccessPolicy.objects.get_or_create(pk=1)
    policy.requests_open = policy.activations_open = True
    policy.save()
    return policy


class LifecycleMixin:
    def initialize(self):
        self.password = secrets.token_urlsafe(24)
        self.policy = open_policy()
        self.owner, _ = services.bootstrap_owner()
        services.confirm_capability(token_for(self.owner, "owner"), "owner", self.password)
        self.owner.refresh_from_db()

    def make_member(self, email="alice@example.test", activate=True):
        identity = services.request_access(
            {"email": email, "nickname": email.split("@")[0], "reason": "Synthetic software verification"}
        )
        services.confirm_capability(token_for(identity, "verify"), "verify")
        identity = services.decide_request(self.owner.user, identity.pk, approve=True)
        if activate:
            services.confirm_capability(token_for(identity, "activate"), "activate", self.password)
        identity.refresh_from_db()
        return identity

    def sign_in(self, identity, client=None):
        client = client or self.client
        response = client.post("/accounts/login/", {"email": identity.email, "password": self.password})
        assert response.status_code == 302
        if identity.is_owner:
            import time

            import pyotp

            from accounts.services import cipher

            client.get("/accounts/security/")
            identity.refresh_from_db()
            secret = cipher().decrypt(identity.totp_secret.encode()).decode()
            step = max(int(time.time()) // 30, identity.totp_last_step + 1)
            reply = client.post("/accounts/security/", {"code": pyotp.TOTP(secret).at(step * 30)})
            assert reply.status_code == 302
        return response

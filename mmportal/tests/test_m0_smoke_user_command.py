from __future__ import annotations

import io
import sys
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase


class EnsureM0SmokeUserCommandTests(TestCase):
    username = "m0-smoke-contract"

    def _run_with_password(self, password: str = "test-only-password") -> None:
        with patch.object(sys, "stdin", io.StringIO(password + "\n")):
            call_command(
                "ensure_m0_smoke_user",
                "--username",
                self.username,
                "--password-stdin",
            )

    def test_requires_reserved_prefix(self):
        with self.assertRaises(CommandError):
            with patch.object(sys, "stdin", io.StringIO("test-only-password\n")):
                call_command(
                    "ensure_m0_smoke_user",
                    "--username",
                    "human-account",
                    "--password-stdin",
                )

    def test_create_rotate_and_disable_are_least_privilege_and_idempotent(self):
        self._run_with_password()
        user = get_user_model().objects.get(username=self.username)
        self.assertTrue(user.is_active)
        self.assertFalse(user.is_staff)
        self.assertFalse(user.is_superuser)
        self.assertEqual(user.email, f"{self.username}@m0-smoke.invalid")
        self.assertTrue(user.check_password("test-only-password"))
        self.assertEqual(user.groups.count(), 0)
        self.assertEqual(user.user_permissions.count(), 0)

        self._run_with_password("rotated-test-only-password")
        self.assertEqual(get_user_model().objects.filter(username=self.username).count(), 1)
        user.refresh_from_db()
        self.assertTrue(user.check_password("rotated-test-only-password"))

        call_command("ensure_m0_smoke_user", "--username", self.username, "--disable")
        user.refresh_from_db()
        self.assertFalse(user.is_active)
        self.assertFalse(user.has_usable_password())
        self.assertFalse(user.is_staff)
        self.assertFalse(user.is_superuser)

    def test_refuses_existing_unmarked_account(self):
        get_user_model().objects.create_user(
            username=self.username,
            email="person@example.invalid",
            password="human-test-password",
        )
        with self.assertRaises(CommandError):
            self._run_with_password()

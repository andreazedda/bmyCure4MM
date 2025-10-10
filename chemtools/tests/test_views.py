from __future__ import annotations

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.test import TestCase
from django.urls import reverse
from unittest.mock import patch

from chemtools import models


class RetryJobViewTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)

    def test_retry_job_shows_info_when_queued(self) -> None:
        job = models.ChemJob.objects.create(kind=models.ChemJob.SIM, input_a="CCC", user=self.user)
        with patch("chemtools.views._enqueue", return_value=(True, False)) as enqueue:
            response = self.client.get(reverse("chemtools:job_retry", args=[job.pk]))
            enqueue.assert_called_once()
        self.assertRedirects(response, reverse("chemtools:tools_home"))
        messages = list(get_messages(response.wsgi_request))
        self.assertTrue(any("queued" in message.message.lower() for message in messages))

    def test_retry_job_shows_error_when_immediate_failure(self) -> None:
        job = models.ChemJob.objects.create(kind=models.ChemJob.PARAM, input_a="CCO", input_b="123", user=self.user)
        with patch("chemtools.views._enqueue", return_value=(False, True)):
            response = self.client.get(reverse("chemtools:job_retry", args=[job.pk]))
        self.assertRedirects(response, reverse("chemtools:tools_home"))
        messages = list(get_messages(response.wsgi_request))
        self.assertTrue(any("failed" in message.message.lower() for message in messages))

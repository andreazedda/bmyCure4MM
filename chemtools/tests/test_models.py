from __future__ import annotations

from django.test import SimpleTestCase

from chemtools import models


class ChemJobStatusLabelTests(SimpleTestCase):
    def test_status_label_failed_when_log_has_error(self) -> None:
        job = models.ChemJob(kind=models.ChemJob.SIM, log="ERROR: something went wrong")
        self.assertEqual(job.status_label(), ("Failed", "danger"))

    def test_status_label_completed_when_outputs_present(self) -> None:
        job = models.ChemJob(kind=models.ChemJob.PARAM)
        job.out_html.name = "chem/1/drug_1_v1.html"
        self.assertEqual(job.status_label(), ("Completed", "success"))

    def test_status_label_default_queued(self) -> None:
        job = models.ChemJob(kind=models.ChemJob.BIND)
        self.assertEqual(job.status_label(), ("Queued", "secondary"))

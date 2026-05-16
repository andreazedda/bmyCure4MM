from __future__ import annotations

import tempfile
from datetime import date
from pathlib import Path

from django.contrib.auth import get_user_model
from django.test import Client, TestCase, override_settings
from django.urls import reverse

from clinic.models import Assessment, Patient
from twin_engine.counterfactual import run_counterfactual
from twin_engine.models import CounterfactualRun
from twin_engine.state_model import initialize_from_assessment


class TwinEngineAccessAndPrivacyTests(TestCase):
    def setUp(self) -> None:
        User = get_user_model()
        self.owner = User.objects.create_user("owner", password="pass1234")
        self.stranger = User.objects.create_user("stranger", password="pass1234")
        self.staff = User.objects.create_user("staff", password="pass1234", is_staff=True)
        self.client = Client()
        self.patient = Patient.objects.create(
            mrn="MM-PRIV-001",
            owner=self.owner,
            first_name="Privacy",
            last_name="Check",
            birth_date=date(1962, 7, 7),
            sex="F",
            diagnosis_date=date(2024, 3, 1),
            notes="Sensitive note that should not be exported.",
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 10),
            m_protein_g_dl=1.1,
            flc_ratio=2.2,
            hemoglobin_g_dl=11.9,
            r_iss="II",
            beta2m_mg_l=3.2,
            ldH_u_l=250,
        )

    def test_non_owner_cannot_access_patient_twin(self) -> None:
        self.client.login(username="stranger", password="pass1234")
        response = self.client.get(reverse("twin_engine:patient_twin_detail", args=[self.patient.id]))
        self.assertEqual(response.status_code, 403)

    def test_owner_can_access_own_patient_twin(self) -> None:
        self.client.login(username="owner", password="pass1234")
        response = self.client.get(reverse("twin_engine:patient_twin_detail", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)

    def test_staff_can_access_patient_twin(self) -> None:
        self.client.login(username="staff", password="pass1234")
        response = self.client.get(reverse("twin_engine:patient_twin_detail", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)

    def test_artifacts_do_not_contain_direct_identifiers_by_default(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.owner)
        intervention = {
            "drug_doses": {"lenalidomide": 25.0},
            "start_day": 1,
            "duration_days": 10,
            "schedule": {},
            "parameter_overrides": {},
            "random_seed": 11,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                run = run_counterfactual(self.patient, state, intervention, 10, user=self.owner)
                report_path = Path(tmpdir) / run.report_artifact.replace("/media/", "")
                content = report_path.read_text(encoding="utf-8")
                self.assertNotIn(self.patient.first_name, content)
                self.assertNotIn(self.patient.last_name, content)
                self.assertNotIn(self.patient.mrn, content)
                self.assertNotIn(self.patient.notes, content)

    def test_logs_do_not_include_direct_identifiers(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.owner)
        intervention = {
            "drug_doses": {"lenalidomide": 25.0},
            "start_day": 1,
            "duration_days": 10,
            "schedule": {},
            "parameter_overrides": {},
            "random_seed": 13,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                with self.assertLogs("twin_engine.research", level="INFO") as captured:
                    run_counterfactual(self.patient, state, intervention, 10, user=self.owner)
        output = "\n".join(captured.output)
        self.assertNotIn(self.patient.first_name, output)
        self.assertNotIn(self.patient.last_name, output)
        self.assertNotIn(self.patient.mrn, output)

    def test_legacy_counterfactual_report_without_predicted_biomarkers_still_renders(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.owner)
        legacy_run = CounterfactualRun.objects.create(
            patient=self.patient,
            base_twin_state=state,
            intervention_definition={"drug_doses": {"lenalidomide": 10.0}},
            simulation_summary={"label": "research simulation"},
            status=CounterfactualRun.STATUS_COMPLETED,
            created_by=self.owner,
        )
        self.client.login(username="owner", password="pass1234")
        response = self.client.get(reverse("twin_engine:counterfactual_report", args=[self.patient.id, legacy_run.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Unavailable for legacy run.")

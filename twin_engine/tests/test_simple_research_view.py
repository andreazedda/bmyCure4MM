from __future__ import annotations

from datetime import date

from django.contrib.auth import get_user_model
from django.test import Client, TestCase
from django.urls import reverse

from clinic.models import Patient, Regimen
from twin_engine.models import CounterfactualRun, LongitudinalLabResult
from twin_engine.state_model import initialize_from_assessment


class SimpleResearchViewTests(TestCase):
    def setUp(self) -> None:
        User = get_user_model()
        self.owner = User.objects.create_user("simple-owner", password="pass1234")
        self.staff = User.objects.create_user("simple-staff", password="pass1234", is_staff=True)
        self.client = Client()
        self.patient = Patient.objects.create(
            mrn="MM-SIMPLE-001",
            owner=self.owner,
            first_name="Private",
            last_name="Story",
            birth_date=date(1964, 1, 1),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = self.patient.assessments.create(
            date=date(2024, 6, 1),
            m_protein_g_dl=1.2,
            flc_ratio=2.6,
            hemoglobin_g_dl=11.4,
        )
        LongitudinalLabResult.objects.create(
            patient=self.patient,
            assessment=self.assessment,
            date=date(2024, 6, 1),
            analyte=LongitudinalLabResult.ANALYTE_M_PROTEIN,
            value=1.2,
            unit="g/dL",
        )
        LongitudinalLabResult.objects.create(
            patient=self.patient,
            assessment=self.assessment,
            date=date(2024, 6, 1),
            analyte=LongitudinalLabResult.ANALYTE_ALT,
            value=44,
            unit="U/L",
        )
        regimen = Regimen.objects.create(name="LEN", line="research", components="lenalidomide")
        self.patient.therapies.create(
            regimen=regimen,
            start_date=date(2024, 6, 2),
            doses={"lenalidomide": 10},
            schedule_notes="daily",
        )
        self.state = initialize_from_assessment(self.assessment, user=self.owner)
        CounterfactualRun.objects.create(
            patient=self.patient,
            base_twin_state=self.state,
            intervention_definition={"drug_doses": {"lenalidomide": 10.0}},
            simulation_summary={
                "label": "research simulation",
                "predicted_biomarkers": {"m_protein_g_dl": 0.9},
                "classification": {"counterfactual_class": "mechanistic"},
            },
            comparison_metrics={
                "research_utility_v2": 1.1,
                "tumor_reduction": {"baseline": 0.2, "alternative": 0.4},
                "healthy_loss": {"baseline": 0.2, "alternative": 0.25},
                "durability_index": {"baseline": 0.2, "alternative": 0.3},
            },
            status=CounterfactualRun.STATUS_COMPLETED,
            created_by=self.owner,
        )

    def test_simple_view_route_renders(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Simple Research View")
        self.assertContains(response, "Patient Research Summary")

    def test_simple_view_contains_required_sections(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        for text in [
            "Research simulation only",
            "Not a treatment recommendation",
            "Causal effect not identified",
            "One-minute summary",
            "What is missing?",
            "What the model actually uses",
            "What has been simulated?",
            "What changed across scenarios?",
            "What is uncertain?",
            "What can I conclude?",
            "What can I NOT conclude?",
            "What should I do next?",
            "Disease markers",
            "Blood / marrow markers",
            "Liver / toxicity markers",
            "Simulated scenarios",
        ]:
            self.assertContains(response, text, html=False)

    def test_simple_view_avoids_raw_json_and_unexplained_technical_terms(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, '"patient_id"')
        self.assertNotContains(response, "counterfactual")
        self.assertNotContains(response, "RMSE")
        self.assertNotContains(response, "MAE")
        self.assertNotContains(response, "do-operator")
        self.assertContains(response, "mathematical starting state")
        self.assertContains(response, "heuristic research score called utility_v2")

    def test_simple_view_contains_data_blocks_and_model_use_table(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Available count")
        self.assertContains(response, "Used by model")
        self.assertContains(response, "Input group")
        self.assertContains(response, "Interpretation risk")

    def test_simple_view_excludes_clinical_overclaims(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        for phrase in [
            " ".join(("recommended", "treatment")),
            " ".join(("best", "therapy")),
            " ".join(("should", "prescribe")),
            " ".join(("clinically", "superior")),
            " ".join(("proven", "causal", "effect")),
            " ".join(("validated", "clinical", "decision", "support")),
        ]:
            self.assertNotContains(response, phrase)

    def test_patient_twin_detail_alias_now_uses_simple_view(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:patient_twin_detail", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Simple Research View")

    def test_missing_data_is_explicit(self) -> None:
        missing_patient = Patient.objects.create(
            mrn="MM-SIMPLE-002",
            owner=self.owner,
            first_name="Missing",
            last_name="Data",
            birth_date=date(1965, 1, 1),
            sex="M",
            diagnosis_date=date(2024, 1, 1),
        )
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[missing_patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "No assessment history is available yet.")
        self.assertContains(response, "No structured treatment schedule is currently available.")

    def test_scenario_details_are_hidden_by_default(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Show technical metrics, including the heuristic research score called utility_v2")
        self.assertNotContains(response, "raw details")

    def test_scenario_table_uses_plain_language_interpretation(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Exposure pattern")
        self.assertContains(response, "Allowed conclusion")
        self.assertContains(response, "Forbidden conclusion")

    def test_developer_details_are_not_shown_by_default_for_non_staff(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertNotContains(response, reverse("twin_engine:developer_console"))

    def test_staff_navigation_includes_developer_console(self) -> None:
        self.client.force_login(self.staff)
        response = self.client.get(reverse("twin_engine:simple_research_view", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, reverse("twin_engine:developer_console"))

    def test_patient_page_points_to_simple_view(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("clinic:patient_detail", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Open Simple Research View")
        self.assertContains(response, reverse("twin_engine:simple_research_view", args=[self.patient.id]))

    def test_cockpit_points_back_to_simple_view(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:research_cockpit", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "This is the scientific cockpit.")
        self.assertContains(response, reverse("twin_engine:simple_research_view", args=[self.patient.id]))
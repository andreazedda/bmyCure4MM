from __future__ import annotations

from datetime import date
from pathlib import Path

from django.contrib.auth import get_user_model
from django.test import Client, TestCase
from django.urls import reverse

from clinic.models import Assessment, Patient
from twin_engine.models import LongitudinalLabResult
from twin_engine.state_model import initialize_from_assessment


class ResearchContentClarityTests(TestCase):
    def setUp(self) -> None:
        User = get_user_model()
        self.owner = User.objects.create_user("clarity-owner", password="pass1234")
        self.staff = User.objects.create_user("clarity-staff", password="pass1234", is_staff=True)
        self.patient = Patient.objects.create(
            mrn="MM-CLARITY-001",
            owner=self.owner,
            first_name="Private",
            last_name="Researcher",
            birth_date=date(1965, 1, 1),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2024, 5, 20),
            m_protein_g_dl=1.2,
            flc_ratio=2.7,
            hemoglobin_g_dl=11.4,
            creatinine_mg_dl=0.8,
        )
        LongitudinalLabResult.objects.create(
            patient=self.patient,
            assessment=self.assessment,
            date=date(2024, 5, 20),
            analyte=LongitudinalLabResult.ANALYTE_M_PROTEIN,
            value=1.2,
            unit="g/dL",
        )
        initialize_from_assessment(self.assessment, user=self.owner)
        self.client = Client()

    def test_cockpit_explains_major_sections_and_formulas(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:research_cockpit", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        for text in [
            "Workflow map",
            "How to use this page",
            "Question answered",
            "Inputs used",
            "Method / computation",
            "Output meaning",
            "What you can do next",
            "Limitations",
            "Developer details",
            "residual = observed value - predicted value",
            "RMSE = root mean square error",
            "MAE = mean absolute error",
            "research_utility = tumor_reduction + (1 - healthy_loss) + durability_index - toxicity_constraint_penalty",
            "research_utility_v2 = research_utility - 0.5 * liver_toxicity_signal_0_1 - 0.5 * neutropenia_signal_0_1",
            "Y_model(a') = f(x_t, theta_hat, a')",
            "E[Y | do(A=a')] - E[Y | do(A=a)]",
            "Causal effect not identified",
            "descriptive_only means observed penalty only",
            "Prototype toxicity signals do not claim clinical validity",
            "Schedule classification",
            "Utility v2",
            "Validation, uncertainty, and robustness",
            "probability-best ranking under aligned uncertainty samples",
            "No peer-reviewed source is attached in the repository for this component",
            "Traceability chain",
        ]:
            self.assertContains(response, text, html=False)

    def test_initialize_page_explains_scope_and_missing_fields(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("twin_engine:initialize_twin_from_assessment", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        for text in [
            "What initialization is",
            "What initialization is not",
            "Assessment(t0) -> risk mapping -> tumor/healthy/immune parameters -> PatientTwinState",
            "Not a fitted calibration result",
            "Field meanings",
            "Recommended because it has the strongest modeled-field completeness",
            "Missing fields reduce model interpretability",
        ]:
            self.assertContains(response, text, html=False)

    def test_patient_detail_explains_workspace_and_next_step(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(reverse("clinic:patient_detail", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        for text in [
            "Patient workspace overview",
            "This page is for data navigation.",
            "Use Simple Research View to understand the data and model.",
            "Use Scientific Cockpit for technical details.",
            "Recommended next step",
            "Research simulation only",
            "not clinically validated",
            "Twin Inputs (over time)",
            "observed inputs, not model predictions",
            "Missing chart lines mean missing structured records",
        ]:
            self.assertContains(response, text, html=False)

    def test_developer_console_and_glossary_are_explanatory(self) -> None:
        self.client.force_login(self.staff)
        console = self.client.get(reverse("twin_engine:developer_console"))
        self.assertEqual(console.status_code, 200)
        for text in [
            "How to use this console",
            "Research navigation layers",
            "Simple Research View",
            "Scientific Cockpit",
            "Question answered",
            "Data checks inspect structured assessments",
            "Model checks inspect twin state",
            "Causal checks inspect assumption-set documentation",
            "Scientific checks inspect model component metadata",
            "Privacy checks inspect tracked/staged paths",
            "raw details",
            "audit signal",
        ]:
            self.assertContains(console, text, html=False)

        glossary = self.client.get(reverse("twin_engine:research_glossary"))
        self.assertEqual(glossary.status_code, 200)
        for text in [
            "Research Glossary",
            "Twin",
            "Calibration",
            "Residual",
            "Counterfactual",
            "Causal effect",
            "Toxicity constraint",
            "Research utility",
            "Exposure bridge",
            "Schedule collapse",
            "Provenance",
        ]:
            self.assertContains(glossary, text)

    def test_forbidden_clinical_claim_phrases_absent_from_templates_and_docs(self) -> None:
        forbidden = [
            " ".join(("recommended", "treatment")),
            " ".join(("best", "therapy")),
            " ".join(("should", "prescribe")),
            " ".join(("clinically", "superior")),
            " ".join(("proven", "causal", "effect")),
            " ".join(("validated", "clinical", "decision", "support")),
        ]
        paths = [
            *Path("twin_engine/templates").rglob("*.html"),
            *Path("clinic/templates").rglob("*.html"),
            *Path("docs/research").rglob("*.md"),
            *Path("docs/development").rglob("*.md"),
        ]
        for path in paths:
            text = path.read_text(encoding="utf-8")
            lowered = text.lower()
            for phrase in forbidden:
                self.assertNotIn(phrase, lowered, msg=f"{phrase!r} found in {path}")

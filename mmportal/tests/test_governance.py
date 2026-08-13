from __future__ import annotations

import json
import re
from pathlib import Path

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from chemtools.pdb_api_client import PDBAPIClient
from mmportal import governance
from simulator.models_help import HelpArticle
from simulator.presets import PRESETS
from simulator.regimen_suggester import suggest_regimens


class ReleaseClaimsContractTests(TestCase):
    def test_machine_readable_claims_match_runtime_constants(self) -> None:
        manifest_path = Path(__file__).resolve().parents[2] / "governance" / "release_claims.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.assertEqual(manifest["intended_use_level"], governance.INTENDED_USE_LEVEL)
        self.assertEqual(manifest["clinical_decision_support"], governance.CLINICAL_DECISION_SUPPORT)
        self.assertEqual(
            manifest["patient_specific_prediction_validated"],
            governance.PATIENT_SPECIFIC_PREDICTION_VALIDATED,
        )
        self.assertEqual(manifest["causal_effect_identified"], governance.CAUSAL_EFFECT_IDENTIFIED)
        self.assertEqual(manifest["model_version"], governance.CURRENT_RESEARCH_MODEL_VERSION)
        self.assertEqual(
            manifest["epistemic_labels"],
            [label.value for label in governance.EpistemicLabel],
        )
        self.assertEqual(
            manifest["forbidden_claim_codes"],
            list(governance.FORBIDDEN_CLAIM_CODES),
        )

    def test_model_relative_diagnostics_are_non_prescriptive(self) -> None:
        result = governance.build_model_relative_diagnostics(
            {
                "tumor_reduction": 0.1,
                "healthy_loss": 0.4,
                "time_to_recurrence": 30,
            }
        )
        self.assertIsNotNone(result)
        assert result is not None
        metadata = result["governance"]
        self.assertFalse(metadata["clinical_decision_support"])
        self.assertFalse(metadata["patient_specific_prediction_validated"])
        self.assertFalse(metadata["causal_effect_identified"])
        serialized = json.dumps(result).lower()
        for phrase in ("best treatment", "recommended therapy", "switch regimen", "reduce dose"):
            self.assertNotIn(phrase, serialized)

    def test_active_sources_contain_no_actionable_clinical_claims(self) -> None:
        root = Path(__file__).resolve().parents[2]
        source_roots = [
            root / "mmportal",
            root / "clinic",
            root / "simulator",
            root / "twin_engine",
            root / "chemtools",
            root / "templates",
            root / "docs",
        ]
        excluded_parts = {"tests", "migrations", "__pycache__", "archive", "governance"}
        allowed_suffixes = {".py", ".html", ".js", ".json", ".md", ".yaml", ".yml"}
        patterns = [
            re.compile(pattern, re.IGNORECASE)
            for pattern in (
                r"\bbest treatment\b",
                r"\brecommended treatment\b",
                r"\brecommended therapy\b",
                r"\bclinically superior\b",
                r"\bpatient should (?:change|increase|decrease|reduce|switch|receive)\b",
                r"\bswitch to .{0,50} regimen\b",
                r"\brequires? (?:aggressive|intensive) (?:treatment|therapy)\b",
                r"\bsafe (?:dose|doses|preset|presets)\b",
                r"\bstrong evidence for mm efficacy\b",
                r"\bstandard of care\b",
                r"\bproven effective\b",
                r"\bbest weapon\b",
                r"\bmiglior trattamento\b",
                r"\btrattamento raccomandato\b",
                r"\bterapia raccomandata\b",
                r"\bclinicamente superiore\b",
                r"\briduci(?:re)? (?:la |le )?dos[ei]\b",
                r"\baumenta(?:re)? (?:la |le )?dos[ei]\b",
                r"\bcambia(?:re)? (?:il )?regime\b",
                r"\brichiede (?:un )?trattamento aggressivo\b",
                r"\bdosi sicure\b",
                r"\bstandard di cura\b",
                r"\bprovato efficace\b",
                r"\barma migliore\b",
            )
        ]

        violations: list[str] = []
        for source_root in source_roots:
            for path in source_root.rglob("*"):
                if not path.is_file() or path.suffix not in allowed_suffixes:
                    continue
                if excluded_parts.intersection(path.relative_to(root).parts):
                    continue
                text = path.read_text(encoding="utf-8", errors="replace")
                for pattern in patterns:
                    if pattern.search(text):
                        violations.append(f"{path.relative_to(root)}: {pattern.pattern}")

        self.assertEqual(violations, [])


class GovernedRuntimeTests(TestCase):
    def setUp(self) -> None:
        user_model = get_user_model()
        self.user = user_model.objects.create_user(username="governance-user", password="pass")
        self.client.force_login(self.user)

    def test_global_banner_and_learning_alias_are_available(self) -> None:
        response = self.client.get(reverse("clinic:dashboard"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "E1 research prototype")
        self.assertContains(response, "Clinical decision support: false")
        self.assertEqual(reverse("learn:scenario_list"), "/learn/")
        self.assertEqual(reverse("twin_engine:research_home"), "/research/")
        research_response = self.client.get(reverse("twin_engine:research_home"))
        self.assertEqual(research_response.status_code, 200)
        self.assertContains(research_response, "lineage-bound workspace")

    def test_major_reference_apis_carry_governance(self) -> None:
        glossary = self.client.get(reverse("api_glossary"))
        self.assertEqual(glossary.status_code, 200)
        self.assertEqual(glossary.json()["governance"]["epistemic_label"], "DERIVED")
        self.assertEqual(glossary.json()["governance"]["model_version"], "research-twin-v1")

        profiles = self.client.get(reverse("api_drugs"))
        self.assertEqual(profiles.status_code, 200)
        self.assertEqual(profiles.json()["governance"]["epistemic_label"], "HYPOTHETICAL")
        self.assertEqual(profiles.json()["governance"]["model_version"], "research-twin-v1")

    def test_runtime_presets_do_not_make_clinical_selection_claims(self) -> None:
        serialized = json.dumps(PRESETS).lower()
        for phrase in (
            "standard of care",
            "best weapon",
            "proven effective",
            "monitor blood pressure",
            "dose reduction",
        ):
            self.assertNotIn(phrase, serialized)

    def test_governed_help_overrides_stale_prescriptive_database_copy(self) -> None:
        HelpArticle.objects.update_or_create(
            slug="kpi_healthy_loss",
            defaults={
                "title_en": "Legacy title",
                "body_en": "Reduce doses when the badge is red.",
                "title_it": "Titolo legacy",
                "body_it": "Riduci la dose.",
            },
        )
        response = self.client.get(reverse("api_help_item", args=["kpi_healthy_loss"]))
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertNotIn("reduce doses", payload["body"].lower())
        self.assertEqual(payload["governance"]["intended_use_level"], "E1_research_prototype")
        self.assertFalse(payload["governance"]["clinical_decision_support"])

    def test_regimen_catalog_output_is_heuristic_and_non_prescriptive(self) -> None:
        result = suggest_regimens(
            transplant_eligible=False,
            ecog=2,
            high_risk_cytogenetics=True,
            line_of_therapy=2,
            refractory_to=["lenalidomide"],
            renal_function="dialysis",
            neuropathy_grade=3,
        )
        self.assertEqual(result["governance"]["epistemic_label"], "HEURISTIC")
        self.assertEqual(result["governance"]["evidence_status"], "NEEDS_EVIDENCE")
        self.assertEqual(result["governance"]["model_version"], "research-twin-v1")
        self.assertEqual(result["governance"]["claims_policy_version"], "claims-policy-v1")
        self.assertFalse(result["governance"]["clinical_decision_support"])
        serialized = json.dumps(result).lower()
        for phrase in (
            "standard of care",
            "switch to",
            "dose reduction",
            "preferred for",
            "treatment decisions",
        ):
            self.assertNotIn(phrase, serialized)

    def test_unsupported_chemistry_estimators_fail_closed(self) -> None:
        client = PDBAPIClient()
        efficacy = client.estimate_mm_efficacy("candidate", {}, [], [])
        survival = client.estimate_survival_impact({}, {})
        toxicity = client.estimate_toxicity_profile({}, {})

        self.assertEqual(efficacy["status"], "INSUFFICIENT_EVIDENCE")
        self.assertIsNone(efficacy["overall_score"])
        self.assertEqual(survival["status"], "PATIENT_SPECIFIC_PREDICTION_NOT_VALIDATED")
        self.assertIsNone(survival["median_pfs_months"])
        self.assertIsNone(survival["median_os_months"])
        self.assertEqual(toxicity["status"], "INSUFFICIENT_EVIDENCE")
        self.assertIsNone(toxicity["risk_score"])
        self.assertIsNone(toxicity["risk_benefit_ratio"])

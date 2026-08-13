from __future__ import annotations

import json
import tempfile
from datetime import date
from pathlib import Path

from django.contrib.auth import get_user_model
from django.core.exceptions import ValidationError
from django.test import TestCase, override_settings
from django.urls import reverse

from clinic.models import Assessment, Patient
from twin_engine.comparability import (
    DIRECTLY_COMPARABLE,
    NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN,
    compare_manifests,
    invalidate_manifest,
)
from twin_engine.contracts import PARAMETERS_SCHEMA_VERSION
from twin_engine.counterfactual import run_counterfactual
from twin_engine.model_registry import (
    get_model_registration,
    load_model_registry,
    validate_registry_entry_points,
)
from twin_engine.models import ResearchRunInvalidation
from twin_engine.provenance import record_simulation_metadata
from twin_engine.report_builder import write_json_artifact
from twin_engine.run_manifest import verify_manifest_artifacts
from twin_engine.state_model import initialize_from_assessment


MANDATORY_VERSION_VECTOR = {
    "app_version",
    "api_version",
    "db_schema_version",
    "dataset_id",
    "dataset_version",
    "dataset_sha256",
    "record_subset_sha256",
    "model_id",
    "model_version",
    "model_card_version",
    "configuration_sha256",
    "evidence_graph_version",
    "validation_protocol_version",
    "report_template_version",
    "git_sha",
    "container_digest",
    "random_seed",
    "created_at",
    "intended_use_level",
    "epistemic_label",
}


class ResearchRunManifestTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("manifest-researcher", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="SYN-MANIFEST-001",
            owner=self.user,
            first_name="Synthetic",
            last_name="Manifest",
            birth_date=date(1970, 1, 1),
            sex="O",
            diagnosis_date=date(2025, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 2, 1),
            m_protein_g_dl=1.2,
            flc_ratio=2.4,
            hemoglobin_g_dl=11.8,
            r_iss="II",
            beta2m_mg_l=3.0,
            ldH_u_l=220,
        )

    def test_registry_contains_all_governed_models_and_rejects_unknown_ids(self) -> None:
        self.assertEqual(
            set(load_model_registry()),
            {
                "patient_twin_state_model",
                "observation_model",
                "lenalidomide_exposure_model",
                "hepatic_signal_model",
                "neutropenia_signal_model",
                "counterfactual_model",
                "therapy_design_toy_model",
            },
        )
        with self.assertRaises(ValidationError):
            get_model_registration("unregistered_model")
        validate_registry_entry_points()

    def test_scientific_run_has_complete_immutable_version_vector(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, override_settings(MEDIA_ROOT=tmpdir):
            state = initialize_from_assessment(self.assessment, user=self.user)
            metadata = state.metadata_records.get(solver_name="initial_state_mapping")
            manifest = metadata.manifest
            self.assertIsNotNone(manifest)
            self.assertTrue(MANDATORY_VERSION_VECTOR.issubset(manifest.version_vector()))
            self.assertEqual(manifest.model_id, "patient_twin_state_model")
            self.assertEqual(manifest.dataset_id, "UNBOUND")
            self.assertEqual(len(manifest.record_subset_sha256), 64)
            self.assertEqual(
                compare_manifests(manifest, manifest).status,
                NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN,
            )
            verify_manifest_artifacts(manifest)
            manifest.status = "failed"
            with self.assertRaises(ValidationError):
                manifest.save()
            metadata.solver_name = "changed"
            with self.assertRaises(ValidationError):
                metadata.save()

    def test_artifact_tampering_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
            state = initialize_from_assessment(self.assessment, user=self.user)
            _, artifact_path = write_json_artifact(
                "tamper_test",
                {"value": 1},
                patient_id=self.patient.id,
                run_id=1,
                folder_name="research_reports",
            )
            metadata = record_simulation_metadata(
                twin_state=state,
                model_id="patient_twin_state_model",
                solver_name="artifact_integrity_test",
                input_payload={"value": 1},
                output_payload={"value": 2},
                artifact_paths=[("test_report", "output", artifact_path)],
            )
            verify_manifest_artifacts(metadata.manifest)
            Path(artifact_path).write_text('{"value": 3}', encoding="utf-8")
            with self.assertRaisesMessage(ValidationError, "Artifact integrity check failed"):
                verify_manifest_artifacts(metadata.manifest)

    def test_schema_validation_rejects_malformed_typed_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, override_settings(MEDIA_ROOT=tmpdir):
            state = initialize_from_assessment(self.assessment, user=self.user)
            state.parameters_schema_version = PARAMETERS_SCHEMA_VERSION
            state.parameters = {"tumor_growth_rate": "not-a-number"}
            with self.assertRaises(ValidationError):
                state.save()

    def test_comparability_and_invalidation_are_explicit_and_append_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir, override_settings(MEDIA_ROOT=tmpdir):
            state = initialize_from_assessment(self.assessment, user=self.user)
            lineage = dict(state.lineage)
            lineage["dataset_binding"] = {
                "status": "bound",
                "dataset_id": "synthetic-comparability",
                "dataset_version": "1.0.0",
                "canonical_dataset_sha256": "d" * 64,
            }
            state.lineage = lineage
            state.save(update_fields=["lineage"])
            first = record_simulation_metadata(
                twin_state=state,
                model_id="patient_twin_state_model",
                solver_name="comparison_test",
                input_payload={"case": "same"},
                solver_parameters={"step": 1.0},
                output_payload={"result": 1.0},
            ).manifest
            second = record_simulation_metadata(
                twin_state=state,
                model_id="patient_twin_state_model",
                solver_name="comparison_test",
                input_payload={"case": "same"},
                solver_parameters={"step": 1.0},
                output_payload={"result": 1.0},
            ).manifest
            self.assertEqual(compare_manifests(first, second).status, DIRECTLY_COMPARABLE)
            invalidation = invalidate_manifest(
                first,
                change_kind=ResearchRunInvalidation.CHANGE_DATASET,
                replacement_identity="f" * 64,
                reason="A corrected immutable dataset version supersedes this binding.",
            )
            self.assertEqual(invalidation.previous_identity, "d" * 64)
            result = compare_manifests(first, second)
            self.assertEqual(result.status, NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN)
            self.assertIn("invalidation", " ".join(result.reasons))
            invalidation.reason = "changed"
            with self.assertRaises(ValidationError):
                invalidation.save()

    def test_report_identity_matches_manifest_and_tampered_report_returns_409(self) -> None:
        intervention = {
            "drug_doses": {"lenalidomide": 5.0},
            "start_day": 0,
            "duration_days": 7,
            "schedule": {},
            "parameter_overrides": {},
            "random_seed": 7,
        }
        with tempfile.TemporaryDirectory() as tmpdir, override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
            state = initialize_from_assessment(self.assessment, user=self.user)
            run = run_counterfactual(self.patient, state, intervention, 7, user=self.user)
            metadata = run.metadata_records.select_related("manifest").get()
            report_path = Path(tmpdir) / run.report_artifact.removeprefix("/media/")
            payload = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["provenance"]["run_identity"]["run_id"], metadata.manifest.run_id)
            self.assertEqual(
                payload["provenance"]["run_identity"]["version_vector"],
                metadata.manifest.version_vector(),
            )
            self.client.force_login(self.user)
            response = self.client.get(reverse("twin_engine:counterfactual_report", args=[self.patient.id, run.id]))
            self.assertEqual(response.status_code, 200)
            report_path.write_text('{"tampered": true}', encoding="utf-8")
            response = self.client.get(reverse("twin_engine:counterfactual_report", args=[self.patient.id, run.id]))
            self.assertEqual(response.status_code, 409)
            self.assertContains(response, "ARTIFACT_INTEGRITY_FAILED", status_code=409)

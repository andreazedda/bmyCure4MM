from __future__ import annotations

import hashlib
import json
import tempfile
from datetime import date
from io import StringIO
from pathlib import Path

from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase

from clinic.models import Patient
from twin_engine.models import (
    AdverseEvent,
    LongitudinalLabResult,
    TherapyInterruption,
)
from twin_engine.toxicity_model import (
    compute_toxicity_constraints,
)


class ResearchDatasetImportCommandTests(TestCase):
    def setUp(self) -> None:
        self.patient = Patient.objects.create(
            mrn="SYNTHETIC-IMPORT-001",
            first_name="Synthetic",
            last_name="ResearchCase",
            birth_date=date(1970, 1, 1),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )

        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

        root = Path(self.tmp.name)

        self.dataset_path = root / "dataset.json"
        self.manifest_path = root / "manifest.json"

        self.dataset = self._make_dataset()
        self._write_bundle(self.dataset)

    @staticmethod
    def _assertion(record_id: str) -> list[dict]:
        return [
            {
                "assertion_id": f"ASSERT_{record_id}",
                "source_id": "SRC_SYNTHETIC",
                "source_locator": "synthetic-test-fixture",
                "source_sha256": None,
                "extraction_method": "manual_verified",
                "technical_validation_status": "validated",
                "clinical_validation_status": "not_required",
                "source_quality": "synthetic_test",
                "epistemic_class": "observed",
                "notes": "Synthetic test evidence only.",
            }
        ]

    def _record(
        self,
        record_id: str,
        record_type: str,
        identity: dict,
        payload: dict,
    ) -> dict:
        return {
            "record_id": record_id,
            "record_type": record_type,
            "identity": {
                "case_ref": "MM_RESEARCH_CASE_999",
                **identity,
            },
            "payload": payload,
            "provenance": self._assertion(record_id),
        }

    def _make_dataset(self) -> dict:
        records = [
            self._record(
                "LAB_20250507_AST",
                "lab_result",
                {
                    "date": "2031-04-01",
                    "analyte": "AST",
                },
                {
                    "date": "2031-04-01",
                    "analyte": "AST",
                    "value": 123,
                    "unit": "U/L",
                    "source_quality": "extracted_document",
                },
            ),
            self._record(
                "LAB_20250507_ALT",
                "lab_result",
                {
                    "date": "2031-04-01",
                    "analyte": "ALT",
                },
                {
                    "date": "2031-04-01",
                    "analyte": "ALT",
                    "value": 87,
                    "unit": "U/L",
                    "source_quality": "extracted_document",
                },
            ),
            self._record(
                "LAB_20250110_NEU",
                "lab_result",
                {
                    "date": "2031-03-15",
                    "analyte": "NEU",
                },
                {
                    "date": "2031-03-15",
                    "analyte": "NEU",
                    "value": 700,
                    "unit": "/uL",
                    "source_quality": "extracted_document",
                },
            ),
            self._record(
                "AE_20250110_NEUTROPENIA",
                "adverse_event",
                {
                    "date": "2031-03-15",
                    "event_type": "neutropenia",
                },
                {
                    "date": "2031-03-15",
                    "event_type": "neutropenia",
                    "grade": "",
                    "suspected_drug": "",
                    "observed_values": {
                        "absolute_neutropenia": True
                    },
                    "action_taken": "Synthetic support event.",
                    "outcome": "Synthetic recovery.",
                    "source_quality": "clinical_record",
                },
            ),
            self._record(
                "AE_20250110_PNEUMONIA",
                "adverse_event",
                {
                    "date": "2031-03-15",
                    "event_type": "pneumonia",
                },
                {
                    "date": "2031-03-15",
                    "event_type": "pneumonia",
                    "grade": "",
                    "suspected_drug": "",
                    "observed_values": {
                        "left_basal_pneumonia": True
                    },
                    "action_taken": "Synthetic anti-infective event.",
                    "outcome": "Synthetic treatment completed.",
                    "source_quality": "clinical_record",
                },
            ),
            self._record(
                "AE_20250507_TRANSAMINASE",
                "adverse_event",
                {
                    "date": "2031-04-01",
                    "event_type": "hypertransaminasemia",
                },
                {
                    "date": "2031-04-01",
                    "event_type": "hypertransaminasemia",
                    "grade": "",
                    "suspected_drug": "lenalidomide",
                    "observed_values": {
                        "AST": 123,
                        "ALT": 87,
                    },
                    "action_taken": "Synthetic interruption.",
                    "outcome": "Synthetic improvement.",
                    "source_quality": "clinical_record",
                },
            ),
            self._record(
                "INT_20250507_LENALIDOMIDE",
                "therapy_interruption",
                {
                    "drug": "lenalidomide",
                    "start_date": "2031-04-01",
                    "end_date": "2031-04-05",
                },
                {
                    "drug": "lenalidomide",
                    "start_date": "2031-04-01",
                    "end_date": "2031-04-05",
                    "reason": "hypertransaminasemia",
                    "evidence": {
                        "AST": 123,
                        "ALT": 87,
                    },
                    "action_taken": "Synthetic interruption.",
                    "source_quality": "clinical_record",
                },
            ),
        ]

        return {
            "schema_version": "0.1.0",
            "dataset_id": "SYNTHETIC_IMPORT_TEST",
            "dataset_version": "0.1.0",
            "case_ref": "MM_RESEARCH_CASE_999",
            "created_at": "2026-07-23T12:00:00Z",
            "sources": [
                {
                    "source_id": "SRC_SYNTHETIC",
                    "source_type": "synthetic_fixture",
                    "sha256": None,
                    "description": "Synthetic importer test source.",
                }
            ],
            "records": records,
        }

    def _write_bundle(self, dataset: dict) -> None:
        self.dataset_path.write_text(
            json.dumps(
                dataset,
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

        canonical = json.dumps(
            dataset,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")

        canonical_sha = hashlib.sha256(
            canonical
        ).hexdigest()

        file_sha = hashlib.sha256(
            self.dataset_path.read_bytes()
        ).hexdigest()

        counts = {}

        for record in dataset["records"]:
            kind = record["record_type"]
            counts[kind] = counts.get(kind, 0) + 1

        self.manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": "0.1.0",
                    "dataset_id": dataset["dataset_id"],
                    "dataset_version": dataset[
                        "dataset_version"
                    ],
                    "case_ref": dataset["case_ref"],
                    "canonical_dataset_sha256": canonical_sha,
                    "dataset_file_sha256": file_sha,
                    "record_counts": {
                        **counts,
                        "total": len(dataset["records"]),
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

    def _run(
        self,
        *,
        dry_run: bool = False,
        expect_error: bool = False,
    ):
        output = StringIO()

        kwargs = {
            "patient_id": self.patient.id,
            "dataset": str(self.dataset_path),
            "manifest": str(self.manifest_path),
            "dry_run": dry_run,
            "stdout": output,
        }

        if expect_error:
            with self.assertRaises(CommandError):
                call_command(
                    "import_research_dataset",
                    **kwargs,
                )

            return output

        call_command(
            "import_research_dataset",
            **kwargs,
        )

        return json.loads(output.getvalue())

    def test_dry_run_classifies_without_mutation(self) -> None:
        result = self._run(dry_run=True)

        self.assertEqual(result["created"], 7)
        self.assertEqual(result["changed"], 0)
        self.assertEqual(result["unchanged"], 0)
        self.assertEqual(result["conflicts"], 0)
        self.assertFalse(result["applied"])

        self.assertEqual(
            LongitudinalLabResult.objects.filter(
                patient=self.patient
            ).count(),
            0,
        )
        self.assertEqual(
            AdverseEvent.objects.filter(
                patient=self.patient
            ).count(),
            0,
        )
        self.assertEqual(
            TherapyInterruption.objects.filter(
                patient=self.patient
            ).count(),
            0,
        )

    def test_second_import_is_semantic_noop(self) -> None:
        first = self._run()
        second = self._run()

        self.assertEqual(first["created"], 7)
        self.assertEqual(first["changed"], 0)
        self.assertEqual(first["conflicts"], 0)
        self.assertTrue(first["applied"])

        self.assertEqual(second["created"], 0)
        self.assertEqual(second["changed"], 0)
        self.assertEqual(second["unchanged"], 7)
        self.assertEqual(second["conflicts"], 0)
        self.assertTrue(second["applied"])

        self.assertEqual(
            LongitudinalLabResult.objects.filter(
                patient=self.patient
            ).count(),
            3,
        )
        self.assertEqual(
            AdverseEvent.objects.filter(
                patient=self.patient
            ).count(),
            3,
        )
        self.assertEqual(
            TherapyInterruption.objects.filter(
                patient=self.patient
            ).count(),
            1,
        )

    def test_import_persists_dataset_provenance(self) -> None:
        self._run()

        lab = LongitudinalLabResult.objects.get(
            patient=self.patient,
            analyte="AST",
        )

        self.assertEqual(
            lab.provenance["dataset_id"],
            "SYNTHETIC_IMPORT_TEST",
        )
        self.assertEqual(
            lab.provenance["dataset_version"],
            "0.1.0",
        )
        self.assertEqual(
            len(lab.provenance["assertions"]),
            1,
        )

        interruption = TherapyInterruption.objects.get(
            patient=self.patient
        )

        self.assertIn(
            "_dataset_provenance",
            interruption.evidence,
        )

        self.assertEqual(
            interruption.evidence[
                "_dataset_provenance"
            ]["record_id"],
            "INT_20250507_LENALIDOMIDE",
        )

    def test_semantic_conflict_blocks_overwrite(self) -> None:
        self._run()

        ast = LongitudinalLabResult.objects.get(
            patient=self.patient,
            analyte="AST",
        )

        ast.value = 999.0
        ast.save(update_fields=["value"])

        output = self._run(expect_error=True)

        payload = json.loads(output.getvalue())

        self.assertEqual(payload["conflicts"], 1)
        self.assertFalse(payload["applied"])

        ast.refresh_from_db()

        self.assertEqual(ast.value, 999.0)

    def test_manifest_detects_physical_dataset_tampering(
        self,
    ) -> None:
        mutated = json.loads(
            self.dataset_path.read_text()
        )

        mutated["records"][0]["payload"]["value"] = 999

        self.dataset_path.write_text(
            json.dumps(
                mutated,
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )

        self._run(expect_error=True)

        self.assertEqual(
            LongitudinalLabResult.objects.filter(
                patient=self.patient
            ).count(),
            0,
        )

    def test_legacy_backfill_command_is_dataset_alias(self) -> None:
        output = StringIO()

        call_command(
            "backfill_research_structured_data",
            patient_id=self.patient.id,
            dataset=str(self.dataset_path),
            manifest=str(self.manifest_path),
            dry_run=True,
            stdout=output,
        )

        payload = json.loads(output.getvalue())

        self.assertEqual(payload["records_total"], 7)
        self.assertEqual(payload["created"], 7)
        self.assertEqual(payload["changed"], 0)
        self.assertEqual(payload["unchanged"], 0)
        self.assertEqual(payload["conflicts"], 0)
        self.assertFalse(payload["applied"])

        self.assertEqual(
            LongitudinalLabResult.objects.filter(
                patient=self.patient
            ).count(),
            0,
        )

    def test_toxicity_layer_uses_imported_records(self) -> None:
        self._run()

        summary = compute_toxicity_constraints(
            self.patient
        )

        self.assertEqual(
            summary["liver"]["max_ast"],
            123.0,
        )
        self.assertEqual(
            summary["liver"]["max_alt"],
            87.0,
        )
        self.assertEqual(
            summary["liver"]["peak_date"],
            "2031-04-01",
        )
        self.assertTrue(
            summary["neutropenia"][
                "absolute_neutropenia_event"
            ]
        )
        self.assertTrue(
            summary["infection"]["pneumonia_event"]
        )

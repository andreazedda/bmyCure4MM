from __future__ import annotations

from datetime import date

from django.db import connection
from django.db.migrations.executor import MigrationExecutor
from django.test import TransactionTestCase


class ResearchRunManifestMigrationTests(TransactionTestCase):
    migrate_from = ("twin_engine", "0003_patienttwinstate_lineage")
    migrate_to = ("twin_engine", "0004_research_run_manifests")

    def setUp(self) -> None:
        super().setUp()
        self.addCleanup(self._migrate_latest)
        executor = MigrationExecutor(connection)
        executor.migrate([self.migrate_from])
        old_apps = executor.loader.project_state([self.migrate_from]).apps

        Patient = old_apps.get_model("clinic", "Patient")
        PatientTwinState = old_apps.get_model("twin_engine", "PatientTwinState")
        CounterfactualRun = old_apps.get_model("twin_engine", "CounterfactualRun")
        ObservationResidual = old_apps.get_model("twin_engine", "ObservationResidual")
        SimulationRunMetadata = old_apps.get_model("twin_engine", "SimulationRunMetadata")

        patient = Patient.objects.create(
            mrn="SYN-MIGRATION-001",
            first_name="Synthetic",
            last_name="Migration",
            birth_date=date(1970, 1, 1),
            sex="O",
            diagnosis_date=date(2025, 1, 1),
        )
        state = PatientTwinState.objects.create(
            patient=patient,
            state_date=date(2025, 2, 1),
            is_current=True,
            state_vector={"historical": True},
            parameters={"historical": True},
            parameter_uncertainty={"historical": True},
            model_version="historical-model",
            config_hash="historical-config",
            lineage={"historical": True},
        )
        run = CounterfactualRun.objects.create(
            patient=patient,
            base_twin_state=state,
            intervention_definition={"historical": True},
            simulation_summary={"historical": True},
            comparison_metrics={"historical": True},
            status="completed",
        )
        ObservationResidual.objects.create(
            patient=patient,
            twin_state=state,
            predicted_values={"historical": 1},
            observed_values={"historical": 2},
            residuals={"historical": 1},
            normalized_residuals={"historical": 1},
            biomarker_weights={"historical": 1},
        )
        SimulationRunMetadata.objects.create(
            counterfactual_run=run,
            twin_state=state,
            model_version="historical-model",
            solver_name="historical-solver",
            input_hash="a" * 64,
        )

    def test_forward_marks_historical_rows_legacy_and_reverse_preserves_cardinality(self) -> None:
        executor = MigrationExecutor(connection)
        executor.migrate([self.migrate_to])
        apps = executor.loader.project_state([self.migrate_to]).apps
        PatientTwinState = apps.get_model("twin_engine", "PatientTwinState")
        CounterfactualRun = apps.get_model("twin_engine", "CounterfactualRun")
        ObservationResidual = apps.get_model("twin_engine", "ObservationResidual")
        SimulationRunMetadata = apps.get_model("twin_engine", "SimulationRunMetadata")
        ResearchRunManifest = apps.get_model("twin_engine", "ResearchRunManifest")

        state = PatientTwinState.objects.get()
        self.assertEqual(state.state_vector_schema_version, "legacy-unversioned-v0")
        self.assertEqual(state.parameters_schema_version, "legacy-unversioned-v0")
        self.assertEqual(state.parameter_uncertainty_schema_version, "legacy-unversioned-v0")
        self.assertEqual(state.lineage_schema_version, "legacy-unversioned-v0")
        run = CounterfactualRun.objects.get()
        self.assertEqual(run.intervention_schema_version, "legacy-unversioned-v0")
        self.assertEqual(run.simulation_summary_schema_version, "legacy-unversioned-v0")
        self.assertEqual(run.comparison_metrics_schema_version, "legacy-unversioned-v0")
        self.assertEqual(ObservationResidual.objects.get().payload_schema_version, "legacy-unversioned-v0")
        metadata = SimulationRunMetadata.objects.get()
        self.assertEqual(metadata.contract_version, "legacy-unversioned-v0")
        self.assertIsNone(metadata.manifest_id)
        self.assertEqual(ResearchRunManifest.objects.count(), 0)

        executor = MigrationExecutor(connection)
        executor.migrate([self.migrate_from])
        old_apps = executor.loader.project_state([self.migrate_from]).apps
        self.assertEqual(old_apps.get_model("twin_engine", "PatientTwinState").objects.count(), 1)
        self.assertEqual(old_apps.get_model("twin_engine", "CounterfactualRun").objects.count(), 1)
        self.assertEqual(old_apps.get_model("twin_engine", "ObservationResidual").objects.count(), 1)
        self.assertEqual(old_apps.get_model("twin_engine", "SimulationRunMetadata").objects.count(), 1)

    def _migrate_latest(self) -> None:
        executor = MigrationExecutor(connection)
        executor.migrate([self.migrate_to])

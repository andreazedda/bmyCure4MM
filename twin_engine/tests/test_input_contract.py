from __future__ import annotations

import json
from datetime import date
from io import StringIO

from django.core.management import call_command
from django.test import TestCase

from clinic.models import (
    Assessment,
    Patient,
    PatientTherapy,
    Regimen,
)
from twin_engine.input_contract import (
    COMPUTATIONAL_INPUT_CONTRACT,
    DATASET_BINDING_SCOPE,
    TWIN_LINEAGE_CONTRACT,
    build_computational_input_manifest,
    build_twin_lineage,
    computational_input_sha256,
)
from twin_engine.contracts import LEGACY_CONTRACT_VERSION
from twin_engine.models import (
    LongitudinalLabResult,
    ObservationResidual,
    PatientTwinState,
    TherapyInterruption,
)
from twin_engine.state_model import (
    initialize_from_assessment,
    set_current_state,
)


class ComputationalInputContractTests(
    TestCase
):
    def _create_case(
        self,
        suffix: str,
    ):
        patient = Patient.objects.create(
            mrn=f"SYN-COMP-{suffix}",
            first_name="Synthetic",
            last_name=f"Case{suffix}",
            birth_date=date(
                1970,
                1,
                1,
            ),
            sex="O",
            diagnosis_date=date(
                2032,
                1,
                1,
            ),
        )

        assessments = [
            Assessment.objects.create(
                patient=patient,
                date=date(
                    2032,
                    2,
                    1,
                ),
                m_protein_g_dl="1.20",
                flc_ratio="2.40",
                hemoglobin_g_dl="11.8",
                beta2m_mg_l="3.00",
                ldH_u_l="220.0",
                r_iss="II",
            ),
            Assessment.objects.create(
                patient=patient,
                date=date(
                    2032,
                    2,
                    15,
                ),
                m_protein_g_dl="1.10",
                flc_ratio="2.20",
                hemoglobin_g_dl="11.9",
                beta2m_mg_l="2.90",
                ldH_u_l="215.0",
                r_iss="II",
            ),
        ]

        regimen = Regimen.objects.create(
            name="SYNTHETIC_REGIMEN",
            line="synthetic",
            components="lenalidomide",
            intent="research_test",
        )

        therapy = (
            PatientTherapy.objects.create(
                patient=patient,
                regimen=regimen,
                start_date=date(
                    2032,
                    2,
                    1,
                ),
                end_date=date(
                    2032,
                    2,
                    15,
                ),
                doses={
                    "lenalidomide": {
                        "dose": 9.0,
                        "unit": "mg",
                        "schedule": {
                            "type": "daily",
                        },
                    }
                },
                cycle_length_days=1,
                days_on=[1],
                source_quality=(
                    PatientTherapy
                    .SOURCE_QUALITY_CURATED_RESEARCH
                ),
                provenance={
                    "source":
                        "synthetic_test_fixture",
                },
            )
        )

        interruption = (
            TherapyInterruption.objects.create(
                patient=patient,
                patient_therapy=therapy,
                start_date=date(
                    2032,
                    2,
                    10,
                ),
                end_date=date(
                    2032,
                    2,
                    11,
                ),
                drug="lenalidomide",
                reason=(
                    TherapyInterruption
                    .REASON_PLANNED_CYCLE_BREAK
                ),
                evidence={
                    "synthetic": True,
                },
            )
        )

        for assessment in assessments:
            for analyte, value, unit in (
                (
                    LongitudinalLabResult
                    .ANALYTE_M_PROTEIN,
                    float(
                        assessment
                        .m_protein_g_dl
                    ),
                    "g/dL",
                ),
                (
                    LongitudinalLabResult
                    .ANALYTE_FLC_RATIO,
                    float(
                        assessment
                        .flc_ratio
                    ),
                    "",
                ),
                (
                    LongitudinalLabResult
                    .ANALYTE_HB,
                    float(
                        assessment
                        .hemoglobin_g_dl
                    ),
                    "g/dL",
                ),
            ):
                (
                    LongitudinalLabResult
                    .objects
                    .create(
                        patient=patient,
                        assessment=assessment,
                        date=assessment.date,
                        analyte=analyte,
                        value=value,
                        unit=unit,
                    )
                )

        state = (
            initialize_from_assessment(
                assessments[0]
            )
        )

        return {
            "patient": patient,
            "assessments": assessments,
            "therapy": therapy,
            "interruption":
                interruption,
            "state": state,
        }

    def test_computational_hash_is_content_addressed_not_pk_addressed(
        self,
    ) -> None:
        first = self._create_case(
            "A"
        )
        second = self._create_case(
            "B"
        )

        def manifest(case):
            return (
                build_computational_input_manifest(
                    assessments=(
                        case["assessments"]
                    ),
                    therapies=[
                        case["therapy"]
                    ],
                    interruptions=[
                        case["interruption"]
                    ],
                    horizon_start_date=(
                        case["state"]
                        .state_date
                    ),
                    horizon_end_date=(
                        case["assessments"][-1]
                        .date
                    ),
                    purpose="calibration",
                    start_state=(
                        case["state"]
                    ),
                    extra={
                        "parameter_bounds": {
                            "synthetic": (
                                0.0,
                                1.0,
                            )
                        }
                    },
                )
            )

        first_manifest = manifest(first)
        second_manifest = manifest(
            second
        )

        self.assertEqual(
            first_manifest,
            second_manifest,
        )

        self.assertEqual(
            computational_input_sha256(
                first_manifest
            ),
            computational_input_sha256(
                second_manifest
            ),
        )

        rendered = json.dumps(
            first_manifest,
            sort_keys=True,
        )

        self.assertNotIn(
            '"id"',
            rendered,
        )

        self.assertEqual(
            first_manifest["contract"],
            COMPUTATIONAL_INPUT_CONTRACT,
        )

    def test_initialization_persists_hash_only_lineage_metadata(
        self,
    ) -> None:
        case = self._create_case(
            "INIT"
        )

        lineage = (
            case["state"].lineage
        )

        self.assertEqual(
            lineage["contract"],
            TWIN_LINEAGE_CONTRACT,
        )

        self.assertEqual(
            lineage["purpose"],
            "initialization",
        )

        computational = (
            lineage[
                "computational_input"
            ]
        )

        self.assertEqual(
            len(
                computational[
                    "sha256"
                ]
            ),
            64,
        )

        self.assertFalse(
            computational[
                "contains_local_database_ids"
            ]
        )

        self.assertFalse(
            computational[
                "raw_manifest_persisted"
            ]
        )

        self.assertEqual(
            lineage[
                "dataset_binding"
            ][
                "coverage_scope"
            ],
            DATASET_BINDING_SCOPE,
        )

        self.assertNotIn(
            "assessments",
            lineage,
        )

        self.assertNotIn(
            "therapies",
            lineage,
        )

    def test_historical_residuals_do_not_complete_current_lineage(
        self,
    ) -> None:
        case = self._create_case(
            "HIST"
        )

        baseline = case["state"]

        historical_parent = (
            PatientTwinState.objects
            .create(
                patient=case[
                    "patient"
                ],
                assessment=case[
                    "assessments"
                ][0],
                state_date=case[
                    "assessments"
                ][0].date,
                is_current=False,
                state_vector=(
                    baseline.state_vector
                ),
                parameters=(
                    baseline.parameters
                ),
                parameter_uncertainty={},
                risk_score=(
                    baseline.risk_score
                ),
                method=(
                    PatientTwinState
                    .METHOD_INITIAL_RISK_MAPPING
                ),
                model_version=(
                    baseline.model_version
                ),
                config_hash=(
                    baseline.config_hash
                ),
                state_vector_schema_version=LEGACY_CONTRACT_VERSION,
                parameters_schema_version=LEGACY_CONTRACT_VERSION,
                parameter_uncertainty_schema_version=LEGACY_CONTRACT_VERSION,
                lineage_schema_version=LEGACY_CONTRACT_VERSION,
            )
        )

        historical_calibrated = (
            PatientTwinState.objects
            .create(
                patient=case[
                    "patient"
                ],
                assessment=case[
                    "assessments"
                ][-1],
                state_date=case[
                    "assessments"
                ][-1].date,
                is_current=False,
                state_vector=(
                    baseline.state_vector
                ),
                parameters=(
                    baseline.parameters
                ),
                parameter_uncertainty={
                    "calibration_status":
                        "usable",
                    "rmse_before": 1.0,
                    "rmse_after": 0.5,
                },
                risk_score=(
                    baseline.risk_score
                ),
                method=(
                    PatientTwinState
                    .METHOD_RESIDUAL_MINIMIZATION
                ),
                model_version=(
                    baseline.model_version
                ),
                config_hash=(
                    baseline.config_hash
                ),
                state_vector_schema_version=LEGACY_CONTRACT_VERSION,
                parameters_schema_version=LEGACY_CONTRACT_VERSION,
                parameter_uncertainty_schema_version=LEGACY_CONTRACT_VERSION,
                lineage_schema_version=LEGACY_CONTRACT_VERSION,
            )
        )

        ObservationResidual.objects.create(
            patient=case["patient"],
            twin_state=historical_parent,
            assessment=case[
                "assessments"
            ][0],
            stage=(
                ObservationResidual
                .STAGE_PRE_CALIBRATION
            ),
        )

        ObservationResidual.objects.create(
            patient=case["patient"],
            twin_state=(
                historical_calibrated
            ),
            assessment=case[
                "assessments"
            ][-1],
            stage=(
                ObservationResidual
                .STAGE_POST_CALIBRATION
            ),
        )

        output = StringIO()

        call_command(
            "audit_patient_research_readiness",
            patient_id=case[
                "patient"
            ].pk,
            stdout=output,
        )

        payload = json.loads(
            output.getvalue()
        )

        self.assertTrue(
            payload[
                "calibration_input_ready"
            ]
        )

        self.assertEqual(
            payload[
                "historical_pre_calibration_residual_count"
            ],
            1,
        )

        self.assertEqual(
            payload[
                "historical_post_calibration_residual_count"
            ],
            1,
        )

        self.assertEqual(
            payload[
                "pre_calibration_residual_count"
            ],
            0,
        )

        self.assertEqual(
            payload[
                "post_calibration_residual_count"
            ],
            0,
        )

        self.assertFalse(
            payload[
                "calibration_lineage_ready"
            ]
        )

        self.assertFalse(
            payload[
                "calibration_ready"
            ]
        )

        self.assertFalse(
            payload[
                "calibrated_prediction_ready"
            ]
        )

    def test_parent_child_calibration_lineage_can_be_ready(
        self,
    ) -> None:
        case = self._create_case(
            "LINEAGE"
        )

        parent = case["state"]

        lineage = build_twin_lineage(
            patient=case[
                "patient"
            ],
            assessments=case[
                "assessments"
            ],
            therapies=[
                case["therapy"]
            ],
            purpose="calibration",
            parent_state=parent,
            extra={
                "parameter_bounds": {
                    "synthetic": (
                        0.0,
                        1.0,
                    )
                }
            },
        )

        calibrated = (
            PatientTwinState.objects
            .create(
                patient=case[
                    "patient"
                ],
                assessment=case[
                    "assessments"
                ][-1],
                state_date=case[
                    "assessments"
                ][-1].date,
                is_current=False,
                state_vector=(
                    parent.state_vector
                ),
                parameters=(
                    parent.parameters
                ),
                parameter_uncertainty={
                    "calibration_status":
                        "usable",
                    "rmse_before": 1.0,
                    "rmse_after": 0.5,
                },
                risk_score=(
                    parent.risk_score
                ),
                method=(
                    PatientTwinState
                    .METHOD_RESIDUAL_MINIMIZATION
                ),
                model_version=(
                    parent.model_version
                ),
                config_hash=(
                    parent.config_hash
                ),
                lineage=lineage,
            )
        )

        calibrated.source_assessments.add(
            *case["assessments"]
        )

        set_current_state(
            calibrated
        )

        ObservationResidual.objects.create(
            patient=case["patient"],
            twin_state=parent,
            assessment=case[
                "assessments"
            ][0],
            stage=(
                ObservationResidual
                .STAGE_PRE_CALIBRATION
            ),
        )

        ObservationResidual.objects.create(
            patient=case["patient"],
            twin_state=calibrated,
            assessment=case[
                "assessments"
            ][-1],
            stage=(
                ObservationResidual
                .STAGE_POST_CALIBRATION
            ),
        )

        output = StringIO()

        call_command(
            "audit_patient_research_readiness",
            patient_id=case[
                "patient"
            ].pk,
            stdout=output,
        )

        payload = json.loads(
            output.getvalue()
        )

        self.assertEqual(
            payload[
                "calibration_parent_state_id"
            ],
            parent.pk,
        )

        self.assertEqual(
            payload[
                "pre_calibration_residual_count"
            ],
            1,
        )

        self.assertEqual(
            payload[
                "post_calibration_residual_count"
            ],
            1,
        )

        self.assertTrue(
            payload[
                "calibration_lineage_ready"
            ]
        )

        self.assertTrue(
            payload[
                "calibration_ready"
            ]
        )

        self.assertTrue(
            payload[
                "calibrated_prediction_ready"
            ]
        )

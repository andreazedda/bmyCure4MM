from __future__ import annotations

import json

from django.core.management.base import (
    BaseCommand,
    CommandError,
)

from clinic.models import Patient
from twin_engine.models import (
    AdverseEvent,
    LongitudinalLabResult,
    ObservationResidual,
    PatientTwinState,
    TherapyInterruption,
)
from twin_engine.state_model import (
    get_current_twin_state,
)
from twin_engine.therapy_schedule import (
    build_therapy_schedule,
)
from twin_engine.toxicity_model import (
    compute_toxicity_constraints,
)
from twin_engine.validators import (
    BANNED_IDENTIFIER_KEYS,
    validate_assessment_minimum_fields,
)


class Command(BaseCommand):
    help = (
        "Audit a patient record for research-twin "
        "readiness and privacy risks."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--patient-id",
            type=int,
            required=True,
        )

    def handle(self, *args, **options):
        patient = (
            Patient.objects
            .prefetch_related(
                "assessments",
                "therapies",
            )
            .filter(
                pk=options["patient_id"]
            )
            .first()
        )

        if patient is None:
            raise CommandError(
                "Patient not found"
            )

        assessments = list(
            patient.assessments.order_by(
                "date"
            )
        )

        therapies = list(
            patient.therapies.order_by(
                "start_date"
            )
        )

        latest_assessment = (
            assessments[-1]
            if assessments
            else None
        )

        current_state = (
            get_current_twin_state(
                patient
            )
        )

        missing_fields = []

        if latest_assessment is not None:
            try:
                missing_fields = (
                    validate_assessment_minimum_fields(
                        latest_assessment
                    )["missing"]
                )
            except Exception as exc:
                missing_fields = [
                    str(exc)
                ]
        else:
            missing_fields = [
                "assessment"
            ]

        if therapies and assessments:
            schedule = (
                build_therapy_schedule(
                    patient,
                    therapies[0].start_date,
                    assessments[-1].date,
                )
            )

            schedule_validation = (
                schedule.get(
                    "validation",
                    {},
                )
            )
        else:
            schedule_validation = {
                "is_valid": False,
                "missing_doses": [
                    {
                        "reason":
                            "no therapy timeline"
                    }
                ],
            }

        privacy_risks = []

        if (
            patient.first_name
            or patient.last_name
            or patient.mrn
        ):
            privacy_risks.append(
                "direct identifiers exist in "
                "clinic.Patient and must not be "
                "exported into research artifacts"
            )

        if patient.notes:
            privacy_risks.append(
                "patient notes require exclusion "
                "from research artifacts by default"
            )

        for banned_key in (
            "notes",
            "schedule_notes",
        ):
            if banned_key in BANNED_IDENTIFIER_KEYS:
                privacy_risks.append(
                    "artifact validator blocks "
                    f"{banned_key}"
                )

        disease_marker_observation_count = (
            LongitudinalLabResult.objects
            .filter(
                patient=patient,
                analyte__in=[
                    LongitudinalLabResult
                    .ANALYTE_M_PROTEIN,
                    LongitudinalLabResult
                    .ANALYTE_FLC_RATIO,
                    LongitudinalLabResult
                    .ANALYTE_HB,
                ],
            )
            .exclude(
                value__isnull=True
            )
            .count()
        )

        if (
            disease_marker_observation_count
            == 0
        ):
            disease_marker_observation_count = sum(
                1
                for assessment in assessments
                if any(
                    getattr(
                        assessment,
                        field_name,
                        None,
                    )
                    not in {
                        None,
                        "",
                    }
                    for field_name in (
                        "m_protein_g_dl",
                        "flc_ratio",
                        "hemoglobin_g_dl",
                    )
                )
            )

        toxicity_lab_count = (
            LongitudinalLabResult.objects
            .filter(
                patient=patient,
                analyte__in=[
                    LongitudinalLabResult
                    .ANALYTE_AST,
                    LongitudinalLabResult
                    .ANALYTE_ALT,
                    LongitudinalLabResult
                    .ANALYTE_NEU,
                ],
            )
            .exclude(
                value__isnull=True
            )
            .count()
        )

        adverse_event_count = (
            AdverseEvent.objects.filter(
                patient=patient
            ).count()
        )

        interruption_count = (
            TherapyInterruption.objects
            .filter(
                patient=patient
            )
            .count()
        )

        historical_pre_count = (
            ObservationResidual.objects
            .filter(
                patient=patient,
                stage=(
                    ObservationResidual
                    .STAGE_PRE_CALIBRATION
                ),
            )
            .count()
        )

        historical_post_count = (
            ObservationResidual.objects
            .filter(
                patient=patient,
                stage=(
                    ObservationResidual
                    .STAGE_POST_CALIBRATION
                ),
            )
            .count()
        )

        parent_state = None
        parent_state_id = None

        if current_state is not None:
            lineage = dict(
                current_state.lineage
                or {}
            )

            parent_payload = dict(
                lineage.get(
                    "parent",
                    {},
                )
                or {}
            )

            raw_parent_id = (
                parent_payload.get(
                    "local_twin_state_id"
                )
            )

            if raw_parent_id is not None:
                try:
                    parent_state_id = int(
                        raw_parent_id
                    )
                except (
                    TypeError,
                    ValueError,
                ):
                    parent_state_id = None

            if parent_state_id is not None:
                parent_state = (
                    PatientTwinState.objects
                    .filter(
                        patient=patient,
                        pk=parent_state_id,
                    )
                    .first()
                )

        pre_calibration_residual_count = 0
        post_calibration_residual_count = 0

        if (
            current_state is not None
            and current_state.method
            == PatientTwinState
            .METHOD_RESIDUAL_MINIMIZATION
            and parent_state is not None
        ):
            pre_calibration_residual_count = (
                ObservationResidual.objects
                .filter(
                    patient=patient,
                    twin_state=parent_state,
                    stage=(
                        ObservationResidual
                        .STAGE_PRE_CALIBRATION
                    ),
                )
                .count()
            )

            post_calibration_residual_count = (
                ObservationResidual.objects
                .filter(
                    patient=patient,
                    twin_state=current_state,
                    stage=(
                        ObservationResidual
                        .STAGE_POST_CALIBRATION
                    ),
                )
                .count()
            )

        toxicity_constraints = (
            compute_toxicity_constraints(
                patient
            )
        )

        disease_marker_whatif_ready = bool(
            current_state
            and disease_marker_observation_count
            >= 3
            and therapies
        )

        toxicity_constrained_whatif_ready = bool(
            disease_marker_whatif_ready
            and (
                toxicity_lab_count > 0
                or adverse_event_count > 0
            )
            and interruption_count > 0
        )

        calibration_input_ready = bool(
            current_state
            and len(assessments) >= 2
            and disease_marker_observation_count
            >= 6
            and therapies
            and schedule_validation.get(
                "is_valid"
            )
        )

        calibration_lineage_ready = bool(
            current_state
            and current_state.method
            == PatientTwinState
            .METHOD_RESIDUAL_MINIMIZATION
            and parent_state is not None
            and pre_calibration_residual_count
            > 0
            and post_calibration_residual_count
            > 0
        )

        calibration_status = dict(
            (
                current_state
                .parameter_uncertainty
                or {}
            )
            if current_state is not None
            else {}
        )

        calibrated_prediction_ready = bool(
            calibration_lineage_ready
            and calibration_status.get(
                "calibration_status"
            )
            == "usable"
            and calibration_status.get(
                "rmse_after"
            )
            is not None
            and calibration_status.get(
                "rmse_before"
            )
            is not None
            and float(
                calibration_status[
                    "rmse_after"
                ]
            )
            < float(
                calibration_status[
                    "rmse_before"
                ]
            )
        )

        readiness = {
            "patient_id": patient.id,
            "missing_fields":
                missing_fields,
            "assessment_count":
                len(assessments),
            "therapy_schedule_completeness":
                schedule_validation,
            "dose_completeness":
                schedule_validation.get(
                    "missing_doses",
                    [],
                ),
            "calibratability_status":
                calibration_input_ready,
            "counterfactual_readiness_status":
                disease_marker_whatif_ready,
            "disease_marker_observation_count":
                disease_marker_observation_count,
            # These are now CURRENT CALIBRATION
            # LINEAGE counts, not patient-global
            # historical counts.
            "pre_calibration_residual_count":
                pre_calibration_residual_count,
            "post_calibration_residual_count":
                post_calibration_residual_count,
            "historical_pre_calibration_residual_count":
                historical_pre_count,
            "historical_post_calibration_residual_count":
                historical_post_count,
            "calibration_parent_state_id":
                (
                    parent_state.pk
                    if parent_state
                    is not None
                    else None
                ),
            "adverse_event_count":
                adverse_event_count,
            "therapy_interruption_count":
                interruption_count,
            "disease_marker_whatif_ready":
                disease_marker_whatif_ready,
            "toxicity_constrained_whatif_ready":
                toxicity_constrained_whatif_ready,
            "calibration_input_ready":
                calibration_input_ready,
            "calibration_lineage_ready":
                calibration_lineage_ready,
            # Compatibility field: this now means
            # that the CURRENT calibration lineage
            # has both pre and post evidence.
            "calibration_ready":
                calibration_lineage_ready,
            "calibrated_prediction_ready":
                calibrated_prediction_ready,
            "causal_inference_ready":
                False,
            "toxicity_constraints":
                toxicity_constraints,
            "privacy_risks":
                privacy_risks,
        }

        self.stdout.write(
            json.dumps(
                readiness,
                indent=2,
            )
        )

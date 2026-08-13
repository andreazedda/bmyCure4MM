from __future__ import annotations

from django.conf import settings
from django.core.exceptions import ValidationError
from django.db import models

from .contracts import (
    COMPARISON_METRICS_SCHEMA_VERSION,
    INTERVENTION_SCHEMA_VERSION,
    LINEAGE_SCHEMA_VERSION,
    OBSERVATION_RESIDUAL_SCHEMA_VERSION,
    PARAMETERS_SCHEMA_VERSION,
    PARAMETER_UNCERTAINTY_SCHEMA_VERSION,
    SIMULATION_SUMMARY_SCHEMA_VERSION,
    STATE_VECTOR_SCHEMA_VERSION,
    validate_counterfactual_contracts,
    validate_residual_contract,
    validate_twin_contracts,
)


class PatientTwinState(models.Model):
    METHOD_INITIAL_RISK_MAPPING = "initial_risk_mapping"
    METHOD_RESIDUAL_MINIMIZATION = "residual_minimization"
    METHOD_MANUAL_RESEARCH_OVERRIDE = "manual_research_override"
    METHOD_BAYESIAN_UPDATE_PLACEHOLDER = "bayesian_update_placeholder"

    METHOD_CHOICES = [
        (METHOD_INITIAL_RISK_MAPPING, "Initial risk mapping"),
        (METHOD_RESIDUAL_MINIMIZATION, "Residual minimization"),
        (METHOD_MANUAL_RESEARCH_OVERRIDE, "Manual research override"),
        (METHOD_BAYESIAN_UPDATE_PLACEHOLDER, "Bayesian update placeholder"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="twin_states",
    )
    assessment = models.ForeignKey(
        "clinic.Assessment",
        on_delete=models.SET_NULL,
        related_name="primary_twin_states",
        null=True,
        blank=True,
    )
    state_date = models.DateField()
    is_current = models.BooleanField(default=False)
    state_vector = models.JSONField(default=dict, blank=True)
    state_vector_schema_version = models.CharField(
        max_length=64,
        default=STATE_VECTOR_SCHEMA_VERSION,
    )
    parameters = models.JSONField(default=dict, blank=True)
    parameters_schema_version = models.CharField(
        max_length=64,
        default=PARAMETERS_SCHEMA_VERSION,
    )
    parameter_uncertainty = models.JSONField(default=dict, blank=True)
    parameter_uncertainty_schema_version = models.CharField(
        max_length=64,
        default=PARAMETER_UNCERTAINTY_SCHEMA_VERSION,
    )
    risk_score = models.FloatField(null=True, blank=True)
    method = models.CharField(
        max_length=64,
        choices=METHOD_CHOICES,
        default=METHOD_INITIAL_RISK_MAPPING,
    )
    model_version = models.CharField(max_length=64)
    config_hash = models.CharField(max_length=128)
    lineage = models.JSONField(default=dict, blank=True)
    lineage_schema_version = models.CharField(
        max_length=64,
        default=LINEAGE_SCHEMA_VERSION,
    )
    source_assessments = models.ManyToManyField(
        "clinic.Assessment",
        related_name="twin_state_sources",
        blank=True,
    )
    created_at = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.SET_NULL,
        related_name="created_twin_states",
        null=True,
        blank=True,
    )

    class Meta:
        ordering = ["-state_date", "-created_at"]
        constraints = [
            models.UniqueConstraint(
                fields=["patient"],
                condition=models.Q(is_current=True),
                name="unique_current_twin_state_per_patient",
            )
        ]

    def __str__(self) -> str:
        return f"Twin state patient={self.patient_id} date={self.state_date}"

    def validate_scientific_contracts(self) -> None:
        validate_twin_contracts(
            state_vector=self.state_vector,
            state_vector_schema_version=self.state_vector_schema_version,
            parameters=self.parameters,
            parameters_schema_version=self.parameters_schema_version,
            parameter_uncertainty=self.parameter_uncertainty,
            parameter_uncertainty_schema_version=self.parameter_uncertainty_schema_version,
            lineage=self.lineage,
            lineage_schema_version=self.lineage_schema_version,
        )

    def save(self, *args, **kwargs) -> None:
        self.validate_scientific_contracts()
        super().save(*args, **kwargs)


class ObservationResidual(models.Model):
    STAGE_PRE_CALIBRATION = "pre_calibration"
    STAGE_POST_CALIBRATION = "post_calibration"
    STAGE_SIMULATION_VALIDATION = "simulation_validation"
    STAGE_UNKNOWN = "unknown"

    STAGE_CHOICES = [
        (STAGE_PRE_CALIBRATION, "Pre calibration"),
        (STAGE_POST_CALIBRATION, "Post calibration"),
        (STAGE_SIMULATION_VALIDATION, "Simulation validation"),
        (STAGE_UNKNOWN, "Unknown"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="observation_residuals",
    )
    twin_state = models.ForeignKey(
        PatientTwinState,
        on_delete=models.CASCADE,
        related_name="residuals",
    )
    assessment = models.ForeignKey(
        "clinic.Assessment",
        on_delete=models.SET_NULL,
        related_name="residual_records",
        null=True,
        blank=True,
    )
    predicted_values = models.JSONField(default=dict, blank=True)
    observed_values = models.JSONField(default=dict, blank=True)
    residuals = models.JSONField(default=dict, blank=True)
    normalized_residuals = models.JSONField(default=dict, blank=True)
    rmse = models.FloatField(null=True, blank=True)
    mae = models.FloatField(null=True, blank=True)
    biomarker_weights = models.JSONField(default=dict, blank=True)
    payload_schema_version = models.CharField(
        max_length=64,
        default=OBSERVATION_RESIDUAL_SCHEMA_VERSION,
    )
    stage = models.CharField(max_length=32, choices=STAGE_CHOICES, default=STAGE_UNKNOWN)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self) -> str:
        return f"Observation residual patient={self.patient_id} twin_state={self.twin_state_id} stage={self.stage}"

    def save(self, *args, **kwargs) -> None:
        validate_residual_contract(self)
        super().save(*args, **kwargs)


class LongitudinalLabResult(models.Model):
    ANALYTE_HB = "HB"
    ANALYTE_WBC = "WBC"
    ANALYTE_NEU = "NEU"
    ANALYTE_PLT = "PLT"
    ANALYTE_CREATININE = "CREATININE"
    ANALYTE_AST = "AST"
    ANALYTE_ALT = "ALT"
    ANALYTE_M_PROTEIN = "M_PROTEIN"
    ANALYTE_KAPPA_FLC = "KAPPA_FLC"
    ANALYTE_LAMBDA_FLC = "LAMBDA_FLC"
    ANALYTE_FLC_RATIO = "FLC_RATIO"
    ANALYTE_LDH = "LDH"
    ANALYTE_BETA2M = "BETA2M"
    ANALYTE_ALBUMIN = "ALBUMIN"
    ANALYTE_CALCIUM = "CALCIUM"
    ANALYTE_OTHER = "OTHER"

    ANALYTE_CHOICES = [
        (ANALYTE_HB, "Hemoglobin"),
        (ANALYTE_WBC, "White blood cells"),
        (ANALYTE_NEU, "Neutrophils"),
        (ANALYTE_PLT, "Platelets"),
        (ANALYTE_CREATININE, "Creatinine"),
        (ANALYTE_AST, "AST"),
        (ANALYTE_ALT, "ALT"),
        (ANALYTE_M_PROTEIN, "M-protein"),
        (ANALYTE_KAPPA_FLC, "Kappa free light chains"),
        (ANALYTE_LAMBDA_FLC, "Lambda free light chains"),
        (ANALYTE_FLC_RATIO, "FLC ratio"),
        (ANALYTE_LDH, "LDH"),
        (ANALYTE_BETA2M, "Beta-2 microglobulin"),
        (ANALYTE_ALBUMIN, "Albumin"),
        (ANALYTE_CALCIUM, "Calcium"),
        (ANALYTE_OTHER, "Other"),
    ]

    SOURCE_QUALITY_UNKNOWN = "unknown"
    SOURCE_QUALITY_EXTRACTED_DOCUMENT = "extracted_document"
    SOURCE_QUALITY_CLINICAL_RECORD = "clinical_record"
    SOURCE_QUALITY_INFERRED_FROM_ASSESSMENT = "inferred_from_assessment"
    SOURCE_QUALITY_MANUALLY_CURATED = "manually_curated"

    SOURCE_QUALITY_CHOICES = [
        (SOURCE_QUALITY_UNKNOWN, "Unknown"),
        (SOURCE_QUALITY_EXTRACTED_DOCUMENT, "Extracted document"),
        (SOURCE_QUALITY_CLINICAL_RECORD, "Clinical record"),
        (SOURCE_QUALITY_INFERRED_FROM_ASSESSMENT, "Inferred from assessment"),
        (SOURCE_QUALITY_MANUALLY_CURATED, "Manually curated"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="longitudinal_lab_results",
    )
    assessment = models.ForeignKey(
        "clinic.Assessment",
        on_delete=models.SET_NULL,
        related_name="longitudinal_lab_results",
        null=True,
        blank=True,
    )
    date = models.DateField()
    analyte = models.CharField(max_length=32, choices=ANALYTE_CHOICES)
    value = models.FloatField(null=True, blank=True)
    unit = models.CharField(max_length=32, blank=True)
    source_quality = models.CharField(
        max_length=32,
        choices=SOURCE_QUALITY_CHOICES,
        default=SOURCE_QUALITY_UNKNOWN,
    )
    provenance = models.JSONField(default=dict, blank=True)
    notes = models.TextField(blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-date", "-created_at"]
        constraints = [
            models.UniqueConstraint(
                fields=["patient", "date", "analyte"],
                name="unique_longitudinal_lab_per_patient_date_analyte",
            )
        ]
        indexes = [
            models.Index(fields=["patient", "date"], name="lab_patient_date_idx"),
            models.Index(fields=["patient", "analyte", "date"], name="lab_pat_an_dt_idx"),
        ]

    def __str__(self) -> str:
        return f"Lab result patient={self.patient_id} {self.analyte} on {self.date}"


class TherapyInterruption(models.Model):
    REASON_HYPERTRANSAMINASEMIA = "hypertransaminasemia"
    REASON_NEUTROPENIA = "neutropenia"
    REASON_INFECTION = "infection"
    REASON_DIARRHEA = "diarrhea"
    REASON_TOXICITY_UNSPECIFIED = "toxicity_unspecified"
    REASON_PLANNED_CYCLE_BREAK = "planned_cycle_break"
    REASON_UNKNOWN = "unknown"

    REASON_CHOICES = [
        (REASON_HYPERTRANSAMINASEMIA, "Hypertransaminasemia"),
        (REASON_NEUTROPENIA, "Neutropenia"),
        (REASON_INFECTION, "Infection"),
        (REASON_DIARRHEA, "Diarrhea"),
        (REASON_TOXICITY_UNSPECIFIED, "Toxicity unspecified"),
        (REASON_PLANNED_CYCLE_BREAK, "Planned cycle break"),
        (REASON_UNKNOWN, "Unknown"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="therapy_interruptions",
    )
    patient_therapy = models.ForeignKey(
        "clinic.PatientTherapy",
        on_delete=models.SET_NULL,
        related_name="interruptions",
        null=True,
        blank=True,
    )
    start_date = models.DateField()
    end_date = models.DateField(null=True, blank=True)
    drug = models.CharField(max_length=64)
    reason = models.CharField(max_length=64, choices=REASON_CHOICES, default=REASON_UNKNOWN)
    evidence = models.JSONField(default=dict, blank=True)
    action_taken = models.TextField(blank=True)
    source_quality = models.CharField(max_length=32, default="clinical_record")
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-start_date", "-created_at"]
        constraints = [
            models.UniqueConstraint(
                fields=["patient", "drug", "start_date", "end_date"],
                name="unique_interruption_per_patient_drug_window",
            )
        ]
        indexes = [
            models.Index(fields=["patient", "start_date"], name="intr_pat_dt_idx"),
            models.Index(fields=["patient", "drug", "start_date"], name="intr_pat_drug_dt_idx"),
        ]

    def __str__(self) -> str:
        return f"Therapy interruption patient={self.patient_id} drug={self.drug} start={self.start_date}"


class AdverseEvent(models.Model):
    TYPE_HYPERTRANSAMINASEMIA = "hypertransaminasemia"
    TYPE_NEUTROPENIA = "neutropenia"
    TYPE_PNEUMONIA = "pneumonia"
    TYPE_UPPER_RESPIRATORY_INFECTION = "upper_respiratory_infection"
    TYPE_DIARRHEA = "diarrhea"
    TYPE_NEUROPATHY = "neuropathy"
    TYPE_PAIN = "pain"
    TYPE_HEPATIC_STEATOSIS = "hepatic_steatosis"
    TYPE_OTHER = "other"

    EVENT_TYPE_CHOICES = [
        (TYPE_HYPERTRANSAMINASEMIA, "Hypertransaminasemia"),
        (TYPE_NEUTROPENIA, "Neutropenia"),
        (TYPE_PNEUMONIA, "Pneumonia"),
        (TYPE_UPPER_RESPIRATORY_INFECTION, "Upper respiratory infection"),
        (TYPE_DIARRHEA, "Diarrhea"),
        (TYPE_NEUROPATHY, "Neuropathy"),
        (TYPE_PAIN, "Pain"),
        (TYPE_HEPATIC_STEATOSIS, "Hepatic steatosis"),
        (TYPE_OTHER, "Other"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="adverse_events",
    )
    date = models.DateField()
    event_type = models.CharField(max_length=64, choices=EVENT_TYPE_CHOICES)
    grade = models.CharField(max_length=32, blank=True)
    suspected_drug = models.CharField(max_length=128, blank=True)
    observed_values = models.JSONField(default=dict, blank=True)
    action_taken = models.TextField(blank=True)
    outcome = models.TextField(blank=True)
    provenance = models.JSONField(default=dict, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-date", "-created_at"]
        constraints = [
            models.UniqueConstraint(
                fields=["patient", "date", "event_type"],
                name="unique_adverse_event_per_patient_date_type",
            )
        ]
        indexes = [
            models.Index(fields=["patient", "date"], name="event_patient_date_idx"),
            models.Index(fields=["patient", "event_type", "date"], name="evt_pat_type_dt_idx"),
        ]

    def __str__(self) -> str:
        return f"Adverse event patient={self.patient_id} type={self.event_type} date={self.date}"


class CounterfactualRun(models.Model):
    STATUS_DRAFT = "draft"
    STATUS_RUNNING = "running"
    STATUS_COMPLETED = "completed"
    STATUS_FAILED = "failed"

    STATUS_CHOICES = [
        (STATUS_DRAFT, "Draft"),
        (STATUS_RUNNING, "Running"),
        (STATUS_COMPLETED, "Completed"),
        (STATUS_FAILED, "Failed"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="counterfactual_runs",
    )
    base_twin_state = models.ForeignKey(
        PatientTwinState,
        on_delete=models.CASCADE,
        related_name="counterfactual_runs",
    )
    actual_therapy = models.ForeignKey(
        "clinic.PatientTherapy",
        on_delete=models.SET_NULL,
        related_name="counterfactual_actual_runs",
        null=True,
        blank=True,
    )
    alternative_regimen = models.ForeignKey(
        "clinic.Regimen",
        on_delete=models.SET_NULL,
        related_name="counterfactual_alternatives",
        null=True,
        blank=True,
    )
    alternative_parameters = models.JSONField(default=dict, blank=True)
    intervention_definition = models.JSONField(default=dict, blank=True)
    intervention_schema_version = models.CharField(
        max_length=64,
        default=INTERVENTION_SCHEMA_VERSION,
    )
    simulation_summary = models.JSONField(default=dict, blank=True)
    simulation_summary_schema_version = models.CharField(
        max_length=64,
        default=SIMULATION_SUMMARY_SCHEMA_VERSION,
    )
    comparison_metrics = models.JSONField(default=dict, blank=True)
    comparison_metrics_schema_version = models.CharField(
        max_length=64,
        default=COMPARISON_METRICS_SCHEMA_VERSION,
    )
    trajectory_artifact = models.TextField(blank=True)
    report_artifact = models.TextField(blank=True)
    status = models.CharField(max_length=16, choices=STATUS_CHOICES, default=STATUS_DRAFT)
    error_message = models.TextField(blank=True)
    created_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.SET_NULL,
        related_name="created_counterfactual_runs",
        null=True,
        blank=True,
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self) -> str:
        return f"Counterfactual run patient={self.patient_id} state={self.base_twin_state_id}"

    def _assert_not_completed_mutation(self) -> None:
        if not self.pk:
            return
        previous = type(self).objects.filter(pk=self.pk).first()
        if previous is None or previous.status != self.STATUS_COMPLETED:
            return
        immutable_fields = (
            "patient_id",
            "base_twin_state_id",
            "actual_therapy_id",
            "alternative_regimen_id",
            "alternative_parameters",
            "intervention_definition",
            "intervention_schema_version",
            "simulation_summary",
            "simulation_summary_schema_version",
            "comparison_metrics",
            "comparison_metrics_schema_version",
            "trajectory_artifact",
            "report_artifact",
            "status",
            "error_message",
            "created_by_id",
        )
        changed = [name for name in immutable_fields if getattr(previous, name) != getattr(self, name)]
        if changed:
            raise ValidationError(
                "Completed counterfactual runs are immutable; create a new run. "
                f"Changed fields: {', '.join(changed)}"
            )

    def save(self, *args, **kwargs) -> None:
        validate_counterfactual_contracts(self)
        self._assert_not_completed_mutation()
        super().save(*args, **kwargs)

    def delete(self, *args, **kwargs):
        if self.status == self.STATUS_COMPLETED:
            raise ValidationError("Completed counterfactual runs cannot be deleted")
        return super().delete(*args, **kwargs)


class CausalAssumptionSet(models.Model):
    IDENT_NOT_IDENTIFIED = "not_identified"
    IDENT_PARTIALLY_IDENTIFIED = "partially_identified"
    IDENT_IDENTIFIED_UNDER_ASSUMPTIONS = "identified_under_assumptions"
    IDENT_MECHANISTIC_ONLY = "mechanistic_only"

    IDENTIFICATION_CHOICES = [
        (IDENT_NOT_IDENTIFIED, "Not identified"),
        (IDENT_PARTIALLY_IDENTIFIED, "Partially identified"),
        (IDENT_IDENTIFIED_UNDER_ASSUMPTIONS, "Identified under assumptions"),
        (IDENT_MECHANISTIC_ONLY, "Mechanistic only"),
    ]

    patient = models.ForeignKey(
        "clinic.Patient",
        on_delete=models.CASCADE,
        related_name="causal_assumption_sets",
        null=True,
        blank=True,
    )
    name = models.CharField(max_length=255)
    graph_definition = models.JSONField(default=dict, blank=True)
    variables = models.JSONField(default=list, blank=True)
    intervention = models.JSONField(default=dict, blank=True)
    outcome = models.JSONField(default=dict, blank=True)
    adjustment_set = models.JSONField(default=list, blank=True)
    assumptions = models.JSONField(default=dict, blank=True)
    identification_status = models.CharField(
        max_length=48,
        choices=IDENTIFICATION_CHOICES,
        default=IDENT_MECHANISTIC_ONLY,
    )
    notes = models.TextField(blank=True)
    created_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.SET_NULL,
        related_name="created_causal_assumption_sets",
        null=True,
        blank=True,
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self) -> str:
        return self.name


class ResearchRunManifest(models.Model):
    STATUS_COMPLETED = "completed"
    STATUS_FAILED = "failed"
    STATUS_CHOICES = [
        (STATUS_COMPLETED, "Completed"),
        (STATUS_FAILED, "Failed"),
    ]

    run_id = models.CharField(max_length=192, unique=True)
    contract_version = models.CharField(max_length=64, default="research-run-manifest-v1")
    status = models.CharField(max_length=16, choices=STATUS_CHOICES)
    app_version = models.CharField(max_length=64)
    api_version = models.CharField(max_length=64)
    db_schema_version = models.CharField(max_length=128)
    dataset_id = models.CharField(max_length=255)
    dataset_version = models.CharField(max_length=128)
    dataset_sha256 = models.CharField(max_length=128)
    record_subset_sha256 = models.CharField(max_length=128)
    model_id = models.CharField(max_length=128)
    model_version = models.CharField(max_length=128)
    model_card_version = models.CharField(max_length=128)
    configuration_sha256 = models.CharField(max_length=128)
    evidence_graph_version = models.CharField(max_length=128)
    validation_protocol_version = models.CharField(max_length=128)
    report_template_version = models.CharField(max_length=128)
    git_sha = models.CharField(max_length=128)
    container_digest = models.CharField(max_length=255)
    dependency_lock_sha256 = models.CharField(max_length=128)
    model_registry_sha256 = models.CharField(max_length=128)
    random_seed = models.IntegerField()
    intended_use_level = models.CharField(max_length=128)
    epistemic_label = models.CharField(max_length=64)
    artifact_manifest = models.JSONField(default=list)
    manifest_artifact_uri = models.TextField()
    manifest_sha256 = models.CharField(max_length=128)
    created_at = models.DateTimeField()

    class Meta:
        ordering = ["-created_at"]

    def __str__(self) -> str:
        return self.run_id

    def version_vector(self) -> dict[str, object]:
        fields = (
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
            "dependency_lock_sha256",
            "model_registry_sha256",
            "random_seed",
            "intended_use_level",
            "epistemic_label",
        )
        payload = {name: getattr(self, name) for name in fields}
        payload["created_at"] = self.created_at.isoformat()
        return payload

    def save(self, *args, **kwargs) -> None:
        if self.pk:
            raise ValidationError("Research run manifests are immutable")
        from .run_manifest import validate_manifest_instance

        validate_manifest_instance(self)
        super().save(*args, **kwargs)

    def delete(self, *args, **kwargs):
        raise ValidationError("Research run manifests cannot be deleted")


class ResearchRunInvalidation(models.Model):
    CHANGE_CLINICAL_VALUE = "clinical_value"
    CHANGE_DATASET = "dataset"
    CHANGE_UNIT_MAPPING = "unit_mapping"
    CHANGE_MODEL = "model"
    CHANGE_PARAMETER_DEFAULT = "parameter_default"
    CHANGE_DEPENDENCY = "dependency"
    CHANGE_CONFIGURATION = "configuration"
    CHANGE_CHOICES = [
        (CHANGE_CLINICAL_VALUE, "Clinical value"),
        (CHANGE_DATASET, "Dataset"),
        (CHANGE_UNIT_MAPPING, "Unit mapping"),
        (CHANGE_MODEL, "Model"),
        (CHANGE_PARAMETER_DEFAULT, "Parameter default"),
        (CHANGE_DEPENDENCY, "Dependency"),
        (CHANGE_CONFIGURATION, "Configuration"),
    ]

    manifest = models.ForeignKey(
        ResearchRunManifest,
        on_delete=models.PROTECT,
        related_name="invalidations",
    )
    change_kind = models.CharField(max_length=32, choices=CHANGE_CHOICES)
    previous_identity = models.CharField(max_length=255)
    replacement_identity = models.CharField(max_length=255)
    change_sha256 = models.CharField(max_length=128)
    reason = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["created_at"]
        constraints = [
            models.UniqueConstraint(
                fields=["manifest", "change_kind", "change_sha256"],
                name="unique_run_invalidation_change",
            )
        ]

    def save(self, *args, **kwargs) -> None:
        if self.pk:
            raise ValidationError("Research run invalidations are immutable")
        super().save(*args, **kwargs)

    def delete(self, *args, **kwargs):
        raise ValidationError("Research run invalidations cannot be deleted")


class SimulationRunMetadata(models.Model):
    contract_version = models.CharField(max_length=64, default="simulation-metadata-v1")
    manifest = models.OneToOneField(
        ResearchRunManifest,
        on_delete=models.PROTECT,
        related_name="simulation_metadata",
        null=True,
        blank=True,
    )
    simulation_attempt = models.ForeignKey(
        "simulator.SimulationAttempt",
        on_delete=models.SET_NULL,
        related_name="research_metadata",
        null=True,
        blank=True,
    )
    counterfactual_run = models.ForeignKey(
        CounterfactualRun,
        on_delete=models.SET_NULL,
        related_name="metadata_records",
        null=True,
        blank=True,
    )
    twin_state = models.ForeignKey(
        PatientTwinState,
        on_delete=models.SET_NULL,
        related_name="metadata_records",
        null=True,
        blank=True,
    )
    model_version = models.CharField(max_length=64)
    code_commit_hash = models.CharField(max_length=128, blank=True)
    twin_config_hash = models.CharField(max_length=128, blank=True)
    drug_preset_hashes = models.JSONField(default=dict, blank=True)
    solver_name = models.CharField(max_length=128)
    solver_parameters = models.JSONField(default=dict, blank=True)
    input_hash = models.CharField(max_length=128)
    output_hash = models.CharField(max_length=128, blank=True)
    random_seed = models.IntegerField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self) -> str:
        return f"Simulation metadata model={self.model_version} created={self.created_at:%Y-%m-%d}"

    def save(self, *args, **kwargs) -> None:
        if self.pk:
            previous = type(self).objects.filter(pk=self.pk).values_list("manifest_id", flat=True).first()
            if previous is not None:
                raise ValidationError("Manifest-bound simulation metadata is immutable")
        super().save(*args, **kwargs)

    def delete(self, *args, **kwargs):
        if self.manifest_id is not None:
            raise ValidationError("Manifest-bound simulation metadata cannot be deleted")
        return super().delete(*args, **kwargs)

# Real-Patient What-if Architecture

## Current Vs Target Architecture

Current implementation:

- `clinic` stores real-patient demographics, assessments, therapies, cytogenetics, and symptom snapshots.
- `simulator` stores educational scenarios and the existing mathematical solver path.
- `simulator.twin.build_patient_twin()` derives a transient patient twin from a single assessment.
- `simulator.models.SimulationAttempt` persists educational scenario runs and writes learner-facing artifacts.

Target implementation in this repository:

- `clinic` remains the source of truth for real-patient longitudinal data.
- `simulator` remains the educational scenario and solver engine.
- `twin_engine` adds the research-only layer for persisted patient twin state, residual calibration, mechanistic what-if runs, causal assumption manifests, and provenance metadata.
- research artifacts are written under `media/research_*` with pseudonymized filenames and no direct identifiers in payloads.

## Data Contract

The research path requires:

- one `clinic.Patient`
- at least one `clinic.Assessment`
- optional `clinic.PatientTherapy` rows with structured `doses`, `cycle_length_days`, and `days_on`
- one or more `twin_engine.PatientTwinState` rows
- optional `twin_engine.ObservationResidual` rows after calibration
- optional `twin_engine.CounterfactualRun` rows for what-if branches
- optional `twin_engine.CausalAssumptionSet` rows for estimand documentation
- `twin_engine.SimulationRunMetadata` for provenance and reproducibility

## Mathematical Model Contract

The research path reuses the current ODE solver in `simulator.models_simulation.MathematicalModel`.

Inputs passed by `twin_engine.simulation_bridge`:

- `baseline_tumor_cells`
- `baseline_healthy_cells`
- `growth_rates.tumor`
- `growth_rates.healthy`
- `drug_doses`
- `interaction_strength`
- `immune_compromise_index`
- optional carrying-capacity overrides
- optional PK/PD schedule functions from `simulator.pharmaco.registry`

Outputs summarized by the bridge:

- tumor reduction
- healthy loss
- time to recurrence
- durability index
- nadir summary
- milestone snapshots
- AUC by drug

All model-derived outputs are labeled `research simulation`.

## Observation Model Formulas

Implemented in `twin_engine.observation_model`.

Predicted biomarkers:

```text
M_hat = alpha_M * tumor_cells + beta_M
FLC_hat = alpha_F * (tumor_cells / max(T_ref, eps)) ** gamma_F
Hb_hat = Hb_baseline * (healthy_cells / max(H_ref, eps)) ** eta_H
```

Default observation parameters are explicit and serialized into the twin state under `parameters.observation`.

## Calibration Objective

Implemented in `twin_engine.calibration` using `scipy.optimize.minimize`.

Initial calibrated parameters:

- `tumor_growth_rate`
- `immune_compromise_index`
- `carrying_capacity_tumor`
- optional `observation.alpha_M`

Objective:

```text
L(theta) = sum_t sum_j w_j * (y_obs[t,j] - y_pred[t,j; theta])^2
```

Residual diagnostics are persisted in `twin_engine.ObservationResidual`.

## Counterfactual Definition

Implemented in `twin_engine.counterfactual`.

Supported intervention definition keys:

- `alternative_regimen_id`
- `drug_doses`
- `start_day`
- `duration_days`
- `schedule`
- `parameter_overrides`
- `random_seed`

Each counterfactual report explicitly separates:

- mechanistic model counterfactual
- causal inference estimand
- unvalidated exploratory hypothesis

## Causal Interpretation Limits

Implemented in `twin_engine.causal`.

Supported classifications:

- `mechanistic_simulation_only`
- `causal_estimand_defined_not_identified`
- `causal_estimand_identified_under_assumptions`
- `insufficient_data`

Important limit:

- a single-patient mechanistic simulation does not establish causal proof
- any identified estimand remains conditional on the explicitly stored graph and assumptions

## Privacy Constraints

Implemented in `twin_engine.validators`.

Rules:

- direct identifiers must not appear in research artifacts
- blocked keys include `first_name`, `last_name`, `mrn`, `notes`, `schedule_notes`, and `birth_date`
- research views enforce owner/demo/staff access
- research logs use `patient_id` and `state_id`, not names or MRNs
- artifact filenames are pseudonymized hashes, not names or MRNs

## Reproducibility And Provenance Contract

Implemented in `twin_engine.provenance`, `twin_engine.run_manifest`, and the
`ResearchRunManifest` / `SimulationRunMetadata` models. See
`docs/research/RUN_MANIFESTS.md` for the full contract.

Every recorded research run captures the complete immutable version vector,
including:

- input hash
- output hash
- code commit hash
- twin config hash
- drug preset hashes
- model version
- solver name
- solver parameters
- random seed when provided
- timestamp via model creation time
- application, API, database-schema, dataset/subset, model-card, evidence-graph,
  validation-protocol, report-template, dependency-lock, registry, intended-use,
  and epistemic identities
- a matching manifest artifact and hashes for every persisted artifact

## How To Run Commands

Initialize a twin from an assessment:

```bash
./venv/bin/python manage.py initialize_patient_twin --patient-id 1 --assessment-id 10
```

Dry-run initialization:

```bash
./venv/bin/python manage.py initialize_patient_twin --patient-id 1 --assessment-id 10 --dry-run
```

Calibrate from history:

```bash
./venv/bin/python manage.py calibrate_patient_twin --patient-id 1
```

Run a what-if counterfactual:

```bash
./venv/bin/python manage.py run_patient_whatif --patient-id 1 --intervention-json /path/to/intervention.json --horizon-days 90
```

Audit research readiness:

```bash
./venv/bin/python manage.py audit_patient_research_readiness --patient-id 1
```

## Test Commands

Focused research suite:

```bash
./venv/bin/python manage.py test \
  twin_engine.tests.test_models \
  twin_engine.tests.test_services \
  twin_engine.tests.test_views \
  twin_engine.tests.test_regressions
```

Existing clinic and simulator regression slice:

```bash
./venv/bin/python manage.py test \
  clinic.tests.test_patient_crud \
  clinic.tests.test_security_and_ux \
  simulator.tests.test_twin_pr1_integration \
  simulator.tests.test_simulation
```

## What Remains Unvalidated

- the solver remains a research simulator, not a validated clinical digital twin
- no ingestion pipeline from external EHR/LIS systems is implemented yet
- imaging, genomics, and MRD are represented only as contract placeholders today
- causal identification is documentation-first and does not infer proof from single-patient simulation

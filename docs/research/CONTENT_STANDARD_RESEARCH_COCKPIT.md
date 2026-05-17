# Research Cockpit Content Standard

This standard keeps the research cockpit readable, auditable, and honest about limits. It applies to cockpit sections, result pages, provenance pages, and developer-facing validation views.

## Required Section Schema

Every major section should include:

- Title: short, task-oriented label.
- One-sentence purpose: what the section helps the user decide or inspect.
- Question answered: the concrete question the section answers.
- Inputs: database records, model state, artifacts, or assumptions used.
- Computation: counts, transformations, formulas, model calls, or checks performed.
- Output interpretation: how to read the displayed values.
- What you can do next: concrete next step, link, command, or data action.
- Limitations: what the section does not prove.
- Developer details: raw JSON, IDs, hashes, artifact paths, or diagnostic payloads collapsed by default.

Default layer is operator-readable. Technical layer may include formulas, model assumptions, data lineage, and definitions. Developer layer contains raw payloads only inside collapsed details blocks.

## Section Examples

### Data Availability

- Question answered: Do we have enough structured longitudinal information to run and interpret simulations?
- Inputs: Assessment, LongitudinalLabResult, PatientTherapy, TherapyInterruption, AdverseEvent, PatientTwinState, ObservationResidual, CounterfactualRun, CausalAssumptionSet.
- Computation: counts and date ranges from the database for the pseudonymized patient.
- Output interpretation: READY means sufficient for the current workflow; PARTIAL means usable but missing relevant fields; MISSING means no usable records.
- What you can do next: add labs, backfill structured records, add therapy schedules, initialize twin, or rerun calibration.
- Limitations: completeness does not imply clinical validity.
- Developer details: record IDs, query payloads, and raw availability metadata.

### Twin Inputs Over Time

- Question answered: What observed values were available before model-generated simulations?
- Inputs: structured labs, assessment fallback values, therapy intervals, interruptions, adverse events.
- Computation: observed values are grouped into disease, hematology, and toxicity panels.
- Output interpretation: lines show observations, not predictions; missing lines mean unavailable structured records.
- What you can do next: backfill LongitudinalLabResult rows, verify assessment fields, or review source quality.
- Limitations: observation plots do not prove disease biology or treatment effect.
- Developer details: chart JSON, source labels, and point IDs.

### Initialize Twin

- Question answered: What mathematical starting state is created from one assessment?
- Inputs: one Assessment and available modeled markers.
- Computation: Assessment(t0) -> risk mapping -> initial tumor/healthy/immune parameters -> PatientTwinState.
- Output interpretation: a reproducible starting state for research simulation.
- What you can do next: initialize from the recommended assessment, then calibrate.
- Limitations: initialization is not treatment simulation, calibration, counterfactual analysis, or diagnosis.
- Developer details: selected assessment ID and model parameters.

### Calibration Quality

- Question answered: After initialization, how well does the model reproduce observed biomarker data?
- Inputs: PatientTwinState, ObservationResidual rows, calibration diagnostics.
- Computation: residual = observed value - predicted value; RMSE and MAE summarize residual size.
- Output interpretation: lower residual after calibration means better fit to available observations.
- What you can do next: add observations, rerun calibration, inspect residuals.
- Limitations: better fit does not prove clinical validity or future predictive accuracy.
- Developer details: optimizer status, objective values, n_observations, n_parameters, residual payload.

### What-if Scenarios

- Question answered: What simulated outcomes change when the intervention schedule changes?
- Inputs: current twin state, recorded baseline therapy, intervention definition, toxicity constraints, solver outputs.
- Computation: rerun the mechanistic model with alternative intervention input.
- Output interpretation: tumor reduction, healthy loss, durability, predicted biomarkers, toxicity penalty, and research utility are model outputs or heuristics.
- What you can do next: compare scenarios, inspect trajectories, review toxicity and causality panels.
- Limitations: research utility is a heuristic ranking, not a treatment recommendation.
- Developer details: run ID, report artifact, trajectory artifact, raw comparison metrics.

### Trajectory Comparison

- Question answered: How do simulated trajectories evolve over the prediction horizon?
- Inputs: trajectory artifacts from completed runs.
- Computation: plot baseline and alternative tumor/healthy compartments over simulation time.
- Output interpretation: curves show model state evolution under solver assumptions.
- What you can do next: inspect schedule-collapse warnings or improve exposure bridge resolution.
- Limitations: identical curves are a model-resolution limitation, not evidence of biological equivalence.
- Developer details: artifact paths and trajectory arrays.

### Toxicity Constraints

- Question answered: What safety-related constraints affect interpretation of simulated schedules?
- Inputs: AST, ALT, neutropenia, infection, interruptions, adverse events.
- Computation: descriptive summaries and heuristic penalties.
- Output interpretation: descriptive_only means observed toxicity context, not predicted toxicity dynamics.
- What you can do next: add missing toxicity observations, review interruptions, improve toxicity model.
- Limitations: lower exposure penalty does not prove lower toxicity.
- Developer details: toxicity summary payload and event IDs.

### Causality Status

- Question answered: Is this a causal effect estimate or a mechanistic model branch?
- Inputs: CausalAssumptionSet and run classification metadata.
- Computation: classify outputs as mechanistic-only unless a separate identification design exists.
- Output interpretation: mechanistic simulation may be present while causal effect estimation is absent.
- What you can do next: define estimand, assumptions, data design, adjustment strategy, and sensitivity analysis.
- Limitations: single-patient model branches do not identify population causal effects.
- Developer details: graph, variables, adjustment set, identification status.

### Scientific And Mathematical Basis

- Question answered: Which implemented components, assumptions, and references support the displayed outputs?
- Inputs: model_references.json and component metadata.
- Computation: display repository-declared implementation and evidence status.
- Output interpretation: missing citation means no peer-reviewed source is attached in the repository.
- What you can do next: attach reviewed references, update assumptions, or mark validation status.
- Limitations: documented assumptions do not prove external validity.
- Developer details: component key, variables, validity domain, next review action.

### Provenance

- Question answered: Which data, model state, intervention file, parameters, and artifacts produced this result?
- Inputs: pseudonymized patient, observations, twin state, diagnostics, intervention definition, run metadata, artifacts.
- Computation: trace chain and hashes where available.
- Output interpretation: provenance supports reproducibility and audit.
- What you can do next: inspect run metadata, verify artifacts, regenerate missing reports.
- Limitations: provenance proves traceability, not clinical validity.
- Developer details: metadata ID, hashes, artifact paths, raw payloads.

### Developer Console

- Question answered: Is the research workflow internally coherent before commit, demo, or analysis?
- Inputs: database records, model outputs, references, privacy checks, git staging state, artifacts.
- Computation: grouped offline checks with pass/warn/fail status.
- Output interpretation: warnings identify review tasks; failures identify blockers.
- What you can do next: apply suggested fixes, rerun checks, then run the pre-push safety gate.
- Limitations: developer checks are not clinical validation.
- Developer details: raw check payloads and object IDs.

## Forbidden Language

Do not use user-facing language that:

- labels a regimen as the clinically recommended option
- calls one therapy the best clinical choice
- tells a clinician what to prescribe
- claims clinical superiority without validation
- claims validated decision-support status
- claims a causal effect has been proven when the platform only produced a mechanistic branch

## Allowed Language

Use language like:

- "research simulation"
- "mechanistic model"
- "heuristic utility"
- "future decision-copilot prototype"
- "not clinically validated"
- "not a treatment recommendation"
- "causal effect not identified"
- "descriptive constraint"

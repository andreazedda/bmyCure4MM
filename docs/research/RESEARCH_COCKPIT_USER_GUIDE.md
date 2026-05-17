# Research Cockpit User Guide

The research cockpit is available at:

- `/research/patient/<patient_id>/cockpit/`
- `/research/glossary/`

It is a pseudonymized, research-only workspace for inspecting longitudinal data, twin state status, calibration diagnostics, and mechanistic what-if simulations. It is not clinically validated, does not estimate causal effects, and must not be used as a prescribing system.

## Suggested Walkthrough

1. Start with **Workflow map** to see which research steps are ready, partial, or missing.
2. Open **Data availability** and confirm structured observations, therapies, interruptions, adverse events, residuals, runs, and provenance exist.
3. Review **Twin state / Initialize Twin**. Initialization means `Assessment(t0) -> risk mapping -> tumor/healthy/immune parameters -> PatientTwinState`.
4. Review **Calibration quality** before interpreting what-if rows. Residuals use `residual = observed value - predicted value`; RMSE and MAE summarize residual size.
5. Review **Twin Inputs over time** to understand which observed biomarkers and safety markers were available before simulation.
6. Compare **What-if scenarios** only after reading the metric definitions and toxicity constraints.
7. Compare **Trajectory comparison** and watch for schedule-resolution limitation warnings.
8. Read **Causality status** before making any causal interpretation. Current simulations are `Y_model(a') = f(x_t, theta_hat, a')`, not `E[Y | do(A=a')] - E[Y | do(A=a)]`.
9. Read **Scientific and mathematical basis** for each component's source type, assumptions, variables, validity domain, limitations, and next review action.
10. Use **Provenance** and **Developer checks** for traceability before sharing artifacts or pushing code.

## Main Sections

- Data available: counts assessments, structured labs, therapies, interruptions, adverse events, residuals, and completed runs.
- Twin state / Initialize Twin: shows the current twin state and recommends the earliest assessment with the strongest modeled marker completeness.
- Calibration quality: reports persisted residuals and before/after RMSE when calibration diagnostics exist.
- Longitudinal markers: charts structured `LongitudinalLabResult` rows first, with assessment fallback only when a matching structured value is absent.
- What-if scenarios: shows the latest completed run per scenario label, sorted by heuristic `research_utility`.
- Trajectory comparison: plots baseline and alternative trajectories from JSON trajectory artifacts.
- Toxicity constraints: summarizes descriptive toxicity history and documented interruptions.
- Causality status: states whether assumption sets exist and keeps mechanistic simulations separate from causal estimands.
- Scientific and mathematical basis: lists implemented and placeholder model components from `twin_engine/model_references.json`.
- Developer checks: links staff users to the internal console.
- Glossary: explains terms such as twin, initialization, calibration, residual, counterfactual, causal effect, toxicity constraint, research utility, provenance, schedule collapse, exposure bridge, longitudinal lab, adverse event, and therapy interruption.

## Scenario Row Interpretation Example

A scenario row is a mechanistic branch, not an instruction. Read it like this:

- **Scenario**: the intervention label used to rerun the model.
- **Status**: whether the latest run completed.
- **Tumor reduction**: model-estimated disease-burden reduction over the horizon.
- **Healthy loss**: model-estimated healthy-compartment loss over the horizon.
- **Durability**: model-derived persistence term.
- **Predicted biomarkers**: horizon biomarker projections from the model, not observed labs.
- **Toxicity penalty**: descriptive heuristic derived from observed toxicity context.
- **Research utility**: `tumor_reduction + (1 - healthy_loss) + durability_index - toxicity_constraint_penalty`.

## Interpretation Rules

- Scenario rankings are heuristic research utilities, not treatment rankings.
- Near-identical trajectories can indicate schedule-resolution collapse in the exposure bridge.
- Toxicity summaries are descriptive unless a future predictive toxicity model is implemented and validated.
- Raw JSON artifacts are available behind developer details blocks for auditability.
- Missing chart lines mean unavailable structured records, not absence of disease.
- Toxicity penalty is descriptive unless a future predictive toxicity model is implemented and validated.
- Causal effect not identified unless a separate design, estimand, assumptions, data source, and sensitivity analysis support it.

## Patient 4 Local Baseline

The local pseudonymized patient 4 workflow currently has structured assessments, labs, therapies, interruptions, adverse events, residuals, causal assumptions, provenance, and completed counterfactual runs. The cockpit intentionally displays this as `Research Patient 4` rather than direct identifiers.
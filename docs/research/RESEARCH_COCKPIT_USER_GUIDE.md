# Research Cockpit User Guide

The research cockpit is available at:

- `/research/patient/<patient_id>/cockpit/`

It is a pseudonymized, research-only workspace for inspecting longitudinal data, twin state status, calibration diagnostics, and mechanistic what-if simulations. It is not clinically validated, does not estimate causal effects, and must not be used as a prescribing system.

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

## Interpretation Rules

- Scenario rankings are heuristic research utilities, not treatment rankings.
- Near-identical trajectories can indicate schedule-resolution collapse in the exposure bridge.
- Toxicity summaries are descriptive unless a future predictive toxicity model is implemented and validated.
- Raw JSON artifacts are available behind developer details blocks for auditability.

## Patient 4 Local Baseline

The local pseudonymized patient 4 workflow currently has structured assessments, labs, therapies, interruptions, adverse events, residuals, causal assumptions, provenance, and completed counterfactual runs. The cockpit intentionally displays this as `Research Patient 4` rather than direct identifiers.
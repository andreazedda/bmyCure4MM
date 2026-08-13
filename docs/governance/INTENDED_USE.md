# Canonical Intended Use

Status: **current and authoritative**
Claims policy: `claims-policy-v1`

```yaml
intended_use_level: E1_research_prototype
clinical_decision_support: false
patient_specific_prediction_validated: false
causal_effect_identified: false
```

bmyCure4MM is a computational research platform for mechanistic simulation,
data-readiness assessment, reproducible scenario comparison, hypothesis
generation, and falsification in Multiple Myeloma research.

It is not a medical device, prescribing system, clinical decision-support
system, or validated patient-specific predictor. Its outputs cannot authorize a
treatment, regimen, dose, schedule, interruption, or duration change.

## Current surfaces

- `/learn/`: synthetic educational scenarios and tutorials.
- `/research/`: mechanistic twin, lineage, data-readiness, and counterfactual
  research views.
- `/admin/`: platform administration.
- Clinic records provide data navigation; their presence does not confer
  clinical decision authority.

## Permitted uses

- reproduce a version-bound scientific run;
- inspect observed, derived, simulated, heuristic, and hypothetical content;
- compare mechanistic scenarios under explicit assumptions;
- detect missing data or model-readiness gaps;
- report negative or inconclusive results;
- formulate hypotheses that require independent testing.

## Not permitted by current evidence

- choosing or ranking treatment for a patient;
- issuing patient-specific dose or schedule instructions;
- predicting patient benefit or individual prognosis;
- claiming that a regimen is clinically superior;
- interpreting a mechanistic counterfactual as an identified causal effect.

The canonical implementation constants live in `mmportal/governance.py` and the
machine-readable release declaration lives in `governance/release_claims.json`.

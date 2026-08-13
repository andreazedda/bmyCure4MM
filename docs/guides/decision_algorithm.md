# Model-relative diagnostic thresholds

Source: `simulator/decision_algorithm.py`.

The declarative catalog classifies simulated KPIs and emits research workflow
flags. It is not a treatment-selection algorithm.

## Semantics

- Thresholds are `HEURISTIC` model-navigation rules.
- Inputs to the evaluator are `SIMULATED` metrics.
- Rule outputs request inspection of uncertainty, provenance, or alternative
  hypothetical configurations.
- R-ISS and cytogenetic entries provide population context only.
- No rule estimates patient benefit or identifies a causal treatment effect.

The machine-readable response includes algorithm version, intended-use level,
epistemic label, model version, and the three false clinical-claim flags.

Threshold changes alter output-language semantics and require explicit review,
tests, and versioning. They do not change the underlying ODE/PK/PD equations.

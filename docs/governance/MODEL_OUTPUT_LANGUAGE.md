# Model Output Language

Status: **current and authoritative**

## Required output context

Every major rendered or machine-readable scientific output must expose:

- `intended_use_level`;
- `epistemic_label`;
- `model_version`;
- the relevant input/run lineage when available;
- a limitation that distinguishes mechanistic simulation from clinical and
  causal conclusions.

## Preferred phrases

Use model-relative descriptions:

- `constraint_flag`;
- `constraint_watch_zone`;
- `simulated_high_impact_zone`;
- `simulated_low_impact_zone`;
- `model_regrowth_signal`;
- `short_model_durability_signal`;
- `insufficient_evidence`.

State the metric, threshold, units, and epistemic status. For example:

> `constraint_flag` (`HEURISTIC`): simulated healthy-cell loss is at or above
> the configured 0.30 threshold in model `research-twin-v1`.

## Prohibited transformations

Do not translate a threshold into a patient instruction. Do not turn simulated
chronology into causality. Do not call a mechanistic scenario “benefit,”
“success,” “failure,” “safe,” or “effective” without a precise model-relative
definition and the appropriate evidence gate.

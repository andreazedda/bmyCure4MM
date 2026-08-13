# Research Simulator

Status: current under `claims-policy-v1`.

The simulator compares mechanistic scenarios under explicit assumptions. Its
tumor, healthy-cell, exposure, toxicity, and recurrence outputs are `SIMULATED`.
Thresholds attached to them are `HEURISTIC` model-relative diagnostic flags.

## What a run can show

- trajectories produced by the selected model, parameters, and schedule;
- whether a configured model-relative threshold is crossed;
- sensitivity to schedules, parameters, and uncertainty assumptions;
- reproducible artifacts when model/input/version identity is recorded.

## What a run cannot show

- clinical therapy selection;
- whether a dose or schedule should change;
- validated patient benefit or prognosis;
- a causal effect of treatment.

Use `/learn/` for synthetic tutorials and `/research/` for lineage-bound
research workflows. See [Canonical Intended Use](../governance/INTENDED_USE.md)
and [Model Output Language](../governance/MODEL_OUTPUT_LANGUAGE.md).

# Simulator Walkthrough

## Overview
- Coupled PK/PD + logistic growth ODE solved via `MathematicalModel`.
- HTMX form groups: Baseline, Doses, Horizon, Advanced.
- Novice mode surfaces bilingual hints, range badges, optimization shortcut.

## Workflow
1. Validate labs (ANC, platelets, pregnancy) – hard blocks prevent unsafe runs.
2. Confirm preset defaults; use **Reset to preset** if values drift.
3. Enter optional **Random seed** for reproducible cohorts.
4. Choose cohort size and time horizon; warning when >300 days.
5. Submit – solver writes CSV/Plot/Twin JSON artifacts for traceability.

## Outputs
- KPIs with explainable tooltips (tumor reduction, healthy loss, AUC, time to recurrence).
- Artifacts stored under `media/sim_data` and `media/sim_plots`.
- Warnings highlight toxicity or regrowth.

## Tips
- Interaction strength badge turns yellow (0.1–0.2) or red (>0.2).
- Always log rationale in the free-text note for audit trails.

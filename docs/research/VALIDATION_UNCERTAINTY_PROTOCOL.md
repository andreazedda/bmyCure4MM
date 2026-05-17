# Validation and Uncertainty Protocol

This protocol defines the internal research diagnostics used to answer five questions in the mechanistic simulation cockpit:

1. How stable are scenario outputs under parameter uncertainty?
2. How well do held-out observed biomarker points agree with the model?
3. Which scenario ranking is robust versus fragile?
4. Which parameters drive the output most?
5. Which outputs should be treated as exploratory only?

## Guardrails

- Research simulation only.
- Do not modify real patient source data while running these diagnostics.
- Do not weaken privacy protections or export direct identifiers into artifacts.
- Do not claim clinical validity.
- Do not claim causal effect estimation.
- Keep wording in cockpit and reports limited to mechanistic prototype, research simulation, or future decision-copilot prototype framing.

## Diagnostic Layers

### 1. Scenario uncertainty

- Input: current twin state parameters, intervention definition, and heuristic parameter uncertainty scales when calibrated distributions are unavailable.
- Method: truncated-normal Monte Carlo perturbations of model parameters, observation parameters, and prototype toxicity coefficients.
- Output: point estimate, mean, median, `p05`, `p25`, `p75`, `p95`, interval width, coefficient of variation, and narrow/moderate/wide classification for scenario outputs.
- Interpretation: heuristic intervals stress-test ranking stability but do not claim posterior calibration unless an explicit calibrated parameter distribution is stored.

### 2. Rolling-origin backtesting

- Input: ordered historical assessments.
- Method: fit the mechanistic model in memory on training observations only, then predict the next held-out assessment.
- Output: fold rows with training window, held-out prediction, observed values, residuals, RMSE, and MAE; aggregate RMSE/MAE by biomarker and overall.
- Interpretation: backtests quantify held-out agreement at observed biomarker points only. They do not establish clinical validity.

### 3. One-at-a-time sensitivity

- Input: current twin state parameters and selected intervention knobs such as dose and duration when available.
- Method: perturb one driver at a time by small positive and negative fractions while holding other quantities fixed.
- Output: ranked drivers by absolute `utility_v2` movement, direction labels, elasticity, and unstable-parameter flags.
- Interpretation: a large sensitivity response marks fragile exploratory ranking and should lower trust in a point estimate.

### 4. Robust ranking

- Input: stored uncertainty sample traces across multiple scenarios.
- Method: align sample indices, compute probability-best, pairwise dominance, interval overlap, and compare point rank to robust rank.
- Output: per-scenario point rank, robust rank, probability-best, interval overlap, and robust/contested/fragile classification.
- Interpretation: robust rank is still an internal research diagnostic. It does not imply treatment recommendation or causal effect.

## Trust Labels

- `internally_checked_research_signal`: uncertainty stored, backtest present, and robust ranking not fragile.
- `fragile_research_signal`: scenario rank shifts under uncertainty overlap or sensitivity.
- `research_signal_pending_validation`: some diagnostics exist, but backtest or robust ranking is missing.
- `exploratory_only`: uncertainty is wide or scenario uncertainty is missing.

## Operational Notes

- Commands write summaries into existing JSON-bearing records only.
- Scenario-level uncertainty and sensitivity live on `CounterfactualRun.comparison_metrics`.
- Patient-level backtesting and robust-ranking summaries are attached through `SimulationRunMetadata.solver_parameters.diagnostic_summary`.
- Developer checks warn when these summaries are missing or when uncertainty is wide.
# Optimization Lab

## Objectives
- Maximize **tumor_reduction**
- Minimize **healthy_loss** (constraint healthy_loss ≤ 0.25)
- Minimize **total_drug_exposure** (Σ AUC)

## Algorithm
- MOTPE (Optuna) sampler with configurable `n_trials` and seed.
- Penalties applied for toxicity violations or schedule breaches.
- Pareto front sorted by efficacy desc, safety asc.

## Reproducibility
- Pass a custom seed in the Optimization Lab form or Simulation form.
- Use **Export study** (coming soon) to capture params + KPIs as JSON.
- Artifacts stored with attempt metadata for auditing.

## Tutorial: R-ISS III scenario
1. Select high-risk scenario, enable Patient Twin.
2. Run Optimization Lab with 50 trials, seed 2025.
3. Inspect Pareto table: pick 3 non-dominated solutions (green badges).
4. Click *Simulate solution* to push params into the main form and rerun.
5. Compare KPIs vs baseline attempt; keep artifact links for review.

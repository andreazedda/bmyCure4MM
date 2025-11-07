# Predictive Lab Quickstart

## Goal & Constraints
- Enable non-experts to simulate MM regimens safely under preset guardrails
- Respect preset dose ranges, toxicity limits (healthy_loss ≤ 0.25), reproducible seeds

## Guided Steps
1. **Select scenario** – review patient summary, R-ISS, labs.
2. **Enable Patient Twin** – keeps biology synced with LDH/β2M/FLC.
3. **Pick preset** – VRd, Dara-Rd, KRd; presets lock safe dose bounds.
4. **Adjust doses/schedules** – stay within displayed ranges and scheduling patterns.
5. **Choose horizon + cohort** – default 84 days, cohort 1; increase cohort (10/50/200) to expose uncertainty bands.
6. **Run simulation** – watch spinner + toast, then inspect KPIs with tooltips.
7. **Interpret KPIs** – tumor reduction ≥90% and healthy loss ≤25% = strong response.
8. **Launch Optimization Lab** – MOTPE search with toxicity constraint, then re-simulate Pareto candidates.
9. **Iterate** – tweak inputs, rerun; export CSV/Plot for audit.

## FAQ
- **AUC vs toxicity?** Higher AUC (drug exposure) often correlates with immunosuppression; keep low while maintaining efficacy.
- **Healthy loss threshold?** Guardrail at 0.25 (25%). Yellow badges warn at 0.2; red at 0.3.
- **When to increase cohort size?** Use 10 for sensitivity analysis, 50+ when presenting uncertainty bands to stakeholders.

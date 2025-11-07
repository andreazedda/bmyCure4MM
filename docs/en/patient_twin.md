# Patient Twin

## Inputs
- Assessment fields: R-ISS, LDH, β2M, FLC ratio.
- Optional immune markers feed into `twin_risk.yaml` weights.

## Mapping
- Weighted risk score → `_lerp` mapping for tumor/healthy growth, carrying capacity, immune compromise.
- Coefficients editable in `simulator/presets/twin_risk.yaml`.

## Tutorial
1. Create Assessment A (R-ISS I, LDH normal, β2M 3, FLC 0.9) → expect low risk.
2. Create Assessment B (R-ISS III, LDH 600, β2M 12, FLC 12) → high risk.
3. Run simulations with `use_twin=True`; observe tumor growth increases and healthy capacity shrinks from A→B.
4. Inspect saved `twin_params.json` artifact for provenance.

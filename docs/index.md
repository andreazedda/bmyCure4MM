# bmyCure4MM documentation

bmyCure4MM is an `E1_research_prototype` for reproducible computational
multiple-myeloma research. It stores structured research records, constructs
versioned model state, and runs mechanistic hypothetical scenarios.

It is not clinical decision support. Patient-specific prediction is not
validated, and causal treatment effects are not identified.

## Start with governance

1. [Canonical intended use](governance/INTENDED_USE.md)
2. [Claims policy](governance/CLAIMS_POLICY.md)
3. [Epistemic labels](governance/EPISTEMIC_LABELS.md)
4. [Model output language](governance/MODEL_OUTPUT_LANGUAGE.md)

## Current research surfaces

- **Clinic** records pseudonymized structured inputs and dated observations.
- **Research Twin** stores model state and lineage to its input snapshot.
- **Simulator** evaluates hypothetical mechanistic configurations.
- **Scientific Cockpit** exposes diagnostics, assumptions, and provenance.
- **ChemTools** is an optional exploratory cheminformatics surface.

Every displayed simulation must be read as `SIMULATED` or `HYPOTHETICAL` and
must retain its stated model version and limitations. A successful process exit
means computation completed; it does not establish scientific or clinical
validity.

## Contributor path

- [Architecture](guides/architecture.md)
- [Database](guides/database.md)
- [Operations](guides/operations.md)
- [Development](development/DEVELOPMENT.md)
- [Research simulator](en/simulator.md)

Historical material with superseded claims is isolated under
[`archive/`](archive/LEGACY_PRESCRIPTIVE_OUTPUTS.md) and is not a source of
current product claims.

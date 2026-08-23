# Changelog

All current release claims are governed by [`docs/governance/CLAIMS_POLICY.md`](docs/governance/CLAIMS_POLICY.md). A changelog entry records implementation history; it does not establish scientific validity.

## Unreleased — M0-R Canonical Scientific Baseline

### Documentation and planning

- Added a canonical mission: virtual patients, current and experimental therapy hypotheses, falsification, reproducible protocols and an interdisciplinary Knowledge Commons.
- Separated current verified architecture from the target virtual-laboratory architecture.
- Added a current model-status catalog with formulas, units, examples, code mappings and falsification boundaries.
- Added source-of-truth and documentation-governance policies plus a corpus inventory.
- Reconciled GitHub roadmap/tracker state and Notion project state at canonical `master` `bf097810...`.
- Added strict MkDocs validation to documentation pull requests without deploying PR content.

### Scientific identity and reproducibility

- Added the Structured Research Dataset Contract v0.1 and generic idempotent importer.
- Added content-addressed computational-input and Patient Twin lineage contracts.
- Added immutable research run manifests and a governed model registry.
- Added typed scientific payload contracts, artifact integrity checks, direct-comparability rules and append-only invalidation.
- Locked the dependency graph with `uv.lock`, Python 3.11/3.12 support and numerical/repository/dependency checks.
- Recorded newly disclosed Django/sqlparse advisories as a blocking current condition under issues `#69` and `#70`; no dependency fix or audit bypass is claimed by this documentation change.

### Intended use and language

- Declared `E1_research_prototype` as the current intended-use level.
- Declared clinical decision support, validated patient-specific prediction and identified causal effects as false.
- Replaced prescriptive simulation interpretation with model-relative diagnostic flags carrying epistemic and model-version context.
- Added machine-readable release claims and semantic regression tests.
- Separated learning, research and administration navigation.
- Invalidated pre-policy chemistry efficacy, survival, toxicity/risk-benefit and patient-prognosis outputs for scientific interpretation; numerical equations and parameters were not changed by that language migration.

### Operational identity

- Added a governed temporary least-privilege M0 smoke-user command. It does not constitute a normal product login or reviewed RBAC contract.

## Historical 1.0.0 claim — 2024-11

A 2024 document described the repository as ready for a public `1.0.0` release and incorrectly equated passing 30 of 32 tests with 94% code coverage. Those statements are historical and are not current release, coverage or security evidence.

The archived record is retained at [`docs/archive/releases/RELEASE_PREPARATION_2024.md`](docs/archive/releases/RELEASE_PREPARATION_2024.md). Legacy prescriptive claims are summarized in [`docs/archive/LEGACY_PRESCRIPTIVE_OUTPUTS.md`](docs/archive/LEGACY_PRESCRIPTIVE_OUTPUTS.md).

---
title: bmyCure4MM Documentation
status: CANONICAL
owner: Andrea Zedda
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# bmyCure4MM documentation

bmyCure4MM is an `E1_research_prototype` for Multiple Myeloma computational research. Its target is a virtual and reproducible laboratory for virtual patients, current and experimental therapy hypotheses, falsification, reproducible protocols and an interdisciplinary Knowledge Commons.

It is not clinical decision support. Patient-specific prediction is not validated and causal treatment effects are not identified.

## First reading

1. [Mission and north star](product/MISSION_AND_NORTH_STAR.md)
2. [Current verified state](product/CURRENT_STATE.md)
3. [Capabilities and limitations](product/CAPABILITIES_AND_LIMITATIONS.md)
4. [Canonical intended use](governance/INTENDED_USE.md)
5. [Source-of-truth hierarchy](governance/SOURCE_OF_TRUTH.md)

## Architecture

- [Current system](architecture/CURRENT_SYSTEM.md)
- [Target virtual laboratory](architecture/TARGET_VIRTUAL_LAB.md)
- [Current model registry and formulas](models/REGISTRY.md)

Current and target architecture are deliberately separated. A target component is not an implemented capability.

## Research contracts

- [Structured dataset contract](research/STRUCTURED_DATASET_CONTRACT_V0_1.md)
- [Twin lineage](research/TWIN_LINEAGE.md)
- [Immutable run manifests](research/RUN_MANIFESTS.md)
- [Comparability and invalidation](research/COMPARABILITY_AND_INVALIDATION.md)
- [Real-patient Research Loop](research/RESEARCH_LOOP.md)
- [Validation and uncertainty](research/VALIDATION_UNCERTAINTY_PROTOCOL.md)

## Governance and operations

- [Claims policy](governance/CLAIMS_POLICY.md)
- [Epistemic labels](governance/EPISTEMIC_LABELS.md)
- [Model output language](governance/MODEL_OUTPUT_LANGUAGE.md)
- [Documentation policy](governance/DOCUMENTATION_POLICY.md)
- [Authentication and access](operations/AUTHENTICATION_AND_ACCESS.md)
- [Testing](operations/TESTING.md)
- [Dependency operations](operations/DEPENDENCIES.md)

A successful process exit means computation completed. It does not establish scientific validity, external validity or clinical utility.

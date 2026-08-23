---
title: Current Verified State
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: all contributors
last_verified_at: 2026-08-24
last_verified_git_sha: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
source_of_truth: code, registries, schemas, migrations and current validation evidence
---

# Current verified state

## Canonical identity

```yaml
repository: andreazedda/bmyCure4MM
branch: master
verified_head: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
intended_use_level: E1_research_prototype
clinical_decision_support: false
patient_specific_prediction_validated: false
causal_effect_identified: false
```

This page describes the repository state at the SHA above. Later changes require a new verification record.

## Current components

| Component | Current responsibility | Status |
|---|---|---|
| `clinic` | Structured patients, assessments, therapies, longitudinal navigation and current authentication entry points | CURRENT_PARTIAL |
| `simulator` | Mechanistic tumour/healthy-cell and PK/PD simulation, schedules, scenarios and experimental therapy-design code | CURRENT_HEURISTIC_PROTOTYPE |
| `twin_engine` | Twin state/lineage, observation residuals, calibration infrastructure, counterfactual research, run identity, provenance and research views | CURRENT_PARTIAL |
| `chemtools` | Optional exploratory cheminformatics and molecular utilities | CURRENT_PARTIAL |
| `docs_viewer` | In-application documentation rendering and feedback | CURRENT_EXECUTABLE |
| `governance` | Machine-readable release claims and intended-use controls | CURRENT_EXECUTABLE |
| `scripts` | Repository hygiene, numerical baseline, dependency audit and research safety tooling | CURRENT_EXECUTABLE |

## Completed M0-R foundations

- Structured Research Dataset Contract v0.1 and generic idempotent importer.
- Content-addressed computational-input and Twin-lineage contracts.
- Canonical E1 intended use and non-prescriptive model-output language.
- Deterministic `uv.lock` dependency graph using Django 5.2.17 LTS, sqlparse 0.6.0 and supported Python 3.11/3.12 tooling.
- Immutable research run manifests and version vectors.
- Governed model registry and versioned scientific JSON contracts.
- Artifact hashing, direct-comparability rules and append-only invalidation semantics.
- Numerical baseline, repository hygiene and dependency-audit infrastructure.
- Synthetic Django deployment checks, migration-drift checks and disposable migration evidence.
- A stable protected `required` result over locked quality, Python matrix and Docker build.

## Current registered models

The machine-readable source is `twin_engine/model_registry.json`.

| Model ID | Version | Current status |
|---|---|---|
| `patient_twin_state_model` | `research-twin-v1` | CURRENT_HEURISTIC_PROTOTYPE |
| `observation_model` | `observation-model-v1` | CURRENT_HEURISTIC_PROTOTYPE |
| `lenalidomide_exposure_model` | `exposure-bridge-v1` | CURRENT_EXECUTABLE with research-only interpretation |
| `hepatic_signal_model` | `toxicity-prototype-v1` | CURRENT_HEURISTIC_PROTOTYPE |
| `neutropenia_signal_model` | `toxicity-prototype-v1` | CURRENT_HEURISTIC_PROTOTYPE |
| `counterfactual_model` | `research-twin-v1` | CURRENT_HEURISTIC_PROTOTYPE |
| `therapy_design_toy_model` | `therapy-design-toy-v1` | CURRENT_HYPOTHETICAL_PROTOTYPE |

Detailed status and formula mappings are in [Current model registry](../models/REGISTRY.md).

## Current research loop capability

The repository can persist structured longitudinal entities, construct lineage-bound Twin states, execute mechanistic simulations, create immutable run identity and compare hypothetical configurations. It also has residual/calibration, uncertainty, sensitivity, robustness and temporal-diagnostic infrastructure.

It has **not** yet demonstrated a source-verified, leakage-free, identifiable and out-of-sample-valid patient-specific Research Loop. M1 remains blocked until the M0-R documentation, final CI-governance and minimum QA baseline closes.

## Current authentication and authorization boundary

Verified current behavior includes public documentation and public simulator surfaces plus protected clinic, research, simulator-management and API routes. The configured login target is the Django admin login. A reviewed normal product login and least-privilege product-role contract is not yet complete.

Known authorization gaps are tracked in GitHub issue `#8`. Admin access or temporary smoke-user creation must not be interpreted as successful normal-user product authentication.

## Current CI evidence

GitHub Actions cover dependency/reproducibility checks, documentation, secret scanning, image workflows, Django 5.2 compatibility, UI screenshots and redeployment.

Merged PRs `#71` and `#72` established:

```text
Django 5.2.17 LTS and sqlparse 0.6.0 identity: PASS
dependency audit: PASS
ordinary and synthetic deploy checks: PASS
migration drift and disposable migrations: PASS
Ruff, format, mypy and repository hygiene: PASS
numerical baseline: PASS / unchanged
Python 3.11 core check/tests: PASS
Python 3.12 core check/tests: PASS
authorization-sensitive compatibility subset: PASS
Secret Scan: PASS
Docker locked build: PASS
protected required result: PASS
```

The workflows are not yet fully rationalized: external Actions remain tag-pinned rather than immutable-SHA-pinned, the research safety script lacks verified base/head PR mode, and PR templates/CODEOWNERS/final branch-protection evidence remain under issue `#13`.

## Privacy boundary

Private patient documents, direct identifiers, source excerpts and private dataset payloads must remain outside Git under ignored `local_private/` storage. Git may contain schemas, generic import logic, synthetic fixtures and explicitly validated de-identified artifacts.

## Immediate execution path

```text
#14 documentation/source-of-truth baseline
→ #13 mandatory CI aggregate gate and PR governance
→ #23 minimum M0-R scientific/privacy invariants
→ #26 close M0-R
```

Private read-only source verification for `#32` may prepare in parallel. New canonical patient-specific calibration or published what-if results may not be promoted until M0-R closes.

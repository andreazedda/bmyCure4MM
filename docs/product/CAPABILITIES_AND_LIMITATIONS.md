---
title: Capabilities and Limitations
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: users, reviewers and contributors
last_verified_at: 2026-08-24
last_verified_git_sha: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
---

# Capabilities and limitations

## Interpretation rule

A capability is usable only inside its stated evidence and intended-use boundary.

```text
computation completed
≠ model verified
≠ model validated
≠ externally validated
≠ clinically useful
```

## Current capability matrix

| Capability | What exists | What it supports | What it does not establish |
|---|---|---|---|
| Structured research data | Versioned schema, content hashes, idempotent import infrastructure | Reproducible data snapshots and semantic identity | Source-level verification of every private-case value |
| Patient Twin state | Persistent state, parameters, uncertainty fields and lineage | Versioned research state and parent/child history | A validated physiological digital twin |
| Observation residuals | Predicted/observed payloads, RMSE and MAE | Calibration diagnostics | Identifiability or future predictive validity |
| Mechanistic simulation | Logistic tumour/healthy dynamics, PK/PD and schedule inputs | Exploratory trajectories and model-relative comparisons | Clinical outcome prediction |
| Exposure bridge | Day-resolved administered-dose vectors, cumulative dose and hashes | Distinguishing schedules with equal average dose | Patient-specific pharmacokinetics or bioavailability |
| Toxicity signals | Heuristic normalized liver and neutropenia signals | Stress testing and research ranking | AST, ALT or ANC prediction; causal attribution |
| Counterfactual research | Baseline/alternative mechanistic runs with immutable identity | Hypothetical configuration comparison | Identified causal treatment effects |
| Therapy-design toy model | Experimental target/modality/logic and Pareto scaffolding | Testing software abstractions and hypotheses | Validated ADC, CAR-T, bispecific or optimal-treatment design |
| Uncertainty/sensitivity | Prototype sampling and diagnostics | Inspecting model fragility | Calibrated population uncertainty |
| Provenance and run identity | Version vector, hashes, comparability and invalidation | Reproduction and audit | Correctness of the biological assumptions |
| Knowledge/documentation | Governance, research contracts and technical guides | Learning and source navigation | A complete curated Multiple Myeloma evidence graph |

## Allowed current conclusions

- A run was executed using a stated code/model/data/configuration identity.
- Two runs are or are not directly comparable under declared compatibility rules.
- A scenario produces a particular trajectory under the current model and assumptions.
- A schedule differs temporally from another schedule.
- Data or provenance are missing for a requested research operation.
- A calibration attempt improved or failed a declared metric, with all other validity questions still open.
- A hypothesis is unsupported, non-identifiable or unstable under sensitivity analysis.

## Forbidden current conclusions

- one regimen is clinically best or superior for a patient;
- a specific dose or schedule should be used;
- a patient will benefit from a simulated intervention;
- a model-relative trajectory is an identified causal effect;
- MRD negativity or a negative imaging result proves zero disease;
- a B-like or stem-like reservoir is an established universal mechanism;
- heuristic toxicity signals predict clinical laboratory values;
- successful optimization validates the selected strategy.

## Main scientific limitations

1. **Observation model:** current mappings compress heterogeneous biomarkers into simplified functions of tumour and healthy-cell states.
2. **Latent-state interpretation:** current Twin states are research abstractions, not directly measured physiological states.
3. **Identifiability:** multiple parameter sets may explain the same sparse longitudinal observations.
4. **External validity:** current models are not externally benchmarked for patient-level prediction.
5. **Disease representation:** the runtime core does not yet represent validated subclonal, spatial, niche, immune, bone and organ dynamics.
6. **Toxicity attribution:** current risk signals do not separate drug effect from disease, infection, organ reserve, co-medications and measurement error.
7. **Causality:** mechanistic counterfactuals remain distinct from identified causal estimands.

## Main engineering and platform limitations

- normal product authentication and least-privilege roles are incomplete;
- object-level authorization remains inconsistent in some patient-derived web/aggregate paths;
- a protected `required` dependency/reproducibility gate exists, but the complete #13 workflow rationalization, immutable Action pinning, exact base/head safety mode, PR templates and CODEOWNERS remain open;
- documentation baseline exists in draft PR `#68` and is not yet canonical on `master`;
- production runtime, PostgreSQL, storage, observability, quotas and backup/restore remain governed by M0-S rather than assumed complete.

The supported dependency baseline is now Django 5.2.17 LTS with sqlparse 0.6.0. Merged PRs `#71` and `#72` passed the repository dependency audit, Python 3.11/3.12 suites, numerical baseline, Secret Scan and Docker build. This does not by itself establish production security.

## Falsification triggers

A current capability claim must be downgraded if:

- the referenced entry point no longer exists;
- the registered model version does not match runtime output;
- a documented command or route fails under the supported environment;
- numerical identity changes without a version/output-diff decision;
- an artifact hash cannot be verified;
- a current document describes a target model as implemented;
- source-level evidence does not support a patient-derived input;
- a required security/dependency gate is red or covered only by an expired risk decision.

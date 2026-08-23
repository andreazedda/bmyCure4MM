---
title: Project Roadmap
status: CANONICAL_STRATEGY
owner: Andrea Zedda
audience: maintainers, contributors and collaborators
last_verified_at: 2026-08-24
last_verified_git_sha: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
canonical_execution_tracker: https://github.com/andreazedda/bmyCure4MM/issues/31
---

# Project roadmap

GitHub issues and milestones are the execution control plane. This page explains the scientific sequence; issue `#31` is the mutable canonical tracker.

## Current path

```text
COMPLETED:
#11 → #12 → #15 → #19 → #69 → #70

CURRENT:
#14 documentation/source of truth
→ #13 mandatory CI and PR governance
→ #23 minimum QA invariants
→ #26 close M0-R

NEXT:
#32 → #33 → #34 → #35 → #36 → #18
```

## Current dependency and framework baseline

The deterministic lock now uses Django 5.2.17 LTS and sqlparse 0.6.0. Merged PRs `#71` and `#72` removed the current sqlparse findings and all Django 4.2 advisory exceptions while preserving Python 3.11/3.12, numerical and container evidence. The remaining exact development-tool advisory decision is documented in the dependency triage record.

## Milestones

| Milestone | Purpose | Exit emphasis |
|---|---|---|
| M0-R — Canonical Scientific Baseline | One technical and scientific source of truth | Documentation, supported/security-triaged dependencies, protected CI, minimum contracts and privacy invariants |
| M0-S — Safe Shared Research Platform | Safe multi-user/external research operation | Product auth, object authorization, RBAC, security, quotas, runtime, API and observability |
| M1 — Real-Patient Research Loop v0.1 | First source-verified lineage-bound longitudinal loop | Observation, calibration decision, exposure, toxicity attribution, temporal backtest and immutable report |
| M2 — Measurement & Evidence Layer v0.2 | Define what inputs and outputs mean | Evidence registry, biomarker semantics, genomic risk, MRD, imaging, functional outcomes and cards |
| M3 — External Validation & Benchmark v0.3 | Test beyond the private case | Governed datasets, frozen tasks, simple baselines and external evidence |
| M4 — Multiscale MM Model v0.4 | Represent competing disease mechanisms | Clones, residual-disease hypotheses, niche, normal lineage, immunity, bone and organ states |
| M5 — Therapy Design & Optimization v0.5 | Research therapeutic strategies | Resistance, escape, multi-drug PK/PD, modalities, robust Pareto, causal protocol and control |
| M6 — Research v1 — E2 | Auditable reproducible-research release | Complete evidence bundle, release artifacts and explicit E2 intended use |

## Controlled parallelism

Private read-only source verification for M1 issue `#32` may prepare while M0-R closes. It may not promote a new canonical dataset, calibration, what-if result or project claim before protected documentation, CI and QA gates exist.

## Decision rule

Priorities maximize scientific validity, falsifiability, evidence, reproducibility and risk reduction while penalizing implementation cost, migration risk, model complexity, validation burden, non-identifiability and maintenance cost.

## Release boundary

M6 targets `E2_reproducible_research` with `clinical_decision_support: false`. Clinical pilot, SaMD, regulated use, automated treatment selection and prescriptive dosing require a separate roadmap and evidence base.

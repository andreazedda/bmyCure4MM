# Codex v2 — microtask specifications and verified code anchors

**Planning status**: PR-review candidate, not merged into `master` and not executed in production.

## Overview and counts

- 63 goals (54 pre-existing M0–M6 issues + 9 clinician-demo goals).
- 476 currently planned atomic tasks: 54 × 7 = 378 scientific/engineering legacy items, 9 × 8 = 72 principal E1 items and 26 supplementary source-level E1 items.
- 10 machine-readable microtask manifests (all task IDs unique within group) plus five source-anchor maps.
- 63 verified Python source files and 647 top-level or class-level function/class declarations with actual source lines, at code SHA `4a62cf7f1b6f05724faec203f2f3b8847f88868b`.
- Source anchor maps are **point-in-time navigational data**, not proof that every task already has the exact target function. Files labeled `PROPOSED` or not in maps require an existence check and precise line discovery before editing.
- All task states are PLANNED. These manifests do NOT authorize deployments, privileged access, patient imports or clinical interpretation.

## Canonical documents

- `docs/operations/CODEX_MICROTASK_EXECUTION.md`: short prompt and task-execution state machine.
- `docs/operations/CODEX_GOAL_INDEX.yaml`: 63-goal overview, dependency and gates.
- `docs/operations/E1_RELEASE_IMPLEMENTATION_V2.md`: source-line-specific release audit, job lifecycle, localization, UX, security, scientific audit and signoffs.
- `docs/research/M0_M6_SCIENCE_ACCEPTANCE_V2.md`: all 54 historical scientific/software goals mapped to function-level code candidates, measurements, false hypotheses and validation criteria.

## Task manifest groups

| Group | Goals | Tasks | File |
| --- | ---: | ---: | --- |
| E1 demo G00–G08 | 9 | 98 | `microtasks_v2/G_DEMO.json` |
| Canonical M0-R | 3 | 21 | `microtasks_v2/M0_R.json` |
| Shared M0-S, batch A | 7 | 49 | `microtasks_v2/M0_S_A.json` |
| Shared M0-S, batch B | 6 | 42 | `microtasks_v2/M0_S_B.json` |
| M1 longitudinal loop | 7 | 49 | `microtasks_v2/M1.json` |
| M2 evidence/measurement | 7 | 49 | `microtasks_v2/M2.json` |
| M3 external benchmark | 5 | 35 | `microtasks_v2/M3.json` |
| M4 multiscale biology | 7 | 49 | `microtasks_v2/M4.json` |
| M5 therapy/drug discovery | 9 | 63 | `microtasks_v2/M5.json` |
| M6 release/roadmap | 3 | 21 | `microtasks_v2/M6.json` |

## Verified code anchor maps

- `anchors/SCIENCE_CORE_20261110.json` — mechanistic model, observation, PK, counterfactual, calibration, exposure and backtesting: 11 files, 107 symbols.
- `anchors/SCIENCE_STATE_20261110.json` — lineage, uncertainty, registry and research state: 12 files, 156 symbols.
- `anchors/CHEM_DESIGN_20261110.json` — job/view/model chemistry and design prototypes: 14 files, 131 symbols.
- `anchors/ACCESS_CLINIC_UI_20261110.json` — accounts, patient views/API and main UI: 14 files, 199 symbols.
- `anchors/API_DEPLOY_HELP_20261110.json` — simulator/API, help and scripts: 12 files, 54 symbols.

Each anchor-map file records exact file blob SHA, source URL, function/class name and *line at inspected commit*.

## How Codex consumes a task

1. Read the GitHub issue as canonical. Do NOT trust the manifest over a changed issue checkbox or changed master SHA.
2. Resolve dependencies: to target a *synthetic clinician presentation* work G00 first, then G01/G02/G03/G04/G05; G07 and G08 follow verified prerequisites. G06 patient-data scope is a separately gated clinical/legal project.
3. Select first unchecked and not-blocked microtask; read one goal's manifest and only its relevant source-map files. Avoid loading all 476 tasks into context for a single 10-line change.
4. Refresh checked-out `git rev-parse HEAD`; rerun `rg -n` and AST/grep for exact source target. If declared function moved, update anchor evidence rather than blindly editing stale lines.
5. Record before/expected condition; add a synthetic regression and negative/forbidden-case test. Implement a bounded diff, run focused and required quality/numerical/privacy checks.
6. Open a reviewed PR. Keep status PENDING until required CI and any manual review complete. Do not autonomously make grants, invite real people, ingest protected cases or deploy production.
7. Record actual SHA, outputs, evidence provenance, scientific interpretation limits and explicit FAIL/BLOCKED states. Only then mark the corresponding GitHub checkbox complete and sync Asana/Notion/Trello.

## Data → Model → Decision

- DATA: measurement, assay, units, source, missingness, privacy & legal basis.
- MODEL: state equations, observation mapping, priors/identifiability, calibrations, external tests and provenance.
- DECISION: whether this microtask has met its **software research-readiness** gate. A simulation does NOT establish a validated clinical treatment decision.

## Acceptance and audit priorities

A. **Production defects**: Binding Visualizer fails at `pkg_resources`, persistent ChemJob status can be misclassified on retry, Research Glossary may return 500; all must be reproduced with synthetic tests on the real release image.

B. **Scientific contradictions**: possible PK dosage/input concentration unit mismatch, initial-dose double-count, cross-unit biomarker RMSE, calibration readiness inferred from in-sample improvement; these are audit hypotheses until tested and versioned, and must not be silently patched as cosmetic fixes.

C. **Security/privacy**: owner privilege is explicit authenticated operator action; unrelated registered users never see case data or artifacts; no private identifiers in public Git or ordinary SaaS trackers.

D. **Language and UX**: persistent personal settings, full IT/EN surface, corrected stacking and contrast, keyboard/mobility acceptance; measuring five synthetic clinician tasks is usability evidence only.

## Source documentation

Clinical model evidence and named public literature are cataloged in `docs/research/M0_M6_SCIENCE_ACCEPTANCE_V2.md`, including the revision dates of living guidelines. The source IDs cannot be replaced with generic unsourced 'standard of care' wording. Never invent effect sizes/EC50/PK constants from a DOI.

## Readback validations already performed on planning documents

The 10 manifest files collectively contain exactly 476 records spanning 63 goals; IDs were checked unique within each manifest and source URL candidates were classified as verified or requiring inspection. This is documentation consistency, **not** a successful production smoke test, full test suite or scientific validation.

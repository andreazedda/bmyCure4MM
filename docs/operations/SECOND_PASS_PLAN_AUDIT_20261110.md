# Second-pass planning quality audit — bmyCure4MM

Date: 2026-11-10. Scope: public-safe planning data in draft PR #97. This document **does not verify the running production application**, install dependencies, change clinical data or certify a model.

## Audited artifacts

| Check | Result | Method | Limit |
| --- | --- | --- | --- |
| Goal coverage | **63/63** open goals mapped | Enumerated GitHub issues and grouped manifests | Existing closed issues are not fully remapped here |
| Microtask identity coverage | **476/476** IDs exact match, zero missing/extras/duplicate IDs | Compared ten manifest task IDs to all unchecked ID patterns across 63 issue bodies | This is cross-document consistency, not a functional test |
| Scientific goals | 54 pre-existing open goals × 7 = 378 | Existing issue bodies and M0–M6 manifest groups | Scientific hypotheses still unvalidated |
| Clinician demo work | G00–G08: 72 principal + 26 extra code microtasks = 98 | Existing issue bodies and E1 manifest | None marked completed |
| Source navigation | **63** Python files, **647** file/class/function symbols | Read GitHub file contents at SHA `4a62cf7` and extracted verified source-line positions | Not every task has an exact verified target symbol, all anchors stale if HEAD changes |
| Public sensitive-token spot scan | No matches for selected supplied private person-name tokens or email-pattern strings in 5 core documents, 5 source anchor maps and 10 manifest files | Case-insensitive literal/token scan of planned new GitHub artifacts | NOT a comprehensive PHI leak audit; future uploads/comments must undergo independent checks |
| Execution outcome | 0 new application fixes/deploys/clinical data imports | Planning writes only | No actual synthetic end-to-end run was done in this pass |
| PR state | **DRAFT; not merged** | GitHub PR #97 metadata | No new instructions currently on deployed master |

## Verified code-source observations requiring execution

1. `modules/binding_visualizer/binding_visualizer.py:30,328–330` imports and uses deprecated `pkg_resources`; Setuptools documented removal in version 82.0.0. The user-reported failure is consistent with these source lines. Repair not applied in this planning PR.
2. `chemtools/models.py:61–65` derives status from log prefix/output-file existence; `chemtools/tasks.py:122–166` appends logs and may leave stale artifact paths on retries. This creates a credible misclassification scenario; write failing synthetic tests to confirm exact behavior.
3. `chemtools/viewsets.py:8` defines an all-jobs queryset. Effective DRF permissions may block some users; do not claim unauthorized disclosure without adversarial API test. Explicit per-user object scoping is the design acceptance.
4. `simulator/models_simulation.py:59–70,122–132` sets initial concentration as dose divided by Vd but adds subsequent input rate without explicit Vd scaling. Verify dose-rate units and day-0 administration in a controlled numerical experiment. Do not silently change scientific runtime semantics.
5. `twin_engine/observation_model.py:134–179` computes raw/weighted RMSE from heterogeneous biomarker units; scientific comparison must use justified dimensionless/likelihood or assay-aware scales, and missingness.
6. `twin_engine/calibration.py:466–479` checks in-sample RMSE improvement and a minimum observation count. This does NOT establish parameter identifiability, future predictive validity or clinical benefit.
7. `chemtools/templates/chemtools/tools_home.html:62` labels a static 2025 MM treatment algorithm. Versioned clinical content needs review against current living guidance; this is not automatic prescribing.

## Remaining gaps and exact next action

- **Gate R1 — reviewer**: review draft PR #97 for scope, privacy and source-bound accepted claims, run protected docs/CI checks, then merge only if authorized.
- **Gate R2 — G00 forensic**: fetch actual cluster runtime SHA/digest, active routes, role matrix, HTTP 500 and Celery state using synthetic probes. No cross-user/private data screenshots.
- **Gate R3 — G01/M53**: separate software job fix from scientific PK/PD numerical audit; each must have failing tests and independent review.
- **Gate R4 — G02/G03/G05**: introduce testable bilingual preferences, accessibility and object-level authorization; no direct role promotion through public UI.
- **Gate R5 — G04/G07**: evidence-to-claim audit and controlled synthetic release; zero unresolved P0/P1 on demo critical path.
- **Gate R6 — G06**: restricted clinical case behind data-controller/processor, legal/ethical, encryption, access and human source-verification approvals; never infer clinical efficacy/prescribe.
- **Gate R7 — G08**: human specialist usability acceptance, distinct from scientific/external validation.

## Deterministic Codex mode

Task selection uses GitHub checkbox state and dependency; program loads one goal's JSON and only relevant source-map records, not the entire project backlog. The action loop is source verification → failing synthetic fixture → minimal code diff → focused/negative tests → privacy & scientific gates → PR evidence → GitHub checkbox → cross-tool sync. Explicit human-controlled gates stop progression; asynchronous/background execution is not implied.

## Evidence sources used to define research claims

- IMWG 2014 diagnostic criteria: https://pubmed.ncbi.nlm.nih.gov/25439696/
- IMWG 2016 response/MRD: https://pubmed.ncbi.nlm.nih.gov/27511158/
- R2-ISS (2022): https://pubmed.ncbi.nlm.nih.gov/35605179/
- ASCO living MM guideline v2026.1.3: https://ascopubs.org/doi/10.1200/JCO-26-02399
- TRIPOD+AI: https://pubmed.ncbi.nlm.nih.gov/38626948/
- PROBAST+AI: https://www.bmj.com/content/388/bmj-2024-082505
- Antigen escape BCMA/GPRC5D: https://www.nature.com/articles/s41591-023-02491-5
- EU software qualification MDCG 2019-11 rev.1: https://health.ec.europa.eu/medical-devices-sector/new-regulations/guidance-mdcg-endorsed-documents-and-other-guidance_en

These papers and standards support model definitions, benchmark design and limitations; they do not determine missing patient-level parameters or certify the platform.

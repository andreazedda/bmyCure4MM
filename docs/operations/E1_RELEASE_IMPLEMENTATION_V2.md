# E1 clinician demonstration — code-level execution and acceptance v2

**Status:** planning / draft PR only. **Source baseline:** `4a62cf7f1b6f05724faec203f2f3b8847f88868b`. For any proposed change, re-resolve the line in the current checkout and fail closed on drift. No changes to production, patient records or user permissions are authorized by this document.

## Release boundaries and dual approval

**GO_SYNTHETIC_DEMO**: functioning platform on existing cluster with synthetic inputs, both IT/EN, accessible on mobile/desktop, with claims precisely limited to E1 and supported by source/DOI. **GO_PRIVATE_CASE**: distinct authorization/privacy/source-verification gate; real data or specialist invitations NEVER follow automatically from a synthetic demo success. **GO_E2_RESEARCH**: separate M0–M6 evidence requirements; none of the E1 software tests establish prediction validity or causal effects.

## P0 triage order and path

`G00 #88` forensic read-only → `G01 #89` job/exceptions and `G05 #93` security → `G02 #90` language + `G03 #91` responsive navigation + `G04 #92` scientific audit → `G07 #95` staging/production synthetic verification → `G08 #96` expert acceptance. `G06 #94` real-case lane remains BLOCKED until explicit legal, security and source-data signoffs. Each goal has 8 primary and 2–3 supplemental checklist IDs; do not collapse unrelated tasks into one PR.

## G00 — forensic runtime inventory (#88)

| Microtask | Source / probe | Evidence and acceptance |
| --- | --- | --- |
| G00.T01/T02 | `mmportal/urls.py`; `accounts/middleware.py:10–75`; deployed image labels and readiness | Capture source SHA, digest, `kubectl get` service/pod/worker state, dependency-lock hash and 403/404/200 matrix for each role; NEVER dump user data or environment secrets. |
| G00.T03/T04 | `simulator/views.py:19,54,60,108`; `twin_engine/views.py:47,86`; `chemtools/views.py:89,226`; browser synthetic login | Reproduce each user-observed 500 or dead link with the exact URL template, sanitized exception type, request ID and image SHA. A single HTTP 200 is not proof of function. |
| G00.T05/T06 | Django URL resolver + `docs/development/CODE_ATLAS.md` | Define a capability table: reachable, permission denied as designed, bug, stale documentation, hidden navigation, unsupported scientific claim. Do not mark missing just because its navbar item is hidden. |
| G00.T07/T08 | existing Django tests + new `tests/release/` (PROPOSED) | Each P0/P1 failure has a reproducible test or a documented inaccessible third-party dependency; baseline SHA/fixture and issue link retained. Do not release from a dirty working tree. |

## G01 — binding, ChemJob state, downstream HTTP 500 (#89)

1. **Dependency import** (`modules/binding_visualizer/binding_visualizer.py:30,328–330`): replace `pkg_resources` and its working set with `importlib.metadata.version`, catch `PackageNotFoundError`; test empty installation metadata and actual locked container import. Setuptools documents removal in 82.0.0; avoid a setuptools downgrade as permanent solution.
2. **Command runner** (`chemtools/utils.py:95–120,228`): assert bounded `subprocess.run` timeout, exit-code propagation and stdout/stderr capture. Add retry/backoff only for idempotent external transient errors, not parser failures. Do not log raw clinical or user-supplied molecular payloads.
3. **Job state** (`chemtools/models.py:10,61–69`; `chemtools/tasks.py:122–166`): design schema migration adding explicit PENDING/RUNNING/SUCCEEDED/FAILED with start/end times, canonical failure code and attempt identifier; migrate existing rows conservatively as UNKNOWN if artifact validity is uncertain. A previous HTML file cannot turn a failed retry into COMPLETED. A previous log prefix cannot mislabel a successful later retry.
4. **Concurrency** (`chemtools/views.py:21–35,89,174,204`): detect simultaneous retry and original jobs; enforce a transition with DB compare-and-set/locking and attempt ownership; no duplicate queue jobs; subsequent progress polling reflects attempt-specific status.
5. **Artifact integrity** (`chemtools/views.py:226–305`; `chemtools/viewsets.py:8–15`): require authorized job owner and verified complete output; redirect/link must not expose arbitrary media. Existing `ChemJobViewSet.queryset` is unscoped — verify effective global permission class then add owner-filtered queryset and direct-ID negative tests.
6. **Research Glossary** (`twin_engine/views.py:86–87`): create a failing authenticated and unauthorized URL test; inspect `build_research_glossary` and template dependency. Never classify a generic 500 as fixed merely because menu navigation was changed.
7. **Test matrix**: synthetic correct/incorrect PDB and SMILES, absent target, external 429/503/timeout, missing system package, Celery broker down, duplicate submit, retry after failed, retry after successful, success with missing output, artifact access from another user.
8. **Gate**: all selected ChemTools tasks terminally SUCCEEDED or FAILED in bounded time, all denied-user paths return a safe response, no pending action runs indefinitely and no raw traceback shown in user UI; scientific binding/efficacy claims remain explicitly unsupported unless independently validated.

## G02 — global localization and preferences (#90)

1. Inventory **every** user-visible string in templates, Django forms/exceptions, JavaScript, Plotly labels, tool progress/status, help, exports, login/verification emails and navigation. Preserve stable machine-readable error codes separate from localized labels.
2. Verify `mmportal/settings.py` current language middleware and supported locales before proposing new settings; use Django `gettext_lazy`/`gettext`, `LocaleMiddleware` and locale catalogs where compatible, rather than only CSS hiding `t-en/t-it` spans.
3. Add a server-side `UserPreference` model or compatible field on `AccessProfile` (PROPOSED), with `language` in `{it,en}`, `onboarding_seen_version`, `novice_mode` and optional accessibility flags. Migrate without resetting accounts. Anonymous preference stored only with consent-compatible nontracking cookie and deterministic fallback.
4. Move controls out of persistent top strip (`templates/base.html:158–168`) into account Settings, retain a discoverable compact language entry; keep E1 research warning accessible but nonoccluding.
5. Test login/logout/relogin, new device, invalid/expired session, browser reload, back/forward, markdown render, validation error, e-mail and long Italian text; language persists and fallback is coherent.
6. Acceptance: no mixed-language primary UI across designated demo flows, all form errors and notifications translated, stable API codes preserved, documentation references may remain in original publication language with clear labeling.

## G03 — interface information architecture and accessibility (#91)

1. Product navigation target: **Home → Workspace/Simulations → Research (advanced) → Molecular Tools → Knowledge/Help**; explain within each page the relation between *Learn*, technical *Docs* and the existing *Tutorial* before consolidating links or redirects.
2. For each role define `landing_route`, `available_actions`, `denied_actions`, contextual explanations, and back navigation. Avoid exposing disabled clinical features as active controls to demo members.
3. Audit stacking contexts of header, E1 banner (`templates/base.html:122–136`), nav dropdowns (`:29`), Plotly/3Dmol canvases and error popovers. Root cause may be `overflow:hidden`/containing blocks, not only `z-index`.
4. Test screenshots and actual DOM bounds on 320, 360, 390, 768, 1280px, with keyboard-only focus, 200% text scaling, reduced-motion preference and expanded long-content panels. Validate WCAG 2.2 reflow/contrast/focus not merely 'no horizontal overflow'.
5. Change color tokens for all contrast failures (white-on-white, warning over light cards) and ensure hover/focus/active labels remain legible.
6. Human specialist usability sample (synthetic only): find a case, identify observed vs simulated, inspect evidence/model card, run a what-if with declared parameters, explain uncertainty and download report. Measure completion and interpretation errors; do not call this clinical efficacy validation.

## G04 — scientific validity before visualizing treatment claims (#92)

1. **Prescriptive-looking content:** `chemtools/templates/chemtools/tools_home.html:62` contains a static 2025 treatment algorithm. Freeze/remove from E1 landing until reconciled with the versioned [ASCO living guideline v2026.1.3](https://doi.org/10.1200/JCO-26-02399); guideline summaries are context, NOT automated patient instructions.
2. **Units audit:** `simulator/models_simulation.py:59–70,122–132` initializes `C(0)=dose/Vd` but adds schedule input `u(t)` without explicit `/Vd`. Verify `u` contract (mg/day vs concentration/day); if mg/day, flag dimensional mismatch before any scientific model promotion.
3. **Duplicate initial exposure:** if schedule contains dose at day 0 and solver starts at positive `C(0)` from that dose, explicitly test whether the same dose enters twice. Do not modify numerical semantics without version change.
4. **Schedules:** `simulator/pharmaco/registry.py:134–188` uses `dose/window` with distinct continuous/weekly/cycle semantics; integrate each schedule numerically and compare administered mass across equal-horizon schedules to declared dosing contract.
5. **PD interpretation:** `simulator/models_simulation.py:94–117` maps dimensionless Emax to growth/kill RHS; document whether `Emax` represents effect rate (1/day) or dimensionless fraction. `interaction @ effects` is not automatically drug synergy.
6. **Observation model:** `twin_engine/observation_model.py:19–68` fits parameters from baseline assessment and simulated states. Test data leakage, missing FLC absolute values, renal confounding and influence of chosen latent state; raw RMSE across heterogeneous units is not interpretable without a declared scale.
7. **Calibration:** `twin_engine/calibration.py:466–479` uses optimizer/improved in-sample RMSE/point count; do not accept prediction-readiness without frozen out-of-time holdout, identifiability, baseline comparison, calibration and independent data where available.
8. **Evidence registry:** every biomarker, formula and claim points to a source DOI/PMID/version, associated source scope, observed/derived/hypothetical status and counterevidence. Use existing `twin_engine/model_registry.json` as version authority.
9. Gate: claims matrix with SUPPORTED_AS_STATED / NEEDS_SOURCE / MODEL_HEURISTIC / CONTRADICTED / OUT_OF_SCOPE; discrepancies block only their claim/feature, not unrelated educational paths. Never infer causality from a mechanistic scenario.

## G05 — full owner authority with specific specialist grants (#93)

1. Verified owner profile and MFA, never displayed nickname, are the principal for unrestricted administration. `accounts/middleware.py:42–48`, `accounts/permissions.py:9–13`, `simulator/access.py` and `templates/base.html:32,69` are independent checks to reconcile.
2. Define `PermissionGrant(subject,resource_type,opaque_resource_id,capabilities,expires_at,granted_by,revoked_at)` or equivalent (PROPOSED), with explicit RBAC roles and deny-by-default policy.
3. Owner can manage all local research data and accounts, whereas approved specialist sees exactly owner-granted private case records and synthetic research, not all clinic records. Full owner access must be provisioned through authenticated operator controls, not user-side self-assignment.
4. Invitations are owner-initiated to verified private emails, single-use, expiring, MFA-aware and governed. No auto-approval from a public registration request.
5. Run `owner, clinical reviewer with grant, same role no grant, ordinary demo, anonymous, suspended` × `patient list/detail, research, ChemJob, API, media, exports, search, notifications`; cross-case 403/404 and leak-free counts.
6. Test revocation of existing sessions, caching, already queued jobs and downloadable artifacts; invalidate authorization on each request, not only at login.

## G06 — real-case data is an independent delivery lane (#94)

1. Formal assessment of legal roles (controller/processor), lawful basis and special-category data processing; DPIA and storage location needed; no automatic import based solely on being a family member or platform owner.
2. Define privacy architecture: encrypted transit/storage/backups, encryption key rotation, offline recovery, retention/deletion, case-scoped audit, and access revocation. Threat model explicitly includes static assets, HTTP cache, logs, email, Celery, object store, CI.
3. Import only through private provenance-verified workflow `twin_engine/research_dataset_import.py`, `twin_engine/input_contract.py`, `clinic/models.py` and `twin_engine/management/commands/import_research_dataset.py`; no real records in GitHub, Notion, Asana or Trello.
4. Map laboratory measurements with **date, sample, analyte, units, lab/assay method, detection limit, interpretation status, source assertion**; therapy with actual daily exposure, pauses and adverse-event source; imaging with modality and historical lesion status.
5. Specialist-facing viewer must distinctly label `OBSERVED`, `DERIVED`, `HEURISTIC`, `SIMULATED`, `UNAVAILABLE`. Prognosis and therapy-ranking views remain disabled for clinical interpretation without separate validation/regulatory approval.
6. Negative IDOR/authorization tests, backup restoration with restricted case, search/cache/media scan and independent source verification; explicit GO_PRIVATE_CASE requires human signoff. No appointment/invitation message sent as part of technical automation.

## G07 — qualification, release and rollback (#95)

1. Freeze full route/role matrix, language and representative viewport cases; function-dependent tests plus all critical access-control combinations.
2. Preserve existing `uv.lock` scientific dependency graph; run `uv lock --check`, Ruff/mypy, migration drift, full Django tests, secret/PHI scan, numerical golden, SCA and Docker build.
3. Browser E2E: synthetic login (2 users), owner with MFA, ChemTools success/failure, Research glossary, synthetic simulation and report export, English/Italian and mobile sizes.
4. Load tests: bounded Celery concurrency, broker down, external API timeout, disk pressure and read-only data store, without uncontrolled cluster traffic or real patient access.
5. Controlled backup/restore and immutable web/worker/static SHA/digest parity. Record Kubernetes resource budgets and cluster readiness. Never treat CI pass as deployed production acceptance.
6. Deployment, privilege changes and real data import require human authorized gate; release evidence must include actual image digests, timing, route statuses, 500 count, rollback rehearsal and P0/P1 defect count.
7. Human signoff: demo only if zero unresolved P0/P1 for demonstration critical paths; residual P2 explicitly listed, not silently called perfect.

## G08 — specialist walk-through and acceptance (#96)

1. Build 10–15 minute bilingual, no-PHI walkthrough with declared E1 intended use and evidence cards.
2. Five task questions: can the clinician find a synthetic patient, identify a model assumption, distinguish observed from simulated, interpret uncertainty, and export a report?
3. Track objective task completion, critical interpretation errors, navigation misclicks, and time-on-task (without promising thresholds as clinical validation).
4. User feedback becomes a reproducible defect linked to the precise goal; independent human usability acceptance can fail or request rework.

## Executable test command baseline

```bash
uv lock --check
uv run python manage.py check
uv run python manage.py makemigrations --check --dry-run
uv run python manage.py test chemtools.tests accounts.tests simulator.tests twin_engine.tests
uv run python -m scripts.check_numerical_baseline
uv run python -m scripts.check_repository_hygiene
bash scripts/pre_push_research_safety_check.sh
```

These commands are candidates based on the repository's current tooling. Execute with disposable test settings; record actual PASS/FAIL and environment before treating them as evidence.

## Evidence template for every microtask

```yaml
task_id: G01.T02
baseline_sha: REPLACE_WITH_ACTUAL_SHA
code_file: modules/binding_visualizer/binding_visualizer.py
symbol: package_version_lookup
line_range: VERIFY_AT_BASELINE_SHA
before: observed synthetic failure with request ID
after: expected outcome on locked image
tests: [failing_regression, targeted_suite, privacy_checks]
science_status: NOT_APPLICABLE_FOR_IMPORT_FIX
artifact: reviewed_PR
release_gate: pending
state: PLANNED
```

**No symbolic name or line in the template authorizes skipping source verification.** Completion is evidenced only by actual test logs, an accepted PR, and the appropriate runtime or scientific gate.

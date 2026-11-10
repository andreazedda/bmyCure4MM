# Codex goal/microtask execution contract — E1 specialist demonstration and M0–M6 science

Status: **PLANNING / REVIEW REQUIRED**. This document must be merged through protected PR gates before being treated as authoritative local instructions. Source baseline: `4a62cf7f1b6f05724faec203f2f3b8847f88868b`, 2026-11-10.

## Scope

- E1 scientific research, synthetic demonstration, evidence transparency, precision-engineered access control and user-friendliness. It is **not** validated patient-specific prognosis, intervention selection or identified causal therapy effects.
- GitHub issues **#88–#96** are the nine clinician-demo goals, each with code-level checklists, dependency gates and test requirements.
- Existing open scientific/software issues **M0-R/M0-S/M1–M6** remain their canonical goals; the appended section `Code-executable microtask plan — 2026-11-10` specifies seven atomic numbered tasks per issue. Issue #31 remains master scientific roadmap.
- Do NOT duplicate existing issues into a competing roadmap. Do NOT edit scientific equations or numeric thresholds for a UX-only task.
- If a source path is proposed, confirm absence/presence before creating it. All line numbers are commit-dependent. Locate an AST symbol with `rg -n` at the checked-out HEAD; put exact line range and SHA in task evidence.

## One-line invocation

```text
Continua goal G01 (o GitHub #<number>): esegui il primo microtask non completato e non bloccato, prova l'esito, registra evidenze su GitHub e continua nello stesso turno fino al gate o al primo blocco che richiede azione umana.
```

The agent works synchronously within its active execution. Do not promise background execution, do not invent completion. A separate new invocation can resume from the persisted issue checkboxes.

## Goal selector

1. Fetch current `master` and current issue bodies, check state and all `- [ ]` task IDs. Never rely on cached summaries or historical closed milestones as proof.
2. If a goal is explicitly selected, choose its first open, unblocked microtask. Otherwise prioritize `#88 G00` forensic baseline; then `#89 G01`, `#90 G02`, `#91 G03`, `#92 G04`, `#93 G05`. Each may proceed on separate reviewed PRs after dependency checks.
3. `#95 G07` only after G01/G02/G03/G04/G05 relevant gates are passed on synthetic test environment.
4. `#96 G08` synthetic presentation after G07; its private-case lane requires G06 separately.
5. `#94 G06` (real health data) may ONLY proceed past design/synthetic tests after legal basis, privacy/security authorization, clinical data-controller agreement and explicit owner controlled signoff. Do not ingest private medical data from a chat into production or GitHub.
6. M1–M6 science research tasks may be worked on synthetic or duly licensed datasets without falsifying external validation. Do not block pure prototype work on hypothetical features, but do not promote unsupported outputs to production claims.

## Per-microtask source and evidence contract

Write a short task execution record in a GitHub issue comment (no personal information), for example:

```yaml
issue: 89
microtask: G01.T02
baseline_git_sha: "<actual sha>"
source_anchor: "modules/binding_visualizer/binding_visualizer.py::<import block>"
line_range_verified_at_sha: "<Lx-Ly; only after inspection>"
classification: "SOFTWARE_FIX / not a biological model change"
expected_behavior: "Binding executable imports successfully with locked runtime"
observed_failure: "ModuleNotFoundError: pkg_resources"
test_first: "named synthetic regression"
code_diff: "<PR link>"
verification: "TESTED_PASS | TESTED_FAIL | BLOCKED"
privacy_scan: "PASS | FAIL | NOT_RUN"
scientific_reference: "N/A for code import; setuptools deprecation guide"
model_version_change: "none"
review_gate: "CI; prod deployment separate"
```

All seven conditions must be addressed: exact source + SHA, failing behavior or explicit gap, minimal patch, focused+relevant full tests, scientific source when relevant, no-PHI privacy checks, and updated status/PR evidence. A PASS string without actual logs is prohibited.

## Code-level first defects (verified inspection candidates)

### ChemTools crash #89
- In `modules/binding_visualizer/binding_visualizer.py`, replace `import pkg_resources`; the file uses `pkg_resources.get_distribution(...)` / `working_set` for dependency metadata. Official setuptools guidance replaces this with `importlib.metadata`. A small proposed helper is:

```python
from importlib.metadata import PackageNotFoundError, version

def installed_distribution_version(distribution_name: str) -> str | None:
    try:
        return version(distribution_name)
    except PackageNotFoundError:
        return None
```

- Inspect `chemtools/tasks.py::run_binding_viz_job` and `chemtools/views.py::_enqueue` for subprocess exit failure; inspect `chemtools/models.py::ChemJob.status_label` because output-file existence alone is currently used as a completed signal. Introduce explicit persisted status only with migration/testing, if needed; never redefine existing records silently.
- Tests: `chemtools/tests/test_tasks.py`, `test_views.py`, `test_interface.py`, plus actual target-image synthetic job after controlled release. Source: https://setuptools.pypa.io/en/stable/deprecated/pkg_resources.html .

### Account authorization #93 / #8 / #75
- `accounts/middleware.py::ProductAccessMiddleware.process_view` has a special branch for `profile.is_owner and request.user.is_staff`; not all owners are Django staff. Other modules include `accounts/permissions.py`, `simulator/access.py`, `twin_engine/validators.py`, `clinic/views.py`.
- Read actual account flags via authenticated private operator before proposing owner changes. Centralize permissions by principal+action+resource. Owner may receive complete platform access only after verified identity and MFA; no public role escalation, no identification by screen display nickname. Apply negative cross-case IDOR to URLs, API and artifacts.
- Never include private account emails in public issues or CI.

### Localization and responsive UI #90 / #91
- `templates/base.html` contains bilingual `t-en/t-it` spans, a novice switch, and a persistent intended-use banner. Review DOM stacking context with nav/dropdowns, responsive widths and reading order.
- Full page language must cover server validation strings, JS, charts, static docs shell, accessibility attributes, job results and exports, not only navbar text. Persist user preferences with explicit migration and tests.
- Test browser screenshots at 320/360/390/768/1280 px. WCAG 2.2 AA reflow criterion: https://www.w3.org/TR/WCAG22/ .

### Scientific claim audit #92 / M2
- `chemtools/templates/chemtools/tools_home.html` shows a hardcoded "Multiple Myeloma Treatment Algorithm (2025)" with source labels. It is old relative to ASCO living guideline v2026.1.3 (September 2026). Mark as a historical educational exhibit or remove until independently reviewed; do not silently substitute treatment advice.
- `simulator/models_simulation.py::MathematicalModel.simulate`, `simulator/pharmaco/registry.py::resolve`, `twin_engine/model_registry.json`, `twin_engine/observation_model.py::predict_biomarkers`: compare implementation equations/units to docs and model card. Simulated model outputs remain HYPOTHETICAL/HEURISTIC.
- Cite IMWG diagnosis `10.1016/S1470-2045(14)70442-5`; IMWG response/MRD `10.1016/S1470-2045(16)30206-6`; R2-ISS `10.1200/JCO.21.02614`; ASCO MM living guideline `10.1200/JCO-26-02399`; NCI CTCAE v6 current; TRIPOD+AI `10.1136/bmj-2023-078378`; PROBAST+AI `10.1136/bmj-2024-082505`.
- If a calculation is not identifiable, report `NOT_IDENTIFIABLE`; do not infer individual prognosis.

### Molecular & multiscale hypothesis work M4/M5
- `simulator/design/types.py`, `simulator/design/targeting.py`, `simulator/design/evolution.py`, `simulator/design/reporting.py`, `twin_engine/counterfactual.py`.
- Separately version evidence for antigen targets, clonal evolution, candidate chemistry, PK/PD assumptions and observable outcomes. The current toy Pareto selector is not a validated therapy optimizer. Different candidate mechanisms and healthy compartments require explicit risk/uncertainty.
- Reference CoMMpass molecular profiling `10.1038/s41588-024-01853-0`, antigen-escape study `10.1038/s41591-023-02491-5`, PubChem update `10.1093/nar/gkae1059`, RCSB PDB `10.1093/nar/gkae1091` and FDA/ICH M15 June 2026. These are **sources to validate relevance**, not software parameter values.
- Do not assume progenitor-like reservoir causality, candidate binding efficacy, or treatment superiority from a theoretical model.

## Closure states

`TODO`, `IN_PROGRESS`, `VERIFIED`, `BLOCKED`, `CLOSED`. Mark checkbox complete ONLY when code, test, documentation, evidence and source SHA are independently checked.

- `GO_SYNTHETIC_DEMO` requires zero unresolved P0/P1 user-facing defects; functional end-to-end tests on deployed SHA; role/locale/browser grid; transparent E1 claims.
- `GO_PRIVATE_CASE` additionally requires a lawful authorized clinical data pathway, private storage, explicit per-case grants and negative access validation.
- `GO_E2_RESEARCH_RELEASE` requires separate M0–M6 research validation/benchmark evidence, not merely a polished UI.

## Scientific/legal provenance

- As of 2026-11-10: official ASCO Living Guideline, version 2026.1.3; NCI CTCAE v6.0 (2025); GDPR Art 9 and MDR intended-use boundary.
- Official source links: https://ascopubs.org/doi/10.1200/JCO-26-02399 ; https://dctd.cancer.gov/research/ctep-trials/for-sites/adverse-events ; https://eur-lex.europa.eu/eli/reg/2016/679/ojv ; https://eur-lex.europa.eu/eli/reg/2017/745/oj ; https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m15-general-principles-model-informed-drug-development .

## Public privacy rule

No real clinician/patient names, email addresses, medical history, identifiers, private screenshots, raw labs, protected cohort rows, tokens, credentials or source excerpts in this repository, issues, PRs, CI or public docs. Private SaaS trackers also store only operational metadata. Preserve no-PHI payload boundaries across all releases.

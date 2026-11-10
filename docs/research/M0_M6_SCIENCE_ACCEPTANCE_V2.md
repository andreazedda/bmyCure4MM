# M0–M6 scientific code and evidence acceptance contracts v2

Status: **PLANNING / NOT SCIENTIFIC VALIDATION**. This file is a supplement to existing GitHub issues; the issues retain per-task checkboxes. All links, models and synthetic test examples are public-safe and contain no private clinical records.

## Epistemic separation (required in every page and artifact)

1. OBSERVED is a dated measurement with assay/units/source and known error; DERIVED is deterministic transformation and data/algorithm version; INFERRED is model-dependent estimation of an unobserved state; SIMULATED is output conditional on assumptions and input; HYPOTHETICAL is unvalidated mechanism/parameter; PREDICTED requires a frozen task and measured generalization performance; CAUSAL_IDENTIFIED requires a valid identification argument and suitable data.
2. Never infer that a measurement like FLC ratio is a universal tumor-cell count, that a logistic-growth output is a disease-specific validated forecast, or that a drug-docking/similarity score estimates patient efficacy.
3. For any predicted outcome, define target population, index time, horizon, event, censoring, true comparator, admissible observations, calibration, decision threshold and required external benchmark BEFORE fitting. A single time series or synthetic test is not external validation.

## Mathematical model contract

Latent state \(x_t\), observed biomarker \(y_{j,t}\), treatment/exposure \(u_t\), observation function \(h_j\), and measurement error \(\epsilon_{j,t}\):

\[y_{j,t}=h_j(x_t,u_{\leq t},z_t;\phi_j)+\epsilon_{j,t}.\]

Every \(h_j\) must declare its assay/unit/source/version, confounders \(z_t\), limitations and missingness. A model trajectory is conditional on its ODE and parameter assumptions, not a measured count. An outcome assigned to a different patient without a robust estimation/validation design is not a patient prediction.

For pharmacokinetic concentration \(C(t)\) in mg/L and administered mass rate \(u_{\mathrm{mg}}(t)\) in mg/day, a one-compartment IV-input *illustration* is \(\dot C=-kC+u_{\mathrm{mg}}/V_d\), where \(k\) has 1/day and \(V_d\) has L. This is a **dimensional counterexample** to audit the existing solver, not a universal model for oral lenalidomide, biologics or complex PK. Bioavailability, absorption, multi-compartment kinetics and formulation require their own evidence and assumptions.

## M0-R — canonical scientific baseline

| Issue | Source anchor(s) verified at planning SHA | Minimum falsifiable acceptance |
| --- | --- | --- |
| #13 CI | `.github/workflows/dependency-quality.yml`, `scripts/check_numerical_baseline.py` | One intentionally failing numeric/privacy/permission fixture causes mandatory CI FAIL; no release from stale green evidence. |
| #23 invariants | `simulator/models_simulation.py:31`, `twin_engine/run_manifest.py:126,264` | Check deterministic seeded output, finite state, unit contracts, artifact hashes and forbidden-identifier rejection with synthetic values. |
| #26 M0-R tracker | `docs/models/REGISTRY.md`, `twin_engine/model_registry.json` | Every dependent source/code/dependency/QA condition is evidenced at exact Git SHA; fail if a child is unresolved. |

## M0-S — shared but case-isolated platform

| Issue | Source anchor(s) | Acceptance and failure test |
| --- | --- | --- |
| #8 RBAC | `clinic/views.py:48,142,164`; `twin_engine/validators.py::validate_patient_access` | Patient list/detail/aggregate and API restrict object visibility; direct ID enumeration and exports by unrelated identity deny without metadata leakage. |
| #9 TLS/CSP | `mmportal/settings.py`; `templates/base.html` | HTTPS/cookie/HSTS/CSP/proxy headers verified on actual image; Plotly/3Dmol allowed only in evaluated policy; unsafe framing blocked. |
| #10 quotas | `chemtools/views.py:21,89`; `chemtools/tasks.py:122`; `simulator/views_manage.py:316` | Per-role concurrency, timeout, cancellation and queue limits prove bounded external and simulator workload under synthetic stress. |
| #16 services | `simulator/models.py:333`; `clinic/views.py:142`; `twin_engine/simulation_bridge.py:76` | Refactor cannot alter ODE golden or hide object-permission failure; route/DB behavior comparison prior/after. |
| #17 runtime | `Dockerfile`, `deploy/k8s/DEPLOY_K3S.md`, `mmportal/settings.py` | Production runtime images/worker/storage/backup/restore and rollback proven; do not assert PostgreSQL where target uses SQLite/NFS. |
| #20 API | `clinic/api.py:21,35`; `chemtools/viewsets.py:8` | Serializer field allowlist and per-user job queryset; unauthorized object IDs and media URLs do not disclose. |
| #21 observability | `accounts/services.py:50`; `chemtools/tasks.py:62` | Correlate 500 and retries by request/job IDs; redact query strings, inputs, patient identifiers and tokens. |
| #24 governance | `SECURITY.md`, `docs/governance/INTENDED_USE.md` | Model/dataset/code licenses and human review boundaries explicit; no public private-case data. |
| #27 tracker | `docs/architecture/CURRENT_SYSTEM.md` | Each M0-S child backed by current code/deploy evidence, not historical app-registration smoke. |
| #75 owner | `accounts/services.py:197`; `accounts/middleware.py:42` | Verified owner with MFA has complete administration; demo role cannot self-elevate; sessions invalidated on permission changes. |
| #76 registration | `accounts/views.py`, `accounts/services.py:310,359` | Owner-issued invitations, expiry, revocation, case scope and tenant isolation tested; public approvals remain separate. |
| #78 specialist pilot | `accounts/tests/browser_journey.py` | Human researcher completes synthetic task sequence; no clinical records or unreviewed access expansion. |
| #79 manuals | `docs/development/CODE_ATLAS.md` | Independent maintainer re-runs developer workflow and rollback, with source SHA and zero secrets. |

## M1 — source-verified longitudinal research

| Issue | Code anchor | Evidence gate |
| --- | --- | --- |
| #18 integrated loop | `twin_engine/research_dataset_import.py`, `twin_engine/calibration.py:32`, `twin_engine/counterfactual.py:27` | Reproduce data→observation→fit→held-out→hypothetical scenario→report with hashes; `INSUFFICIENT_DATA`/`NOT_IDENTIFIABLE` are allowed. |
| #28 tracker | `twin_engine/run_manifest.py:126` | Parent closes only when #32–#36 verified without real-case identifier exposure. |
| #32 source provenance | `twin_engine/research_dataset_import.py`, `twin_engine/input_contract.py` | Every imported value has immutable source assertion, locator in private store, unit and validation flag; conflicts stop import. |
| #33 observation | `twin_engine/observation_model.py:19,41,134` | Explicit biomarker-specific observation functions; compare M-protein/FLC/hemoglobin predicted vs observed with units and missingness, not latent-cell identity. |
| #34 lenalidomide exposure | `twin_engine/exposure_bridge.py:55,227`; `twin_engine/therapy_schedule.py:31` | Different daily patterns with same mean must yield distinct time-resolved exposure manifests; interruptions must affect actual scheduled days. |
| #35 toxicity attribution | `twin_engine/toxicity_dynamics.py:29,154`; `twin_engine/models.py` | Temporal association, confounders, dechallenge/rechallenge and alternative disease explanation; do not state a drug caused toxicity without identification. |
| #36 calibration | `twin_engine/calibration.py:466–479`, `twin_engine/backtesting.py:30` | Minimum observations and in-sample RMSE cannot certify prediction; frozen future split, comparator, identifiability and uncertainty must pass independently. |

## M2 — measurement validity, response and prognostic staging

| Issue | Code anchor | Literature-specific measurable test |
| --- | --- | --- |
| #22 model/data cards | `twin_engine/model_registry.json`, `twin_engine/run_manifest.py` | Every model has target, endpoint, units, version, calibration/validation scope, censoring/missingness if relevant, and proper evidence status. |
| #29 tracker | `docs/models/REGISTRY.md` | No closed M2 without traceable classification of every biomarker and evidence issue. |
| #37 evidence registry | `twin_engine/model_references.json`, `docs/models/REGISTRY.md` | Each claim resolves DOI/PMID/version/source date; outdated guideline text is historical context not a clinical instruction. |
| #38 biomarkers | `clinic/models.py:80` (Assessment); `twin_engine/observation_model.py:126` | Map kappa/lambda, ratio, M-protein, Hb, renal state and missingness separately; auditable units and assay method. |
| #39 R2-ISS | `clinic/models.py:60,80`; `twin_engine/state_model.py` | Compute staging only with required ISS, LDH, cytogenetic features and their assay dates; compare with published R2-ISS 2022 synthetic cases; not tumor doubling-time. |
| #40 MRD/imaging | `clinic/models.py`; `twin_engine/observation_model.py` | Record assay LOD, flow/NGS method, sample date, imaging modality/lesion state; MRD-negative != zero residual disease. |
| #41 function and organ | `clinic/models_symptoms.py`; `twin_engine/toxicity_dynamics.py` | Use versioned validated measurement instruments where available; no conflation of mobility/pain/organ function with tumor count. |

## M3 — true external generalization

| Issue | Code anchor | Test design and counterexample |
| --- | --- | --- |
| #42 data license | `twin_engine/research_dataset_import.py` | Dataset access/consent/terms verified; no actual clinical cohort in public CI. |
| #43 CoMMpass adapter | `twin_engine/research_dataset_import.py`; `twin_engine/schemas/` | Versioned schema mapping and patient-grouped train/validation partition; leakage attempts rejected. |
| #44 baseline benchmark | `twin_engine/backtesting.py:30`, `twin_engine/uncertainty.py` | Prespecify horizon, outcome, time origin, censoring and baseline; no gain vs persistence/simple model is a valid negative finding. |
| #45 evidence bundle | `twin_engine/run_manifest.py:126,264` | Frozen external evaluation, calibration and subgroup uncertainty, licensing/completeness with model card; reject unsupported `validated` badge. |
| #59 tracker | `docs/product/ROADMAP.md` | External validation can close with documented failure, but cannot confer validated-prediction claims when failed. |

## M4 — competing multiscale biological hypotheses

| Issue | Code anchor | Required falsification |
| --- | --- | --- |
| #46 state vector | `simulator/design/types.py`, `twin_engine/model_registry.json` | Every new component has definition/unit/time scale/source/measurement proxy; unidentifiable states flagged not silently filled. |
| #47 clonal evolution | `simulator/design/evolution.py` | Fractions sum to one and remain nonnegative where model demands; compare with single-population baseline and respect uncertainty. |
| #48 reservoir hypotheses | `simulator/design/types.py::TumorState` | H0 mature cells, H1 quiescent, H2 plasmablast, H3 B-like, H4 niche-protected compete; no universal stem-like reservoir assertion. |
| #49 niche | `simulator/design/simulator.py` | Coupling enabled/disabled comparator and evidence-based parameter ranges; fail if unsupported or non-identifiable. |
| #50 normal lineages | `simulator/design/types.py::NormalCompartment` | Normal plasma cells, stem/progenitor and immune states distinct from malignant states; toxicity outputs do not become infection predictions without data. |
| #51 bone/organs | `clinic/models_symptoms.py` | Historic bone damage != current active disease; organ proxy changes require explicit measurement models. |
| #60 tracker | `docs/architecture/TARGET_VIRTUAL_LAB.md` | M4 remains hypothesis implementation until sufficient observational identifiability and benchmark evidence. |

## M5 — molecule discovery, mechanism, immunotherapy and optimization

| Issue | Code anchor | Evidence boundary |
| --- | --- | --- |
| #52 resistance | `simulator/design/evolution.py`; `simulator/design/targeting.py` | Separate target antigen escape, subclone selection, pathway bypass with source-linked evidence, not one unexplained resistance scalar. |
| #53 multidrug PKPD | `simulator/pharmaco/registry.py:16,134`; `simulator/models_simulation.py:31` | Test units, schedules, interactions and sequence; diagonal amplification is not mechanistic synergy. |
| #54 CAR-T/ADC | `simulator/design/types.py::ADCModel,CARTModel,LogicGate`; `simulator/design/toxicity.py` | Modality-specific assumptions, target antigen binding and normal-cell impacts, without patient-level safety/efficacy claims. |
| #55 Pareto | `simulator/design/reporting.py::_pareto_front,run_design_report` | Return a set of non-dominated hypothetical candidates; any single winner requires explicit weights, robustness and applicability. |
| #56 causal | `twin_engine/causal.py`; `twin_engine/counterfactual.py:27` | Define estimand, DAG, assumptions, identifiable backdoor/longitudinal regime; simulations alone do not identify effects. |
| #57 dependency data | `chemtools/pdb_api_client.py`, `chemtools/utils.py` | Dependency/drug-response assay-specific and MM-domain evidence, not PubChem structural similarity alone. |
| #58 adaptive control | `twin_engine/therapy_schedule.py`; `twin_engine/robustness.py` | Closed-loop synthetic control compared with fixed schedule under uncertainty; cannot prescribe. |
| #61 tracker | `docs/product/ROADMAP.md` | M5 closes only as a research target with all model/evidence contracts; no therapeutic superiority. |
| #77 mechanism→molecule | `chemtools/models.py`, `simulator/pharmaco/registry.py`, `twin_engine/simulation_bridge.py` | Molecule identity→target evidence→assay potency→PK model→Twin scenario must be separate validated transitions; missing EC50 or clearance blocks pipeline. |

## M6 — reproducible research release distinct from synthetic demo

| Issue | Acceptance |
| --- | --- |
| #25 E2 release | Immutable code/data/model/evidence manifest, out-of-sample and external benchmark where claimed, security+privacy QA, deployment digest and model/data cards. |
| #30 tracker | Link all M0–M5 prerequisites and explicitly mark any unachieved gate; do not conflate E1 demo readiness with validated E2 research. |

## Priority validation experiments (all synthetic/public-safe)

| Test ID | Synthetic counterexample | Falsifying observation |
| --- | --- | --- |
| SCI-01 | (V_d=10\,L), mass-rate (u=10\,mg/day), no elimination, zero initial concentration | A resulting concentration slope other than (1\,mg/(L\cdot day)) signals unit mismatch under this declared one-compartment model. |
| SCI-02 | Same time horizon and mass but daily vs alternate-day doses | Identical time-resolved concentration vectors despite different inputs falsify schedule sensitivity. |
| SCI-03 | A calibration model with one observation and several free parameters | A READY status is incorrect; return INSUFFICIENT_DATA or NON_IDENTIFIABLE. |
| SCI-04 | A perfect training fit with deliberately shifted future observation | Labeling future prediction validated is unsupported; held-out loss must be reported. |
| SCI-05 | MRD-negative synthetic assay with LOD (10^{-5}) | Returning zero malignant cells is an invalid inference. |
| SCI-06 | Structural similarity high to known active molecule but assay potency missing | Automatically populating `EC50` as if measured is unsupported. |
| SCI-07 | Antigen detected by one assay but epitope mutation disrupts binding | Claim that detectable antigen guarantees a particular drug engages is contradicted. |
| SCI-08 | Same patient ID and time point appears in training and validation after join | Benchmark must detect and block leakage. |

## Verified scientific references (not direct parameter lookup tables)

- IMWG diagnosis (Rajkumar et al., Lancet Oncology 2014): https://pubmed.ncbi.nlm.nih.gov/25439696/ ; DOI 10.1016/S1470-2045(14)70442-5.
- IMWG response and MRD (Kumar et al., Lancet Oncology 2016): https://pubmed.ncbi.nlm.nih.gov/27511158/ ; DOI 10.1016/S1470-2045(16)30206-6.
- R2-ISS (D'Agostino et al., JCO 2022): https://pubmed.ncbi.nlm.nih.gov/35605179/ ; DOI 10.1200/JCO.21.02614.
- ASCO living guidance MM v2026.1.3 (last updated 2026-09-14, article 2026-09-30): https://ascopubs.org/doi/10.1200/JCO-26-02399 . Version and review date mandatory; guidance is not patient-specific automated prescription.
- TRIPOD+AI (BMJ 2024): https://pubmed.ncbi.nlm.nih.gov/38626948/ ; DOI 10.1136/bmj-2023-078378.
- PROBAST+AI (BMJ 2025): https://www.bmj.com/content/388/bmj-2024-082505 ; DOI 10.1136/bmj-2024-082505.
- CoMMpass molecular profiling (Nature Genetics 2024, correction 2025): https://www.nature.com/articles/s41588-024-01853-0 ; evidence for genomic heterogeneity, not a validated individual predictor.
- Antigen escape mechanisms (Nature Medicine 2023): https://www.nature.com/articles/s41591-023-02491-5 ; mechanism may be treatment-context specific.
- EU MDR software classification guidance MDCG 2019-11 rev.1 (June 2025): https://health.ec.europa.eu/medical-devices-sector/new-regulations/guidance-mdcg-endorsed-documents-and-other-guidance_en .

## Research acceptance rule

Promotion to a stronger claim requires (i) measurable target and sources, (ii) declared model assumptions, (iii) identifiable estimation and uncertainty where applicable, (iv) frozen independent comparator/testing, (v) explicit domain restriction, (vi) versioned artifacts and code, (vii) privacy/legal authorization. If any required item is absent, the output must preserve its actual lower-status label (e.g. HYPOTHETICAL/NOT_IDENTIFIABLE/EXTERNAL_VALIDATION_FAILED). Do not invent outcomes to close a goal.

---
title: Developer Audit Verification — 2026-11-06
language: en
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: maintainers and reviewers
last_verified_at: 2026-11-06
last_verified_git_sha: 313e2e4d8a839b9c14a8d853d1f008071a6cb013
---

# Developer audit verification — 2026-11-06

This record covers the source/documentation audit and executable developer exercises performed during the follow-up UX audit. It is distinct from the earlier access-release acceptance and from independent human review of #79. The application remains E1/research-only. No real-patient data, external mailbox or production secret was needed for these exercises.

## Source, environment and evidence boundary

- Audited baseline: `fa9c35473302debf3c9a289d5f34f0bc7f771dbd`; source anchors in the inventory are pinned to that commit.
- Original corrections and exercises were run in the audit working tree based on that commit. Final protected PR #86 head is `49d1a89b12cd0508e8c6aba2340cb193a31eed11`, merged source `313e2e4d8a839b9c14a8d853d1f008071a6cb013`; immutable inventory source anchors deliberately retain the audited baseline.
- Local environment: macOS ARM64, Python 3.11, frozen project dependencies, disposable SQLite, independent synthetic fixtures, null test logging and local-memory email.
- `UV_NO_SYNC=1` was used while executing `uv run` against the already synchronized environment so concurrent verification would not change shared dependencies. The clean-environment procedure explicitly performs `uv sync --frozen` first.
- The application database, patient media, owner password/MFA and email capabilities were not copied or printed.

## Confirmed documentation defects and corrections

| ID | Severity | Finding from source/manual audit | Correction / evidence |
| --- | --- | --- | --- |
| DX01 | High | Handbook offered a blank database with no explanation of how a visitor gets a usable learning example or a test activation link | Documented the idempotent `ensure_learning_example` command and the explicit isolated browser/mail-sink journey; no default password |
| DX02 | High | Troubleshooting advised sharing complete logs | Replaced with minimal sanitized state, command/version/exit evidence; explicit credential/request/patient exclusions |
| DX03 | High | Troubleshooting claimed migrations fix a Python `AttributeError`, runserver starts Celery, and a shell eager variable is a supported production switch | Replaced with source/image/import diagnosis and the actual Python test-settings behavior |
| DX04 | Medium | Endpoint reference listed nonexistent `/simulator/scenario/<id>/` and `/chem/jobs/<id>/`, and vague POST/export implementations | Rebuilt the route/authorization map from URLConfs; verified three concrete `resolve()` outputs without querying records |
| DX05 | High | Results documentation could be read as permission to expose MEDIA_URL artifacts | Clarified stored locators vs authorized export/embed routes, ownership, current access checks and closed direct media |
| DX06 | Medium | Current architecture cited old source and said protected aggregate CI was unverified despite recorded release gates | Updated baseline and concrete `required` dependency graph; retained unverified infrastructure/governance limits |
| DX07 | Medium | Broad package table did not provide critical signatures, model fields, failure modes or source anchors | Added 31-module/94-symbol AST inventory, callable/side-effect/test map, typed model fields and static drift checker |
| DX08 | Medium | Short tumor equation omitted implementation details useful for debugging PK/PD, interactions and solver sampling | Added exact current RHS, units/assumptions, clipping, tolerances and output-grid distinction; preserved scientific code |
| DX09 | Medium | Troubleshooting lacked an executable expected-failure and recovery exercise | Extracted and ran the documented invalid-horizon failure, analytic zero-effect check and synthetic SQLite backup/restore |
| DX10 | Medium | Developer commands could inherit a non-disposable database override or production redirect settings | Explicit local debug/HTTP settings, `unset BMY_TEST_DATABASE_NAME`, dedicated temporary root and no production credentials |

## Executed verification

These counts are the observed command results for this audit, not the total application test count. Runtime tests and browser evidence are recorded separately in the [UX audit](../operations/UX_AUDIT_2026_11_06.md) and [access release ledger](../operations/AUTHENTICATION_AND_ACCESS.md).

| Check actually executed | Result | Measured outcome / limits |
| --- | --- | --- |
| `python -m scripts.check_developer_atlas` | PASS | 31 modules / 94 top-level symbols, method signatures and field types resolve; all 94 immutable baseline line anchors independently checked; no DB or network |
| Disposable `migrate`, `ensure_learning_example` twice | PASS | First call `INITIALIZED_SYNTHETIC`; second `PRESERVED`; no duplicate example |
| `accounts.tests.test_lifecycle`, `test_concurrency`, private-ticket relogin test and real-engine/authorized-artifact test | PASS | 24 tests in 1.519s; actual local password/session requests and existing numerical kernel, isolated email backend |
| `simulator.tests.test_learning_entry`, `test_parameter_validation`, `twin_engine.tests.test_run_manifests`, `test_exposure_bridge` | PASS | 21 tests in 0.389s; synthetic fixtures, no production data |
| Two Python snippets extracted from the troubleshooting page and executed verbatim | PASS | `EXPECTED_INVALID_HORIZON=PASS`, `RECOVERY=PASS`; logistic endpoint 1020.891279; synthetic backup/restore `quick_check=ok`; total 5.684s |
| Three URL resolutions from endpoint reference | PASS | `simulator:scenario_detail`, `simulator:attempt_export`, `chemtools:job_detail`; route existence, not authorization proof |
| Ruff check/format on the atlas checker | PASS | Static validator conforms to repository rules |
| `mkdocs build --strict` | PASS | Build checks passed in 2.14–3.21s, most recent 2.19s; all new atlas/manual/evidence routes included in navigation |
| `scripts.check_repository_hygiene` and `git diff --check` | PASS | No new hygiene/whitespace failures; two pre-existing invalid notebooks remain explicitly recorded by hygiene tool |

The exact walkthrough labels and inline failure/restore commands are maintained in [manual diagnosis](TROUBLESHOOTING.md). Tests used dedicated `/tmp/bmy-doc-audit-*` roots; path names are examples, not required shared directories. Reproduce with a fresh `mktemp` root and retain only sanitized outcomes.

## Final PR #86 and target-image verification

This section adds PR #86 evidence without rewriting earlier focused test
measurements. The browser cache finding below prevents treating that
intermediate deployment as completed UX acceptance.

| Check / identity | Observed outcome |
| --- | --- |
| Protected PR | [#86](https://github.com/andreazedda/bmyCure4MM/pull/86); tested head `49d1a89b12cd0508e8c6aba2340cb193a31eed11`; merged source `313e2e4d8a839b9c14a8d853d1f008071a6cb013` |
| Whole local suite | PASS: 473 discovered Django tests on disposable SQLite, macOS ARM64 / Python 3.11 |
| Combined account/docs browser suites | PASS: two tests in 38.646s; real local HTTP forms, isolated synthetic users and SMTP sink, desktop/mobile views |
| Frozen numerical baseline | PASS: zero differences; no golden-output update |
| Protected and publication gates | PASS: required CI, post-merge checks, GitHub Pages and GHCR |
| Published image | Index digest `sha256:0d707ef6c20ee3f6cc4650730610ee045af56955130ab28a5bb98b7bdccfc405`; OCI source label and image environment match merged source |
| Actual target disposable test Job | `ux-verify-313e2e4` PASS: 83 tests (34 accounts, 46 docs, three learning-entry tests), 38.559s tests / 96.69s runner, peak RSS 269340 KiB; existing NFS PVC and disposable SQLite, unchanged 1 CPU / 512 MiB ceiling |
| Target SMTP/TLS connectivity | PASS for connectivity/TLS only; no new authenticated send or received message claimed |
| Controlled backup/bootstrap/catalogue/static Job | `accounts-release-313e2e4` PASS: matching source, backup restore `quick_check=ok`, ACTIVE owner/MFA and open request/activation gates preserved, static assets present, empty catalogue initialized with one active synthetic scenario |
| Web/worker runtime identity and route probes | PASS: web and both workers report PR #86 source/digest, 1 CPU / 1 GiB ceilings; owner/MFA/gates/catalogue preserved; 20 public route/asset probes passed |
| Returning-browser static-cache gate | Historical FAIL for intermediate PR #86: unversioned CSS remained cached up to 14400 seconds despite correct collected files; PR #87 adds source-versioned asset URLs |
| Corrected-image acceptance contract | Require protected checks, exact published source/digest, controlled rollout/readback and returning-browser asset verification; record the result and asset hashes in the dated release-evidence comment in [GitHub #79](https://github.com/andreazedda/bmyCure4MM/issues/79), linked to PR #87 |

The dated #79 release-evidence comment is the final operational ledger for the
current deployed version and cache outcome. This source snapshot deliberately
preserves the intermediate PR #86 measurements and operator acceptance contract;
it does not invent a future build SHA or assert a PR #87 merge/deployment before
verification. The historical cache failure is not a permanent statement about
runtime after a subsequent verified release.

Final review also reproduced an internal symlink allowlist bypass with a
synthetic temporary file, then verified its rejection: documentation page/raw
returned 404 and search omitted the link. The 46 documentation tests passed.
Registered-user clinical assessment choices are now excluded before lookup,
and forged Twin context is rejected across the three simulator aliases;
nonstaff owner behavior is covered separately. Revoked accounts cannot be moved
to suspended and renewed. These safeguards are bounded to the existing member
entry path and do not close broad legacy-clinic authorization issue #8.

The controlled release initializes only a wholly empty existing `Scenario`
catalogue using `ensure_learning_example`; reruns preserve all existing and
inactive examples. It does not create another application, patient fixtures,
per-user provisioning or a computation worker. Rollback to the previous image
requires suspending external members and closing registration/activation gates
before reverting security fixes. Preserve schema/PVCs and use a separate
old-image collectstatic Job, because the previous image lacks this command.

## Available checks and acceptance that remain distinct

| Check / acceptance | How to reproduce or obtain it | Current claim |
| --- | --- | --- |
| Full Django suite and migration drift | Handbook verification command matrix | Final whole-tree result above and PR #86; original focused exercises alone are not whole-tree acceptance |
| Real browser forms/CSRF/QR/two users/results/tickets | Explicit `accounts.tests.browser_journey` command with SMTP sink | Browser audit evidence is separate from the local-memory walkthrough |
| Frozen scientific regression baseline | `python -m scripts.check_numerical_baseline` | Preserves governed outputs; not efficacy/calibration validation |
| Strict docs build | `uv run --frozen --no-default-groups --group docs mkdocs build --strict` | Required before merging; a successful build does not prove source statements or UI rendering |
| Python 3.12, frozen image, secret/dependency gates | Protected `required` GitHub workflow and companion docs/secret workflows | Use exact PR head and conclusion; local Python 3.11 is not a substitute |
| Production source/digest/assets/mail/rollout | Governed existing-cluster release runbook | Intermediate PR #86 readback passed; current version and corrected-image cache outcome belong to the dated #79 release-evidence comment, not an anticipated hash in this source snapshot |
| Independent maintainer no-AI exercise | Human follows two walkthroughs, expected failure and recovery, records own result | PENDING_USER_ACTION; no automated run can certify another person's understanding |
| First real researcher acceptance (#78) | Owner explicitly authorizes a person and collects actual feedback | Separate gate; no synthetic user stands in for it |

## Reviewer exercise

A maintainer should be able to start from a clean clone, identify the current image versus checkout identity, find each view/service/model, execute the two walkthroughs, reproduce the invalid-horizon rejection and explain why rollback preserves additive schema and requires token quarantine after an actual restore. Record the exact commit, commands and outcome without an AI session or secret dump. File a concrete mismatch against the named source/document, preserving prior release evidence; do not close #79 solely because this document or the CI exists.

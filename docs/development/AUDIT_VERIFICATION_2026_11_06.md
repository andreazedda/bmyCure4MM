---
title: Developer Audit Verification — 2026-11-06
language: en
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: maintainers and reviewers
last_verified_at: 2026-11-06
last_verified_git_sha: fa9c35473302debf3c9a289d5f34f0bc7f771dbd
---

# Developer audit verification — 2026-11-06

This record covers the source/documentation audit and executable developer exercises performed during the follow-up UX audit. It is distinct from the earlier access-release acceptance and from independent human review of #79. The application remains E1/research-only. No real-patient data, external mailbox or production secret was needed for these exercises.

## Source, environment and evidence boundary

- Audited baseline: `fa9c35473302debf3c9a289d5f34f0bc7f771dbd`; source anchors in the inventory are pinned to that commit.
- Corrections and exercises were run in the audit working tree based on that commit. Final protected-CI and release SHAs belong to the dated PR/release ledger; this page does not invent a hash for uncommitted edits.
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

## Available checks and acceptance that remain distinct

| Check / acceptance | How to reproduce or obtain it | Current claim |
| --- | --- | --- |
| Full Django suite and migration drift | Handbook verification command matrix | Final audit PR records the executed whole-tree result; do not replace it with the 45 focused tests above |
| Real browser forms/CSRF/QR/two users/results/tickets | Explicit `accounts.tests.browser_journey` command with SMTP sink | Browser audit evidence is separate from the local-memory walkthrough |
| Frozen scientific regression baseline | `python -m scripts.check_numerical_baseline` | Preserves governed outputs; not efficacy/calibration validation |
| Strict docs build | `uv run --frozen --no-default-groups --group docs mkdocs build --strict` | Required before merging; a successful build does not prove source statements or UI rendering |
| Python 3.12, frozen image, secret/dependency gates | Protected `required` GitHub workflow and companion docs/secret workflows | Use exact PR head and conclusion; local Python 3.11 is not a substitute |
| Production source/digest/assets/mail/rollout | Governed existing-cluster release runbook | Not inferred from a local test or healthy HTML response |
| Independent maintainer no-AI exercise | Human follows two walkthroughs, expected failure and recovery, records own result | PENDING_USER_ACTION; no automated run can certify another person's understanding |
| First real researcher acceptance (#78) | Owner explicitly authorizes a person and collects actual feedback | Separate gate; no synthetic user stands in for it |

## Reviewer exercise

A maintainer should be able to start from a clean clone, identify the current image versus checkout identity, find each view/service/model, execute the two walkthroughs, reproduce the invalid-horizon rejection and explain why rollback preserves additive schema and requires token quarantine after an actual restore. Record the exact commit, commands and outcome without an AI session or secret dump. File a concrete mismatch against the named source/document, preserving prior release evidence; do not close #79 solely because this document or the CI exists.

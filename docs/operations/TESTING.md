---
title: Testing and Verification
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: contributors and reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
---

# Testing and verification

## Epistemic boundary

```text
software verification
≠ model verification
≠ calibration evidence
≠ external validation
≠ clinical utility
```

A green test suite proves only the tested contracts under the tested environment.

## Canonical current commands

```bash
uv lock --check
uv sync --frozen --extra chemistry
uv run python manage.py check
uv run python manage.py makemigrations --check --dry-run
uv run python manage.py test
uv run ruff check .
uv run ruff format --check .
uv run mypy
uv run python -m scripts.check_repository_hygiene
uv run python -m scripts.check_numerical_baseline
uv run python -m scripts.audit_dependencies
bash scripts/pre_push_research_safety_check.sh
uv run --frozen --no-default-groups --group docs mkdocs build --strict
```

Use the dependency groups documented in `DEPENDENCIES.md`; core-only and optional-chemistry environments have different purposes.

## Current GitHub Actions

Current workflows cover dependency/reproducibility checks, strict documentation validation, secret scanning, Django 5.2 compatibility, image workflows, UI screenshots and redeployment.

Merged PR `#71` added the stable protected `required` result over locked quality, the Python 3.11/3.12 core matrix and Docker build. This repaired a verified mismatch between branch protection and workflow status names.

Issue `#13` remains open because the final CI contract still requires:

- immutable commit-SHA pinning for external Actions;
- exact base/head pull-request safety mode;
- documentation status/path checks;
- PR and issue templates;
- CODEOWNERS;
- privacy-safe artifact retention;
- final workflow rationalization and branch-protection re-read.

## Current recorded baseline evidence

Merged PRs `#71` and `#72` recorded:

```text
Django 5.2.17 and sqlparse 0.6.0 identity: PASS
uv lock --check and frozen sync: PASS
dependency audit: PASS
Ruff, format and mypy: PASS
repository hygiene: PASS
ordinary and synthetic deploy checks: PASS
migration drift: zero
disposable migrations: PASS
Python 3.11 core check/tests: PASS
Python 3.12 core check/tests: PASS
authorization-sensitive compatibility subset: PASS
numerical baseline: PASS / unchanged
Secret Scan: PASS
Docker locked build: PASS
protected required result: PASS
```

Documentation PR `#68` separately established strict MkDocs validation on pull requests. It must be revalidated on the current master baseline before merge.

## Research safety gate

The current local safety script is designed around staged changes. Until issue `#13` adds a verified base/head CI mode, do not assume that merely invoking the script in a pull-request runner checks the exact proposed diff.

## Numerical baseline

The numerical baseline detects unreviewed output drift. A baseline update requires:

- an output-diff report;
- a model/version decision;
- an invalidation decision for prior runs;
- updated documentation when endpoint interpretation changes.

The dependency and Django 5.2 migrations produced no numerical differences.

## Dependency security

- canonical runtime: Django 5.2.17 LTS and sqlparse 0.6.0;
- no Django advisory exception remains in the audit wrapper;
- one exact development-only pytest advisory decision remains documented;
- high findings cannot be globally ignored;
- a green dependency audit does not by itself establish production security.

## Test data

Tests and CI use synthetic or explicitly governed demo data only. Never copy private patient payloads, source documents or direct identifiers into fixtures, logs, snapshots or uploaded artifacts.

## M0-R minimum QA

Issue `#23` Phase A covers:

- dataset, input, Twin-lineage, run-manifest and artifact contracts;
- canonical hashing and semantic idempotency;
- schedule identity and deterministic numerical invariants;
- current high-risk migration preservation;
- core privacy, authorization and forbidden-claim regressions.

PostgreSQL, extensive property testing, browser E2E, coverage ratchets and asynchronous lifecycle testing remain extended scope unless explicitly completed or transferred to linked issues.

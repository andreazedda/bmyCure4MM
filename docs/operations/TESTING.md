---
title: Testing and Verification
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: contributors and reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
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

Current workflows cover dependency/reproducibility checks, docs deployment, secret scanning, image workflows, UI screenshots and redeployment. Issue `#13` will rationalize them into one stable protected aggregate gate and add exact pull-request safety behavior.

Documentation PR `#68` added strict MkDocs validation on pull requests without deploying PR content. This is a bounded prerequisite for `#13`, not the final aggregate CI design.

## Current recorded PR evidence

At PR `#68` head `3bc54788...`:

```text
strict MkDocs: PASS
secret scan: PASS
Ruff: PASS
format: PASS
mypy: PASS
repository hygiene: PASS
numerical baseline: PASS / unchanged
Python 3.11 core check/tests: PASS
Python 3.12 core check/tests: PASS
Docker locked build: PASS
dependency audit: FAIL
```

The audit failure is not caused by the documentation diff. It identifies newly disclosed advisories in the unchanged canonical dependency graph and is governed by `#69` and `#70`.

## Research safety gate

The current local safety script is designed around staged changes. Until issue `#13` adds a verified base/head CI mode, do not assume that merely invoking the script in a pull-request runner checks the exact proposed diff.

## Numerical baseline

The numerical baseline detects unreviewed output drift. A baseline update requires:

- an output-diff report;
- a model/version decision;
- an invalidation decision for prior runs;
- updated documentation when endpoint interpretation changes.

The current documentation PR produced no numerical differences.

## Dependency security

- `#69` owns immediate sqlparse remediation and exact current advisory triage.
- `#70` owns migration from unsupported Django 4.2 to Django 5.2 LTS.
- a temporary advisory exception must be exact, evidence-backed, owned and expiring;
- high findings cannot be globally ignored;
- shared/production promotion remains prohibited on an unsupported or untriaged framework baseline.

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

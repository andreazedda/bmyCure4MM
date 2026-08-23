---
title: Dependency Audit Triage
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: maintainers and security reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Dependency audit triage

## Issue #70 migration candidate

```text
DJANGO_TARGET = 5.2.17 LTS
SQLPARSE = 0.6.0
DJANGO_ADVISORY_EXCEPTIONS = zero
TEST_TOOL_EXCEPTION = GHSA-6w46-j5rx-g56g only
SHARED_OR_PRODUCTION_PROMOTION = pending full #70 validation
```

This branch migrates from unsupported Django 4.2.30 to Django 5.2.17 LTS and removes all Django 4.2 advisory exceptions from `scripts/audit_dependencies.py`.

The candidate is not canonical until its dedicated pull request passes the full compatibility, migration, authorization, numerical and container evidence required by #70.

## Resolved immediate findings

### sqlparse

Issue #69 upgraded the direct and locked security floor to `sqlparse 0.6.0`, addressing:

```text
PYSEC-2026-3696 / CVE-2026-59894
PYSEC-2026-3697 / CVE-2026-71491
PYSEC-2026-3698 / CVE-2026-59893
PYSEC-2026-3699 / CVE-2026-54284
```

### Django

The previous candidate used an exact temporary exception for:

```text
PYSEC-2026-3717 / CVE-2026-15830
```

because the affected GeoDjango path was not found in repository code. That exception was never a supported-framework claim and carried an expiry of 2026-09-30.

The #70 candidate removes the exception entirely by moving to Django 5.2.17 LTS, the patched supported line identified for the advisory at implementation time.

## Remaining exact development-tool decision

```text
GHSA-6w46-j5rx-g56g
```

This pytest temporary-directory advisory requires a malicious local user on a shared UNIX host. The compatible `pytest-playwright==0.7.1` release still constrains pytest below the patched major line. Tests run on isolated CI runners or single-user development machines. This exception remains exact and must be removed when the toolchain supports the patched pytest line.

## Required #70 evidence

```text
Django 5.2.17 and sqlparse 0.6.0 in pyproject and uv.lock
uv lock --check
pip-audit through scripts.audit_dependencies
Django check
synthetic manage.py check --deploy
migration drift
migration graph and disposable forward evidence
Python 3.11 and 3.12 full core suites
authentication and authorization regressions
numerical baseline unchanged
run/artifact/comparability contracts
optional chemistry environment
strict MkDocs
secret and PHI safety
Docker locked build
```

## Historical decisions

Django 4.2.30 was the final 4.2 release and reached end of extended support on 7 April 2026. It is not an acceptable durable M0-R framework baseline.

The earlier exact Django exceptions remain visible in Git history and issue #69 evidence; they are deliberately absent from the #70 audit script.

## Required record fields

Every future exception or fix records:

```text
advisory ID
package/version
severity
applicability status
technical rationale
evidence
owner
expiry/review date
fix issue/PR/commit
```

## Current execution

```text
#69 immediate sqlparse/advisory remediation
→ #70 Django 5.2 LTS migration candidate
→ #13 protected aggregate CI
```

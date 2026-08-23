---
title: Dependency Audit Triage
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: maintainers and security reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Dependency audit triage

## Current decision — issue #69

```text
LOCK_REFRESH_REQUIRED = true
SQLPARSE_TARGET = 0.6.0
DJANGO_42_EXCEPTION = exact_and_temporary
DJANGO_42_EXCEPTION_EXPIRY = 2026-09-30
DURABLE_FRAMEWORK_TARGET = Django 5.2 LTS under issue #70
SHARED_OR_PRODUCTION_PROMOTION = prohibited
```

The deterministic lock at canonical `master` was not changed by documentation PR #68. New advisory data nevertheless made the repository-native audit fail. This issue-69 branch raises the explicit sqlparse floor to 0.6.0 and adds an exact, expiring Django 4.2 exception while the durable 5.2 LTS migration is executed under #70.

## Current findings

### Django

| Package | Advisory | Severity | Affected component | Decision |
|---|---|---:|---|---|
| Django 4.2.30 | `PYSEC-2026-3717` / `CVE-2026-15830` | medium | GeoDjango `GEOSGeometry` / spatial geometry parsing | Temporary exact exception only after no-GIS evidence; expires 2026-09-30. Migrate under #70. |

Initial repository search found no direct `django.contrib.gis`, `GEOSGeometry` or GIS-field use. A clean-checkout scan and settings/runtime review remain part of the PR evidence. Absence of a current code path does not make Django 4.2 supported.

Django 4.2.30 was the final 4.2 release and reached end of extended support on 7 April 2026. The durable target is a current Django 5.2 LTS security patch.

### sqlparse

| Package | Advisory | Severity | Patched version | Decision |
|---|---|---:|---:|---|
| sqlparse 0.5.4 | `PYSEC-2026-3696` / `CVE-2026-59894` | moderate | 0.6.0 | Upgrade required |
| sqlparse 0.5.4 | `PYSEC-2026-3697` / `CVE-2026-71491` | high | 0.6.0 | Upgrade required |
| sqlparse 0.5.4 | `PYSEC-2026-3698` / `CVE-2026-59893` | high | 0.6.0 | Upgrade required |
| sqlparse 0.5.4 | `PYSEC-2026-3699` / `CVE-2026-54284` | high | 0.6.0 | Upgrade required |

Repository search finds no direct sqlparse call outside dependency/lock declarations. Django and tooling may still invoke parsing paths; absence of direct calls is not a risk-acceptance basis.

## Temporary Django 4.2 exception

`scripts/audit_dependencies.py` enumerates exact advisory IDs and refuses to run them after 2026-09-30.

The exception means only:

```text
current affected GeoDjango path is not known to be exposed
```

It does not mean:

```text
Django 4.2 is supported
production security is certified
shared deployment is approved
issue #70 may be deferred indefinitely
```

The exception is removed when #70 merges.

## Previous 13 August triage

The earlier locked baseline removed 155 of 161 broad-environment findings and narrowly triaged six then-current Django advisories plus one development-only pytest advisory:

```text
GHSA-923m-gv2p-w5qp
GHSA-h7pc-vwp9-298g
GHSA-8cjm-8mp7-r2xf
GHSA-3h9f-r86x-qvjx
GHSA-crhf-3pfg-w68w
GHSA-8qcx-xf44-272x
GHSA-6w46-j5rx-g56g
```

Those are retained as exact temporary decisions; they do not cover newly published IDs and now share the hard Django 4.2 expiry.

## Required evidence before #69 can close

```text
pyproject and uv.lock diff
uv lock --check
pip-audit before/after
no-GIS code/settings scan
Python 3.11 and 3.12 core tests
Django check and migration drift
numerical baseline unchanged
Docker build
exact exception/expiry test
linked #70 migration plan
```

## Required record fields

Every exception or fix records:

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
#69 sqlparse 0.6.0 lock refresh and immediate exact advisory decision
→ #70 Django 5.2 LTS migration and removal of Django 4.2 exceptions
→ #13 protected aggregate CI
```

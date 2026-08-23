---
title: Dependency Audit Triage
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: maintainers and security reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Dependency audit triage

## Current status — 24 August 2026

```text
LOCK_DETERMINISTIC = true
DEPENDENCY_AUDIT_GREEN = false
UNTRIAGED_HIGH_FINDINGS = true
SHARED_OR_PRODUCTION_PROMOTION = prohibited
```

The canonical lock was not changed by documentation PR `#68`. Current advisory data nevertheless makes the repository-native audit fail. This is a newly current supply-chain condition.

## Current findings

### Django

| Package | Advisory | Severity | Affected component | Current decision |
|---|---|---:|---|---|
| Django 4.2.30 | `PYSEC-2026-3717` / `CVE-2026-15830` | medium | GeoDjango `GEOSGeometry` and spatial geometry parsing | Applicability must be reproduced. Django 4.2 is unsupported; migrate under `#70`. |

Initial repository code search found no direct `django.contrib.gis`, `GEOSGeometry` or GIS-field use. This is not sufficient for a permanent exception. A clean-checkout scan and runtime/settings review are required.

Django 4.2.30 reached end of extended support on 7 April 2026. The durable target is a current Django 5.2 LTS security patch under issue `#70`.

### sqlparse

| Package | Advisory | Severity | Patched version | Current decision |
|---|---|---:|---:|---|
| sqlparse 0.5.4 | `PYSEC-2026-3696` / `CVE-2026-59894` | moderate | 0.6.0 | Upgrade required |
| sqlparse 0.5.4 | `PYSEC-2026-3697` / `CVE-2026-71491` | high | 0.6.0 | Upgrade required |
| sqlparse 0.5.4 | `PYSEC-2026-3698` / `CVE-2026-59893` | high | 0.6.0 | Upgrade required |
| sqlparse 0.5.4 | `PYSEC-2026-3699` / `CVE-2026-54284` | high | 0.6.0 | Upgrade required |

Repository search finds no direct sqlparse call outside dependency/lock declarations. Django and tooling may still invoke parsing paths; absence of direct calls is not a risk-acceptance basis. Upgrade and regenerate the lock under `#69`.

## Temporary Django exception policy

A temporary exact allowlist for `PYSEC-2026-3717` is allowed only when:

```text
no GIS import/app/model/form/exposed parsing path is proven
owner and rationale are recorded
expiry <= 2026-09-30
issue #70 remains open and P0
shared/production promotion remains prohibited
```

No blanket ignore is permitted. The exception must be removed when `#70` merges.

## Previous 13 August triage

The earlier locked baseline removed 155 of 161 broad-environment findings and narrowly triaged six then-current Django advisories plus one development-only pytest advisory. Those decisions remain historical evidence but do not cover newly published advisory IDs.

Previously triaged exact IDs:

```text
GHSA-923m-gv2p-w5qp
GHSA-h7pc-vwp9-298g
GHSA-8cjm-8mp7-r2xf
GHSA-3h9f-r86x-qvjx
GHSA-crhf-3pfg-w68w
GHSA-8qcx-xf44-272x
GHSA-6w46-j5rx-g56g
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
#69 sqlparse upgrade and immediate exact advisory decision
→ #70 Django 5.2 LTS migration and removal of Django 4.2 exceptions
→ #13 protected aggregate CI
```

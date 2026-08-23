---
title: Dependency Audit Triage
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: maintainers and security reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
---

# Dependency audit triage

## Canonical status

```text
DJANGO = 5.2.17 LTS
SQLPARSE = 0.6.0
DJANGO_ADVISORY_EXCEPTIONS = zero
DEPENDENCY_AUDIT = pass
TEST_TOOL_EXCEPTION = GHSA-6w46-j5rx-g56g only
```

Merged PR `#71` upgraded sqlparse and established exact temporary Django 4.2 triage. Merged PR `#72` then migrated to Django 5.2.17 LTS and removed every Django 4.2 advisory exception.

The canonical dependency evidence includes:

```text
uv lock --check: PASS
frozen chemistry sync: PASS
pip-audit through scripts.audit_dependencies: PASS
ordinary Django check: PASS
synthetic manage.py check --deploy: PASS
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

## Resolved sqlparse findings

The direct and locked version is `sqlparse 0.6.0`, addressing:

```text
PYSEC-2026-3696 / CVE-2026-59894
PYSEC-2026-3697 / CVE-2026-71491
PYSEC-2026-3698 / CVE-2026-59893
PYSEC-2026-3699 / CVE-2026-54284
```

## Resolved Django finding and support state

Django 4.2.30 was the final 4.2 release and was outside upstream extended support. The former exact temporary exception for:

```text
PYSEC-2026-3717 / CVE-2026-15830
```

was removed when Django 5.2.17 LTS became canonical. The audit wrapper contains no Django advisory exception.

## Remaining exact development-tool decision

```text
GHSA-6w46-j5rx-g56g
```

This pytest temporary-directory advisory requires a malicious local user on a shared UNIX host. The compatible `pytest-playwright==0.7.1` release still constrains pytest below the patched major line. Tests run on isolated CI runners or single-user development machines.

This is a narrow development-only decision. It does not authorize private data in test artifacts, shared hostile-host execution or blanket advisory suppression. Remove it when the browser-test toolchain supports the patched pytest line.

## Required record fields for future findings

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
supported dependency baseline complete
→ #14 documentation source of truth
→ #13 protected aggregate CI hardening
→ #23 minimum scientific QA
```

---
title: Authentication and Access
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: maintainers, security reviewers and operators
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Authentication and access

## Current verified boundary

```text
public: documentation and selected simulator surfaces
protected: clinic, research, simulator-management and DRF/API surfaces
configured login target: Django admin login
normal product login: absent/incomplete
```

The repository uses Django authentication and selected `login_required` controls; DRF defaults to `IsAuthenticated`. This is not a complete product authentication architecture.

## Temporary smoke identity

The `ensure_m0_smoke_user` management command creates, rotates or disables a reserved temporary identity with:

```text
staff = false
superuser = false
groups = 0
direct_permissions = 0
```

It is an operational smoke identity, not a product role. It cannot use the Django admin login and must not be promoted to staff to bypass the missing product-login contract.

## Current authorization risk

Authorization logic exists in several locations, but patient-derived list/detail/dashboard and aggregate behavior is not yet proven to use one central selector/policy layer. Issue `#8` is release-blocking for shared use.

## Target role architecture — not current

Candidate least-privilege roles include:

```text
Research Viewer
Data Curator
Clinical Reviewer
Simulator Editor
Administrator
```

Permissions must distinguish view, edit, import, validate, simulate, export and administer. Role names do not become current until migrations, policies, selectors, tests and audit evidence exist.

## Required access invariant

Patient-derived objects should be resolved from a user-scoped selector before rendering or serialization:

```python
patient = get_object_or_404(accessible_patients(request.user), pk=pk)
```

The exact implementation is governed by issue `#8` and may differ, but post-load authorization after disclosure is not acceptable.

## Audit boundary

Future access events should record actor ID/role snapshot, action, pseudonymous resource reference, allowed/denied decision, policy reason, request ID and status without patient names, raw query strings or clinical payloads.

# Security policy

## Current support state

bmyCure4MM is an `E1_research_prototype`.

```yaml
production_security_certified: false
clinical_use_supported: false
product_auth_model_complete: false
object_authorization_complete: false
dependency_security_green: true
```

The canonical research dependency baseline uses Django 5.2.17 LTS and sqlparse 0.6.0. Dependency audit, synthetic deployment checks, migration evidence, Python 3.11/3.12 suites, numerical identity, Secret Scan, Docker build and the protected `required` result passed in merged PRs `#71` and `#72`.

This does not make the application production-certified or suitable for clinical deployment. Production security also depends on authentication, object authorization, settings, TLS/proxy behavior, storage, backup/restore, quotas, observability and incident response.

## Reporting a vulnerability

Do **not** disclose vulnerabilities or patient-derived information in a public issue.

Preferred reporting path:

1. use the repository's private GitHub Security Advisory reporting flow when available;
2. otherwise email the repository owner at **andrea.zedda@outlook.it** with the subject `bmyCure4MM security report`.

Do not attach real patient records, clinical PDFs, direct identifiers, credentials or production secrets. Use a minimal synthetic reproducer or offer the material through an agreed private channel.

Include:

- affected commit/version;
- vulnerability class;
- minimal reproduction steps;
- potential impact;
- whether patient-derived data, authentication, authorization or artifact access may be involved;
- suggested mitigation, if known.

Acknowledgment and remediation timing are targets, not a contractual SLA. Disclosure will be coordinated according to severity, exploitability and available mitigation.

## Verified current controls

- supported Django 5.2 LTS dependency baseline and deterministic `uv.lock`;
- repository-native dependency audit with exact documented development-tool triage only;
- environment-based secret-key configuration and weak/default-key rejection;
- Django password hashing and CSRF protections;
- selected authenticated views and global DRF `IsAuthenticated` default;
- secret-scanning workflow;
- repository hygiene and research safety tooling;
- ignored `local_private/` boundary for private research material;
- temporary M0 smoke identities constrained to non-staff, non-superuser accounts with no groups or direct permissions;
- synthetic CI deployment check and disposable migration evidence;
- protected `required` dependency/reproducibility result.

These controls do not prove complete application security.

## Known open security gaps

### Product authentication

Documentation and selected simulator surfaces are public, while clinic/research/simulator-management/API surfaces are protected. `LOGIN_URL` currently targets the Django admin login and there is no reviewed normal-user product-login and role contract.

A temporary smoke identity cannot authenticate through the admin login by design and must not be promoted to staff merely to make a smoke test pass.

### Object authorization and RBAC

GitHub issue `#8` tracks inconsistent patient-derived queryset scoping, dashboard aggregation, duplicated policy logic, least-privilege roles and privacy-safe access auditing.

### Production security

GitHub issue `#9` tracks fail-closed production settings, proxy/TLS behavior, HSTS, CSP, cookies, deployment checks and resource-specific framing policy.

### Rate and resource controls

GitHub issue `#10` tracks request throttling, job quotas, concurrency, timeout, deduplication and cost controls.

### CI and branch protection

A protected `required` status exists. GitHub issue `#13` still tracks immutable Action pinning, exact base/head PR safety mode, documentation contracts, PR templates, CODEOWNERS, privacy-safe artifact retention, final workflow rationalization and branch-protection re-read.

### Remaining development-tool advisory decision

The dependency audit contains one exact development-only pytest advisory decision, documented in `docs/operations/DEPENDENCY_AUDIT_TRIAGE.md`. It is not a runtime-framework exception and must be removed when the browser-test toolchain supports the patched pytest line.

## Data-safety rules

Never commit:

- names, medical-record numbers or dates of birth;
- clinical PDFs or source excerpts;
- private longitudinal dataset payloads;
- database files, media uploads or generated private artifacts;
- `.env`, keys, credentials or tokens.

Use synthetic/demo fixtures in tests and CI. Security evidence must identify findings without unnecessarily printing sensitive matched content.

## Deployment boundary

The repository may contain deployment automation, but a successful image build or running deployment does not establish production security. A deployment must separately prove configuration, supported dependencies, database, storage, backup/restore, observability, access control and incident-response requirements.

## Security-related issues

Use public issues only for non-sensitive hardening tasks after removing exploit details and private data. Use private reporting for vulnerabilities.

# Security policy

## Current support state

bmyCure4MM is an `E1_research_prototype`.

```yaml
production_security_certified: false
clinical_use_supported: false
product_auth_model_complete: false
object_authorization_complete: false
dependency_security_green: false
```

No current version is represented as production-certified or suitable for clinical deployment. Security fixes are handled on the active `master` research line, subject to maintainer capacity and explicit verification.

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

- environment-based secret-key configuration and weak/default-key rejection;
- Django password hashing and CSRF protections;
- selected authenticated views and global DRF `IsAuthenticated` default;
- secret-scanning workflow;
- repository hygiene and dependency-audit tooling;
- research safety checks for staged sensitive files and PHI-like markers;
- ignored `local_private/` boundary for private research material;
- temporary M0 smoke identities are constrained to non-staff, non-superuser accounts with no groups or direct permissions.

These controls do not prove complete application security.

## Known open security gaps

### Dependency and framework support

The deterministic lock currently pins unsupported Django 4.2.30 and sqlparse 0.5.4. Newly disclosed advisories make the dependency audit red.

- issue `#69` tracks immediate sqlparse upgrade and exact advisory triage;
- issue `#70` tracks durable migration to Django 5.2 LTS;
- shared/production promotion is prohibited while the framework baseline is unsupported or security findings remain untriaged;
- blanket audit suppression is prohibited.

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

GitHub issue `#13` tracks consolidation of existing workflows into one mandatory protected scientific gate, exact-PR safety behavior, immutable Action pinning and PR governance.

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

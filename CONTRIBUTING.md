# Contributing to bmyCure4MM

bmyCure4MM is an `E1_research_prototype`. Contributions must preserve the
distinction between education, research simulation, validated prediction,
causal inference, and clinical decision-making. Passing CI is necessary
software evidence; it is not proof of model validity or clinical utility.

## Privacy and safety boundary

Never commit or attach:

- PHI, direct identifiers, patient records, or source clinical documents;
- private research datasets or payloads from `local_private/`;
- databases, media, logs, generated reports, or raw request payloads;
- credentials, private keys, `.env` files, or production secrets.

Use synthetic fixtures and public sources only. Screenshots must contain
synthetic data. Report security or privacy vulnerabilities through the private
security-advisory link in the issue chooser, not a public issue.

## Supported environment

The supported Python versions are 3.11 and 3.12. `uv==0.12.3` is required,
`uv.lock` is authoritative, and stale or unlocked dependency graphs are not
accepted.

```bash
python3.11 -m pip install uv==0.12.3
uv sync --frozen --extra chemistry
uv run python manage.py migrate
```

Do not install ad hoc packages to make a check pass. Dependency changes must be
declared in `pyproject.toml`, resolved into `uv.lock`, and justified in the pull
request.

## Pull-request-only workflow

`master` is governed through pull requests. Do not push directly to `master`.

```bash
git switch master
git pull --ff-only origin master
git switch -c <type>/<short-description>
```

Keep the branch bounded to one linked issue. The pull-request template requires
the source of truth, scope, risk classification, scientific and operational
impact, migrations, numerical changes, privacy/authorization impact, test
evidence, rollback, and known limitations.

## Local checks

For a fast formatting, lint, and typing pass:

```bash
bash scripts/ci/run_required_checks.sh quality
```

For the complete provider-independent required gate, including Docker:

```bash
bash scripts/ci/run_required_checks.sh
```

Before committing, stage only the intended files and run the local safety mode:

```bash
git add <intended-files>
bash scripts/pre_push_research_safety_check.sh
git diff --check
git status --short
```

The safety script scans the staged diff plus repository-wide tracked-file
invariants. CI supplies explicit base/head SHAs and scans the exact candidate
diff. Invalid refs fail closed. Matched sensitive content is not printed, and
ignored `local_private/` content is never traversed.

Detailed phases, CI-only settings, the immutable action inventory, and
diagnostic commands are documented in
[Testing and M0-R CI Operations](docs/operations/TESTING.md).

## Required GitHub status

The mandatory workflow is `M0-R CI`. Its jobs are:

- `repository-hygiene`;
- `quality`;
- `django-integrity` on Python 3.11 and 3.12;
- `scientific-regression`;
- `dependency-security`;
- `container-build`;
- final `required` aggregation.

The stable branch-protection status is:

```text
M0-R CI / required
```

Do not merge a red, cancelled, skipped, or pending required status.

## Change-specific evidence

### Scientific model or numerical changes

Declare changes to equations, mechanisms, coefficients, parameters, solver
behavior, tolerances, units, observation semantics, epistemic language, and
comparability. Include:

- scientific question and public evidence provenance;
- model/version decision;
- deterministic before/after numerical diff;
- uncertainty, validation, and falsification plan;
- impact on historical runs and invalidation;
- allowed and forbidden conclusions.

Do not update the golden numerical baseline merely to make CI pass.

### Schema and migrations

Declare whether migrations are absent, included, or intentionally deferred.
Every pull request must pass:

```bash
uv run python manage.py makemigrations --check --dry-run --settings=mmportal.settings_ci
```

Schema changes require forward evidence, cardinality/constraint checks where
relevant, and a rollback path or an explicit irreversibility statement. Do not
silently reinterpret stored scientific data.

### Data, importers, and evidence

State dataset identity, schema version, provenance, units, missingness behavior,
idempotence, and invalidation impact. Tests and examples must use synthetic or
public fixtures. Claims and epistemic labels must remain consistent with the
canonical intended-use governance.

### Security and authorization

Describe object ownership, roles, denial behavior, secret handling, and artifact
access. Authorization failures must not disclose whether another user's object
exists. High-severity dependency findings cannot be suppressed without an
explicit reviewed exception and documented mitigation.

## Numerical-drift policy

The governed command is:

```bash
uv run python -m scripts.check_numerical_baseline
```

Unexpected drift is blocking. An intentional change requires reviewed output
evidence and a model-version decision before the fixture can change. CI success
does not convert a simulated, heuristic, literature-informed, or hypothetical
result into an externally validated result.

## Commits and review

Use concise English commit messages such as:

```text
ci: add mandatory M0-R quality gate
test: preserve public ChemTools privacy contract
docs: document numerical drift review
```

Resolve review conversations, keep the branch current with `master`, and rerun
the safety gate after the final staged change. Prefer a reversible rollback and
state residual risk honestly.

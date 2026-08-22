# Testing and M0-R CI Operations

## Purpose and boundary

The mandatory M0-R gate verifies the locked software environment, repository
hygiene, Django integrity, established scientific regression contracts,
dependency security, and the current container build. It uses synthetic test
data and temporary runtime paths only.

Passing this gate does not establish scientific or clinical validity:

```text
software verification
!=
model validation
!=
external validation
!=
clinical utility
```

The repository remains an `E1_research_prototype`. CI does not validate
patient-specific prediction, causal effects, or clinical decision support.

## Canonical environment

Python 3.11 and 3.12 are supported. `uv==0.12.3` and `uv.lock` are
authoritative. Install the complete locked graph used by the required checks:

```bash
uv sync --frozen --extra chemistry
```

Do not regenerate or bypass the lock in CI. Dependency changes must update both
`pyproject.toml` and `uv.lock` and must pass the dependency audit.

## Local commands

A fast local quality pass is:

```bash
bash scripts/ci/run_required_checks.sh quality
```

The full provider-independent required gate is:

```bash
bash scripts/ci/run_required_checks.sh
```

The full command includes the Docker build and therefore requires a running
Docker engine. It creates bounded temporary database, media, log, cache, and
static paths and verifies that the commands do not change Git working-tree
state.

Canonical individual commands remain available for diagnosis:

```bash
uv sync --frozen --extra chemistry
uv run python manage.py check --settings=mmportal.settings_ci
uv run python manage.py check --deploy --fail-level WARNING --settings=mmportal.settings_ci_deploy
uv run python manage.py makemigrations --check --dry-run --settings=mmportal.settings_ci
uv run python manage.py test --settings=mmportal.settings_ci
uv run ruff format --check .
uv run ruff check .
uv run mypy
uv run python -m scripts.check_numerical_baseline
uv run python -m scripts.audit_dependencies
bash scripts/pre_push_research_safety_check.sh
docker build .
```

`mmportal.settings_ci` is a synthetic, CI-only test overlay. It replaces
external secrets and uses local in-memory services and temporary paths.
`mmportal.settings_ci_deploy` adds strict security values solely for Django's
deployment-system checks. Neither module is a production settings architecture.

## Research safety diff modes

Local mode scans the staged candidate diff and repository-wide tracked-file
invariants:

```bash
git add <intended-files>
bash scripts/pre_push_research_safety_check.sh
```

CI mode fails closed unless both refs resolve to commits:

```bash
BMYCURE4MM_SAFETY_BASE_REF=<base-sha> \
BMYCURE4MM_SAFETY_HEAD_REF=<head-sha> \
bash scripts/pre_push_research_safety_check.sh
```

The gate does not traverse ignored `local_private/` content. It reports finding
locations without printing matched sensitive content. It rejects tracked
private/runtime artifacts, high-confidence secrets and private keys, prohibited
direct identifiers, absolute local paths in runtime code, and invalid Python
syntax. The known public `5LF3` structural report is the only tracked PDF
allowlisted by path; source clinical PDFs are forbidden.

## Workflow graph and required status

`.github/workflows/ci.yml` defines:

```text
repository-hygiene
quality
django-integrity (Python 3.11 and 3.12)
scientific-regression
dependency-security
container-build
        |
        v
required
```

Every job has a bounded timeout. The final job runs after all mandatory jobs,
including failure cases, and fails unless each result is `success`. The stable
branch-protection context is:

```text
M0-R CI / required
```

No production secret is provided to pull-request jobs. Workflow token
permissions default to `contents: read`. Superseded runs for the same pull
request or branch are cancelled.

## Immutable GitHub Action inventory

All external actions used by repository workflows are pinned to the following
verified full commit SHAs; release comments are repeated at each use site.

| Action | Release | Commit SHA |
| --- | --- | --- |
| `actions/checkout` | `v4.2.2` | `11bd71901bbe5b1630ceea73d27597364c9af683` |
| `actions/setup-python` | `v5.6.0` | `a26af69be951a213d495a4c3e4e4022e16d87065` |
| `astral-sh/setup-uv` | `v7.6.0` | `37802adc94f370d6bfd71619e3f0bf239e1f3b78` |
| `actions/upload-pages-artifact` | `v3.0.1` | `56afc609e74202658d3ffba0e8f6dda462b719fa` |
| `actions/deploy-pages` | `v4.0.5` | `d6db90164ac5ed86f2b6aed7e0febac5b3c0c03e` |
| `docker/setup-qemu-action` | `v3.6.0` | `29109295f81e9208d7d86ff1c6c12d2833863392` |
| `docker/setup-buildx-action` | `v3.11.1` | `e468171a9de216ec08956ac3ada2f0791b6bd435` |
| `docker/login-action` | `v3.4.0` | `74a5d142397b4f367a81961eba4e8cd7edddf772` |
| `docker/metadata-action` | `v5.8.0` | `c1e51972afc2121e065aed6d45c65596fe445f3f` |
| `docker/build-push-action` | `v6.18.0` | `263435318d21b8e681c14492fe198d362a7d2c83` |
| `gitleaks/gitleaks-action` | `v2.3.9` | `ff98106e4c7b2bc287b24eaf42907196329070c7` |
| `actions/upload-artifact` | `v4.6.2` | `ea165f8d65b6e75b540449e92b4886f43607fa02` |

## Scientific and migration policy

- The golden numerical fixture must not change in routine CI work.
- Any approved numerical change requires an output diff, model-version decision,
  invalidation/comparability assessment, and scientific review.
- Model, schema, data, evidence, and intended-use impacts must be declared in
  the pull request.
- Migration drift is always blocking. Required migrations need forward evidence
  and rollback or irreversibility documentation.
- A passing test process is software evidence only. Negative, inconclusive, or
  invalid scientific results must not be relabelled to make CI green.

Broader property, contract, migration, E2E, and coverage expansion remains in
issue #23.

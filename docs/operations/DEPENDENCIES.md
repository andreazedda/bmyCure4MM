# Dependency Operations

## Canonical environment

`pyproject.toml` declares direct dependencies and `uv.lock` is the only
authoritative transitive lock. The required resolver is `uv==0.12.3`. Python
3.11 and 3.12 are supported; Docker and the numerical reference gate use
Python 3.11.

The default environment contains runtime, test, quality, and security tools:

```bash
python3.11 -m pip install uv==0.12.3
uv sync --frozen
uv run pre-commit install
```

Molecular pipelines are deliberately optional:

```bash
uv sync --frozen --extra chemistry
```

Omitting `chemistry` must leave Django checks and non-chemistry tests usable.
Documentation and screenshot tooling are separate groups:

```bash
uv sync --frozen --no-default-groups --group docs
uv sync --frozen --no-default-groups --group screenshots --extra chemistry
```

The former root, module, and lab `requirements*.txt` files were removed because
they described mutually inconsistent NumPy, Matplotlib, RDKit, and urllib3
graphs. The old lab notebooks remain historical, not a supported scientific
environment.

## Reproducible commands

Use `--frozen` for installation, CI, Docker, tests, and research work. It fails
instead of changing the lock.

```bash
uv lock --check
uv sync --frozen --extra chemistry
uv run python -m scripts.reproducibility_report
```

The report records the Git SHA, lock SHA-256, Python/platform identity, and
scientifically material package versions. Release automation must provide an
OCI digest through `CONTAINER_IMAGE_DIGEST`; local reports fail closed to
`UNAVAILABLE_LOCAL` rather than inventing a digest. Issue #19 will bind these
fields into the immutable run manifest; until then, retain the report beside
every research artifact.

## System dependencies

| Surface | External requirement |
| --- | --- |
| Base runtime | Python 3.11/3.12 and system OpenSSL/CA certificates |
| Worker | Redis service; the Python client is locked |
| Chemistry | Prebuilt RDKit wheel on supported platforms; molecular viewers may require network access to public structure services |
| Container build | Debian `gcc`, `g++`, and `libpq-dev`; the Python graph still comes exclusively from `uv.lock` |
| PostgreSQL deployment | PostgreSQL service/client runtime when selected; SQLite remains the local default |

## Controlled updates

Never run an unconstrained upgrade on a scientific branch. For one dependency:

1. Save `uv run python -m scripts.reproducibility_report` and
   `uv run python -m scripts.check_numerical_baseline` output.
2. Change the exact direct version in `pyproject.toml` if it is direct.
3. Run `uv lock --upgrade-package <name>`.
4. Run every gate below and attach the numerical before/after diff to the PR.
5. If output exceeds the governed tolerance, explain the scientific cause,
   decide whether `MODEL_VERSION` changes, and list invalidated prior runs.

Do not hand-edit `uv.lock`.

## Required gates

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy
uv run python -m scripts.check_repository_hygiene
uv run python -m scripts.check_numerical_baseline
uv run python -m scripts.audit_dependencies
uv run python manage.py check
uv run python manage.py makemigrations --check --dry-run
```

Ruff currently governs the canonical scientific/reproducibility surface named
in `pyproject.toml`. The measured legacy baseline is 496 findings and 220 files
needing format; silently mass-rewriting that unrelated surface is outside issue
#15. `.quality-baseline.json` also content-addresses two pre-existing malformed
lab notebooks so any change or new malformed notebook fails closed.

## Rollback

Revert the dependency PR and run `uv sync --frozen` against the restored lock.
Do not partially restore `pyproject.toml` without its matching `uv.lock`.

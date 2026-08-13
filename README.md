# bmyCure4MM

bmyCure4MM is an **E1 research prototype** for reproducible Multiple Myeloma
computational research.

```yaml
intended_use_level: E1_research_prototype
clinical_decision_support: false
patient_specific_prediction_validated: false
causal_effect_identified: false
```

The platform supports mechanistic simulation, structured research data,
lineage-bound patient-twin experiments, exploratory scenario comparison, and
scientific QA. It does not choose treatment, issue patient-specific dose or
schedule instructions, predict patient benefit, or identify causal treatment
effects.

The canonical policy is [Canonical Intended Use](docs/governance/INTENDED_USE.md).

## Current capabilities

- versioned structured research-dataset import with content hashes;
- mechanistic PK/PD and tumor/healthy-cell simulation;
- lineage-bound research twin initialization and calibration infrastructure;
- model-relative diagnostic flags with explicit epistemic labels;
- time-respecting research diagnostics, uncertainty, sensitivity, robustness,
  and backtesting prototypes;
- molecular research utilities, with RDKit-dependent features optional;
- separate educational (`/learn/`), research (`/research/`), and administrative
  (`/admin/`) surfaces.

All scientific results require interpretation in model and data context.
Negative and inconclusive results are valid research outcomes.

## Local setup

Python 3.11 and 3.12 are supported. `uv.lock` is the sole dependency lock;
`uv` 0.12.3 is required and refuses a stale lock.

```bash
git clone https://github.com/andreazedda/bmyCure4MM.git
cd bmyCure4MM
python3.11 -m pip install uv==0.12.3  # one-time bootstrap
uv sync --frozen --extra chemistry
uv run python manage.py migrate
uv run python manage.py check
uv run python manage.py runserver
```

Omit `--extra chemistry` for core work that does not use RDKit or the molecular
pipelines. Dependency groups, lock updates, audits, and deterministic numerical
checks are documented in [Dependency Operations](docs/operations/DEPENDENCIES.md).

Do not use production secrets or private patient payloads in development
commands, fixtures, tests, or Git.

## Tests and safety gate

```bash
uv run python manage.py check
uv run python manage.py makemigrations --check --dry-run
uv run python manage.py test
uv run ruff check .
uv run ruff format --check .
uv run mypy
uv run python -m scripts.check_numerical_baseline
uv run python -m scripts.audit_dependencies
scripts/pre_push_research_safety_check.sh
```

The safety script checks ignored private paths, staged sensitive files,
PHI-like markers, Django configuration, twin-engine tests, and selected
authorization/simulation tests.

## Scientific identity

The current research model identifier is `research-twin-v1`. Twin states may
bind model/config/Git identity, dataset identity, and a content-addressed
computational-input hash. Historical states without lineage remain explicitly
unbound; they are not silently reinterpreted.

The mandatory cross-run manifest and model registry are tracked in GitHub issue
#19. Until then, do not infer comparability from a successful process exit.

## Epistemic labels

Current canonical labels are:

`OBSERVED`, `DERIVED`, `USER_PROVIDED`, `SIMULATED`, `HEURISTIC`,
`LITERATURE_INFORMED`, `HYPOTHETICAL`, and `VALIDATED_EXTERNAL`.

Definitions and restrictions are in
[Epistemic Labels](docs/governance/EPISTEMIC_LABELS.md) and
[Model Output Language](docs/governance/MODEL_OUTPUT_LANGUAGE.md).

## Repository areas

```text
clinic/         structured clinical-record workspace and data navigation
simulator/      mechanistic simulation and synthetic learning scenarios
twin_engine/    research twin, lineage, diagnostics, and artifacts
chemtools/      molecular research utilities
docs/           current and archived documentation
governance/     machine-readable release claims
local_private/  ignored private research material (never commit)
```

## Privacy

Names, medical-record numbers, dates of birth, clinical PDFs, source excerpts,
direct identifiers, private dataset payloads, and private artifacts must never
enter Git. Keep private research material under ignored `local_private/` paths
and export only explicitly validated, de-identified research artifacts.

## Documentation status

Current governance documents are authoritative. Historical feature narratives
live under `docs/archive/` and must not be interpreted as current behavior. A
full documentation source-of-truth rebuild is tracked in GitHub issue #14.

## License and security

See [LICENSE](LICENSE) and [SECURITY.md](SECURITY.md).

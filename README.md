# bmyCure4MM

bmyCure4MM is an **E1 research prototype** whose north-star is a virtual and reproducible Multiple Myeloma research laboratory.

```yaml
intended_use_level: E1_research_prototype
clinical_decision_support: false
patient_specific_prediction_validated: false
causal_effect_identified: false
```

The platform currently provides governed structured research data, lineage-bound Patient Twin infrastructure, mechanistic PK/PD and schedule-aware simulation, model-relative counterfactual research, immutable run identity, artifact integrity and scientific QA foundations.

It does **not** choose treatment, issue patient-specific dose or schedule instructions, predict validated patient benefit or identify causal treatment effects.

## Start here

- [Mission and north star](docs/product/MISSION_AND_NORTH_STAR.md)
- [Current verified state](docs/product/CURRENT_STATE.md)
- [Capabilities and limitations](docs/product/CAPABILITIES_AND_LIMITATIONS.md)
- [Canonical intended use](docs/governance/INTENDED_USE.md)
- [Source-of-truth policy](docs/governance/SOURCE_OF_TRUTH.md)
- [Current architecture](docs/architecture/CURRENT_SYSTEM.md)
- [Target virtual laboratory](docs/architecture/TARGET_VIRTUAL_LAB.md)
- [Current model registry and formulas](docs/models/REGISTRY.md)
- [Project roadmap](docs/product/ROADMAP.md)

## Current verified baseline

```text
repository = andreazedda/bmyCure4MM
branch = master
verified_head = bf097810b337dc6b766cda04497005670cd96513
```

Current M0-R path:

```text
#14 documentation/source of truth
+ #69 immediate dependency-advisory remediation
→ #70 supported Django 5.2 LTS baseline
→ #13 mandatory CI aggregate gate and PR governance
→ #23 minimum scientific/privacy invariants
→ #26 close M0-R
```

The current lock remains deterministic, but the security audit now detects newly disclosed advisories in Django 4.2.30 and sqlparse 0.5.4. The durable Django migration is tracked in `#70`; no red dependency gate may be relabelled green through blanket suppression.

## Current capabilities

- versioned Structured Research Dataset Contract and idempotent import infrastructure;
- content-addressed computational-input and Twin lineage;
- persistent Twin states, residuals and calibration infrastructure;
- logistic tumour/healthy-cell and PK/PD research simulation;
- day-resolved administered-dose profiles and schedule identity;
- heuristic hepatic and neutropenia risk signals;
- mechanistic baseline/alternative counterfactual runs;
- immutable run manifests, model registry, artifact hashes, comparability and invalidation;
- uncertainty, sensitivity, robustness and backtesting prototypes;
- optional molecular research utilities;
- separate learning, research and administration surfaces.

See [Capabilities and limitations](docs/product/CAPABILITIES_AND_LIMITATIONS.md) before interpreting any result.

## Local setup

Python 3.11 and 3.12 are supported. `uv.lock` is the sole dependency lock and `uv` 0.12.3 is required.

```bash
git clone https://github.com/andreazedda/bmyCure4MM.git
cd bmyCure4MM
python3.11 -m pip install uv==0.12.3
uv sync --frozen --extra chemistry
uv run python manage.py migrate
uv run python manage.py check
uv run python manage.py runserver
```

Omit `--extra chemistry` for core work that does not require RDKit. Dependency operations are documented in [Dependency Operations](docs/operations/DEPENDENCIES.md). The current Django 4.2 lock is unsupported and must not be promoted to a shared/production baseline before issues `#69` and `#70` are resolved.

## Verification and safety

```bash
uv lock --check
uv sync --frozen --extra chemistry
uv run python manage.py check
uv run python manage.py makemigrations --check --dry-run
uv run python manage.py test
uv run ruff check .
uv run ruff format --check .
uv run mypy
uv run python -m scripts.check_repository_hygiene
uv run python -m scripts.check_numerical_baseline
uv run python -m scripts.audit_dependencies
bash scripts/pre_push_research_safety_check.sh
```

Current test boundaries, the active dependency blocker and CI limitations are documented in [Testing](docs/operations/TESTING.md) and [Dependency Audit Triage](docs/operations/DEPENDENCY_AUDIT_TRIAGE.md).

## Privacy

Names, medical-record numbers, dates of birth, clinical PDFs, source excerpts, direct identifiers, private dataset payloads and private artifacts must never enter Git. Keep private research material under ignored `local_private/` paths and export only explicitly validated, de-identified research artifacts.

## Documentation status

The documentation corpus is governed by:

- [Documentation policy](docs/governance/DOCUMENTATION_POLICY.md)
- [Documentation inventory](docs/DOCUMENTATION_INVENTORY.yaml)

Historical feature and release narratives are retained only as archives and are not current product claims.

## License and security

See [LICENSE](LICENSE) and [SECURITY.md](SECURITY.md).

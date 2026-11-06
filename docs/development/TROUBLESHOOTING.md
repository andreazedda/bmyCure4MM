---
title: Manual Diagnosis and Recovery
language: en
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: developers and operators
last_verified_at: 2026-11-06
---

# Manual diagnosis and recovery

These exercises run from the existing checkout with the locked environment. They require no AI account, API token, production credential or patient data. Keep the [code atlas](CODE_ATLAS.md) and [deployment runbook](https://github.com/andreazedda/bmyCure4MM/blob/master/deploy/k8s/DEPLOY_K3S.md) open beside this page.

## Establish which system you are diagnosing

Record `git rev-parse HEAD`, Python version and the failing command. Distinguish your source checkout from the running image SHA/digest, runtime settings, database and storage. `mmportal.settings_test` is an isolated test configuration: its local mail backend, fast password hasher and eager task behavior cannot prove production delivery or broker health.

For a local audit terminal:

```bash
uv sync --frozen --extra chemistry
export DJANGO_SETTINGS_MODULE=mmportal.settings_test
export DJANGO_DEBUG=1
export DJANGO_SECURE_SSL_REDIRECT=0
unset BMY_TEST_DATABASE_NAME
export BMY_TEST_ROOT="$(mktemp -d /tmp/bmy-manual-audit.XXXXXX)"
uv run --frozen python manage.py migrate --noinput
uv run --frozen python manage.py ensure_learning_example
uv run --frozen python manage.py check
```

`BMY_TEST_DATABASE_NAME` is deliberately cleared: otherwise an inherited override can defeat the disposable database path. Do not export production settings or credentials into this terminal. The database and media used by `settings_test` stay under the temporary root; the in-memory email backend never sends external email. Starting `runserver` does not itself start Redis or Celery.

## Walkthrough A: account lifecycle and private ticket

```bash
uv run --frozen python manage.py test \
  accounts.tests.test_lifecycle \
  accounts.tests.test_concurrency \
  accounts.tests.test_existing_site.ExistingSiteTests.test_ticket_reply_is_private_and_survives_relogin \
  --noinput
```

Read the test source while executing. The fixture builds only synthetic identities, consumes distinct verify/activate capabilities, sets a password and makes ordinary HTTP form requests. Expected outcomes include a rejected public privilege field, failed reused/wrong-purpose token, no concurrent duplicate profile, and a ticket reply visible to its author but not the second user. This is an executable contract walkthrough; it is not a real owner approval or a recipient-mailbox check.

To inspect the decision path with Python's debugger:

```bash
rg -n '^def confirm_capability' accounts/services.py
uv run --frozen python -m pdb manage.py test accounts.tests.test_lifecycle --noinput
```

At the debugger prompt, set `b accounts/services.py:<line printed above>`, then `c`. Inspect only safe state such as `purpose` and branch conditions with `n`; never use `pp locals()` on authentication/email code, print `raw`, inspect passwords, decrypt payloads or copy full request headers. Use `q` to leave the debugger. The source atlas links the exact baseline implementation if line numbers have moved.

For the real desktop/mobile browser flow, run the explicit [browser command in the handbook](DEVELOPMENT.md#verification-commands-and-what-they-prove). It supplies the SMTP sink and disposable users automatically. An empty development database plus a locmem email backend does not provide a usable manual owner invitation by itself.

## Walkthrough B: existing calculation and private artifact

```bash
uv run --frozen python manage.py test \
  accounts.tests.test_existing_site.ExistingSiteTests.test_registered_user_runs_existing_engine_and_cannot_export_another_users_run \
  --noinput
```

Follow `simulator/urls.py` → `simulator/views_manage.py:simulate_scenario` → `SimulationAttempt.run_model` → `_resolve_solver_inputs` → `MathematicalModel.simulate`. The test creates an independent synthetic scenario and computes real results. Then it checks an authorized plot/CSV, another user's denial and suspended-user denial. The shared Scenario is the input template; the user's SimulationAttempt owns the saved result. A stored `/media/...` locator must be resolved through the authorized view, not published directly.

## Reproduce a numerical failure and recover

This direct kernel exercise is independent of any database. The synthetic zero-effect input is solely a mathematical check, not a treatment scenario or efficacy comparison.

```bash
uv run --frozen python - <<'PY'
from dataclasses import replace
from math import exp, isclose
import numpy as np
from simulator.models_simulation import MathematicalModel

model = MathematicalModel(
    baseline_tumor_cells=1000.0,
    baseline_healthy_cells=2000.0,
    drug_doses={"synthetic_zero": 0.0},
    pk_params={"synthetic_zero": {"half_life": 24.0, "Vd": 1.0}},
    pd_params={"synthetic_zero": {"Emax": 0.0, "EC50": 1.0}},
    growth_rates={"tumor": 0.023, "healthy": 0.015},
    interaction_matrix=np.zeros((1, 1)),
    time_span=(0.0, 1.0),
    evaluation_points=11,
)
trajectory = model.simulate()
expected = 10000.0 / (1.0 + 9.0 * exp(-0.023))
observed = float(trajectory["tumor_cells"].iloc[-1])
assert isclose(observed, expected, rel_tol=1e-6)
try:
    replace(model, time_span=(0.0, 0.0)).simulate()
except ValueError as error:
    assert str(error) == "time_span end must be greater than start."
    print("EXPECTED_INVALID_HORIZON=PASS")
else:
    raise AssertionError("Invalid horizon was accepted")
assert len(model.simulate()) == 11
print(f"LOGISTIC_CHECK=PASS endpoint={observed:.6f}; RECOVERY=PASS")
PY
```

The recovery is to correct the invalid input; it does not involve altering equations, increasing quotas or overwriting a golden baseline. This sample checks one zero-effect solution and one failure mode. The [numerical baseline](../research/NUMERICAL_REPRODUCIBILITY.md) covers additional governed cases.

## Reproduce a backup and restore without a clinical database

This exercise uses a new in-memory synthetic table and its own temporary backup. It never opens the application database, production volumes or a restored patient dataset.

```bash
uv run --frozen python - <<'PY'
import sqlite3
from contextlib import closing
from pathlib import Path
from tempfile import TemporaryDirectory

with TemporaryDirectory(prefix="bmy-restore-exercise-") as directory:
    backup_path = Path(directory) / "synthetic.sqlite3"
    with closing(sqlite3.connect(":memory:")) as source:
        source.execute("CREATE TABLE exercise (id INTEGER PRIMARY KEY, state TEXT)")
        source.execute("INSERT INTO exercise VALUES (1, 'before-change')")
        source.commit()
        with closing(sqlite3.connect(backup_path)) as backup:
            source.backup(backup)
        source.execute("UPDATE exercise SET state='after-change'")
        source.commit()
        with closing(sqlite3.connect(backup_path)) as backup, closing(sqlite3.connect(":memory:")) as restored:
            backup.backup(restored)
            assert restored.execute("PRAGMA quick_check").fetchone()[0] == "ok"
            assert restored.execute("SELECT state FROM exercise").fetchone()[0] == "before-change"
            assert source.execute("SELECT state FROM exercise").fetchone()[0] == "after-change"
    print("SYNTHETIC_BACKUP_RESTORE=PASS; source mutation preserved; temporary files removed")
PY
```

This proves the documented SQLite API sequence on a disposable fixture; it is not a production disaster-recovery drill. Production application rollback repins the previously verified image while preserving additive schema and existing secrets. Restoring actual account data also requires session/token quarantine **before** reopening traffic; follow the [account recovery contract](../operations/AUTHENTICATION_AND_ACCESS.md#bootstrap-recovery-and-restore).

## Diagnose by symptom

| Symptom | Check first | Correct next action |
| --- | --- | --- |
| Learning catalogue has no examples | Existing Scenario count and activity; source/runtime identity | Controlled `ensure_learning_example` populates only an empty catalogue with one independent synthetic example; inspect existing inactive entries rather than overriding them |
| Page is 200 but unstyled / formulas are raw | Asset requests, MIME type, failed script loading, static files in image and collectstatic output | Fix the specific asset/configuration; an HTML response alone is not success |
| Request shows waiting; no email yet | Product console notification state, attempts and sanitized error code | Retry eligible failed notifications in console; `sent` is provider acceptance only |
| Owner console asks for MFA again | Fifteen-minute privileged-session freshness; device time | Enter a fresh code; do not reset the owner's authenticator automatically |
| Registered user gets 403 at clinical/research tools | Intended role boundary | Use learning workflow; access approval does not grant clinical/staff privileges |
| CSV or plot link fails | Current access/expiry, author, artifact existence and safe resolved path | Use export/embed route; never expose the media volume to bypass denial |
| `Invalid solver input ...` | The exact field and validation branch in `_resolve_solver_inputs` | Correct the structured input; keep scientific defaults and equations unchanged |
| Python `AttributeError` for a model method | Source/image mismatch, import path, actual class definition | Align source/dependencies and restart only the owned local process; migrations do not add Python methods |
| Migration drift warning | `makemigrations --check --dry-run` under the exact settings | Diagnose the actual diff; never fabricate a production migration from a warning alone |
| Molecular job remains queued | Existing broker/worker reachability and job error code; dependency profile | Diagnose the named job and configured broker; do not infer health from account CronJob success |
| RDKit import fails / worker crashes | Locked Python/NumPy/RDKit versions and reproduction on an isolated local process | Use the frozen chemistry profile; `--pool=solo` can isolate local fork issues but is not a production capacity fix |
| Cluster API says not ready or NFS operations wait | Control-plane readiness, eligible node, mount latency and named release Job | Preserve existing workload/storage and report the exact blocker; no global restart, etcd reset or PVC deletion |

If a local molecular worker is needed, start a Redis instance and worker you own using explicit local settings. `redis-cli -h 127.0.0.1 -p 6379 ping` only tests that endpoint. `uv run --frozen celery -A mmportal worker --pool=solo --loglevel=info` runs one worker; it does not validate queue fairness, cancellation or production scheduling. Under `settings_test`, eager execution is already set in Python. The production settings file does not parse a `CELERY_TASK_ALWAYS_EAGER` environment value; setting that shell variable alone is not a supported switch.

## Evidence to retain and share

Record the source SHA, exact command, isolated backend, start/end time, exit result and the smallest useful sanitized error code. Keep private logs under their normal access controls. Share no email capability, MFA seed, password, cookie, authorization header, clinical row or full request dump. A complete log upload is not a safe default troubleshooting step.

The audit executed both walkthrough commands, the expected numerical failure/recovery and synthetic backup/restore. Exact current results are in the [audit verification record](AUDIT_VERIFICATION_2026_11_06.md). A human maintainer's independent no-AI exercise remains a separate acceptance step, not something the automated run can certify.

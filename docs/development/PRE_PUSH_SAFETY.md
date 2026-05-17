# Pre-push Research Safety

Run the safety gate before pushing research cockpit or twin-engine changes:

```bash
scripts/pre_push_research_safety_check.sh
```

The script is offline and performs these checks:

- Confirms `local_private`, `media`, `db.sqlite3`, and `.env` are ignored.
- Blocks tracked sensitive paths such as local DBs, private folders, logs, keys, and credential files.
- Blocks staged sensitive paths such as DBs, private/media/log folders, PDFs, images, and key material.
- Scans staged text for obvious direct identifiers and PHI-like markers.
- Runs `./venv/bin/python manage.py check`.
- Runs `./venv/bin/python manage.py test twin_engine.tests`.
- Includes cockpit content-clarity tests for section explanations, formulas, glossary, provenance, developer-console guidance, and forbidden claim language.
- Runs `./venv/bin/python manage.py test clinic.tests.test_patient_crud clinic.tests.test_security_and_ux simulator.tests.test_twin_pr1_integration simulator.tests.test_simulation`.

Use `PYTHON_BIN=/path/to/python scripts/pre_push_research_safety_check.sh` if the local virtual environment lives somewhere else.

The scanner is intentionally conservative. If it blocks a staged file, unstage or move the sensitive artifact into an ignored local path before pushing.
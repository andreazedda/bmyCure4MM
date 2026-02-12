"""Generate UI screenshots for MkDocs documentation.

This script is meant to be run locally or in CI. It will:

- Create an isolated SQLite DB via DJANGO_SQLITE_PATH
- Run migrations + load demo fixtures
- Create a demo admin user (admin/admin123)
- Start Django dev server
- Use Playwright (Chromium) to capture screenshots

Output images are written under docs/assets/images/screenshots/.
"""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Final

import requests
from playwright.sync_api import sync_playwright


ROOT: Final[Path] = Path(__file__).resolve().parents[3]
OUT_DIR_DEFAULT: Final[Path] = ROOT / "docs" / "assets" / "images" / "screenshots"


def _run(cmd: list[str], env: dict[str, str]) -> None:
    subprocess.run(cmd, cwd=str(ROOT), env=env, check=True)


def _wait_http_ok(url: str, timeout_s: int = 45) -> None:
    deadline = time.time() + timeout_s
    last_err: Exception | None = None
    while time.time() < deadline:
        try:
            r = requests.get(url, timeout=2)
            if 200 <= r.status_code < 500:
                return
        except Exception as e:  # noqa: BLE001
            last_err = e
        time.sleep(0.5)
    raise RuntimeError(f"Server not reachable at {url}. Last error: {last_err}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="http://127.0.0.1:8001")
    parser.add_argument("--out-dir", default=str(OUT_DIR_DEFAULT))
    parser.add_argument("--port", type=int, default=8001)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="bmycure4mm-screens-") as tmp:
        tmp_path = Path(tmp)
        db_path = tmp_path / "db.sqlite3"
        media_root = tmp_path / "media"
        media_root.mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        env.update(
            {
                "DJANGO_DEBUG": "1",
                "DJANGO_SQLITE_PATH": str(db_path),
                "DJANGO_MEDIA_ROOT": str(media_root),
                "ALLOWED_HOSTS": "127.0.0.1,localhost",
                "PYTHONUNBUFFERED": "1",
            }
        )

        _run([sys.executable, "manage.py", "migrate", "--noinput"], env=env)
        _run(
            [
                sys.executable,
                "manage.py",
                "loaddata",
                "clinic/fixtures/demo_patients.json",
            ],
            env=env,
        )

        _run(
            [
                sys.executable,
                "manage.py",
                "shell",
                "-c",
                "from django.contrib.auth import get_user_model; "
                "User=get_user_model(); "
                "User.objects.filter(username='admin').exists() or "
                "User.objects.create_superuser('admin','admin@example.com','admin123')",
            ],
            env=env,
        )

        # Create a minimal scenario + run one attempt so we can screenshot a real output plot.
        _run(
            [
                sys.executable,
                "manage.py",
                "shell",
                "-c",
                "from django.contrib.auth import get_user_model; "
                "from simulator import models; "
                "User=get_user_model(); u=User.objects.get(username='admin'); "
                "s=models.Scenario.objects.create(title='Screenshot scenario', clinical_stage='newly_diagnosed', summary='Auto-generated', risk_stratification='Standard', guideline_notes=''); "
                "a=models.SimulationAttempt.objects.create(scenario=s, user=u, parameters={'baseline_tumor_cells':1e9,'baseline_healthy_cells':5e11,'lenalidomide_dose':25.0,'bortezomib_dose':1.3,'daratumumab_dose':16.0,'time_horizon':60.0,'tumor_growth_rate':0.023,'healthy_growth_rate':0.015,'interaction_strength':0.05,'immune_compromise_index':1.0}); "
                "a.run_model(); print(a.pk)",
            ],
            env=env,
        )

        # Create a ChemTools binding job so we can screenshot a real molecule viewer inside Job detail.
        _run(
            [
                sys.executable,
                "manage.py",
                "shell",
                "-c",
                "from django.contrib.auth import get_user_model; "
                "from chemtools import models; "
                "from chemtools.tasks import run_binding_viz_job; "
                "User=get_user_model(); u=User.objects.get(username='admin'); "
                "job=models.ChemJob.objects.create(kind=models.ChemJob.BIND, input_a='5LF3', input_b='', user=u); "
                "run_binding_viz_job.run(job.pk, '5LF3', ''); "
                "print(job.pk)",
            ],
            env=env,
        )

        # Start dev server
        server = subprocess.Popen(
            [
                sys.executable,
                "manage.py",
                "runserver",
                f"127.0.0.1:{args.port}",
                "--noreload",
            ],
            cwd=str(ROOT),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        try:
            _wait_http_ok(f"{args.base_url}/healthz/")

            with sync_playwright() as p:
                browser = p.chromium.launch()
                context = browser.new_context(viewport={"width": 1440, "height": 900})
                page = context.new_page()

                def snap(path: str, url_path: str) -> None:
                    page.goto(f"{args.base_url}{url_path}", wait_until="networkidle")
                    page.wait_for_timeout(500)
                    page.screenshot(path=str(out_dir / path), full_page=True)

                # Public-facing pages
                snap("public_docs_home.png", "/docs/")
                snap("public_simulator_list.png", "/sim/")
                snap("public_chemtools_home.png", "/chem/")

                # Authenticated clinician view (via Django admin login)
                page.goto(f"{args.base_url}/admin/login/", wait_until="domcontentloaded")
                page.fill("#id_username", "admin")
                page.fill("#id_password", "admin123")
                page.click("input[type=submit]")
                page.wait_for_load_state("networkidle")

                snap("clinic_patient_list.png", "/patients/")
                snap("clinic_dashboard.png", "/")

                # ChemTools (authenticated) - show the tools landing and a real binding viewer job detail.
                snap("chemtools_home.png", "/chem/")
                snap("chemtools_binding_job.png", "/chem/job/1/")

                # Real simulation plot output (generated above)
                snap("simulation_plot.png", "/media/sim_plots/attempt_1.html")

                context.close()
                browser.close()

        finally:
            if server.poll() is None:
                # Try graceful stop, then hard kill
                server.send_signal(signal.SIGINT)
                try:
                    server.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    server.kill()
                    server.wait(timeout=5)

            if server.stdout is not None:
                # If something failed, print server output to help debugging.
                output = server.stdout.read().strip()
                if output:
                    print("\n--- Django server output ---\n" + output, file=sys.stderr)

    # Ensure files exist (helps CI fail loudly)
    expected = [
        "public_docs_home.png",
        "public_simulator_list.png",
        "public_chemtools_home.png",
        "clinic_patient_list.png",
        "clinic_dashboard.png",
        "simulation_plot.png",
        "chemtools_home.png",
        "chemtools_binding_job.png",
    ]
    missing = [p for p in expected if not (out_dir / p).exists()]
    if missing:
        raise RuntimeError(f"Missing screenshots: {missing}")

    # Remove empty dirs in case the script was pointed elsewhere
    for parent in [out_dir, out_dir.parent]:
        try:
            if parent.exists() and not any(parent.iterdir()):
                shutil.rmtree(parent)
        except OSError:
            pass

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

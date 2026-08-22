from __future__ import annotations

import os
import shutil
import stat
import subprocess
import tempfile
from pathlib import Path
from unittest import TestCase


ROOT = Path(__file__).resolve().parents[1]
SAFETY_SCRIPT = ROOT / "scripts" / "pre_push_research_safety_check.sh"
REQUIRED_WRAPPER = ROOT / "scripts" / "ci" / "run_required_checks.sh"


class SafetyDiffModeTests(TestCase):
    def setUp(self) -> None:
        self._temp = tempfile.TemporaryDirectory(prefix="bmycure4mm-safety-test-")
        self.repo = Path(self._temp.name)
        self._git("init", "-q")
        self._git("config", "user.name", "CI Safety Test")
        self._git("config", "user.email", "ci-safety@example.invalid")
        (self.repo / ".gitignore").write_text(
            ".env\n.env.*\n!.env.example\nlocal_private/\nmedia/\nlogs/\n*.sqlite3\n",
            encoding="utf-8",
        )
        (self.repo / "README.md").write_text("Synthetic fixture repository.\n", encoding="utf-8")
        self._git("add", ".gitignore", "README.md")
        self._git("commit", "-qm", "initial synthetic fixture")
        self.base = self._git("rev-parse", "HEAD").stdout.strip()

    def tearDown(self) -> None:
        self._temp.cleanup()

    def _git(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", *args],
            cwd=self.repo,
            check=True,
            capture_output=True,
            text=True,
        )

    def _run_safety(self, *, base: str | None = None, head: str | None = None):
        env = os.environ.copy()
        env.update(
            {
                "BMYCURE4MM_SAFETY_SCAN_ONLY": "1",
                "PYTHON_BIN": shutil.which("true") or "/usr/bin/true",
            }
        )
        if base is not None:
            env["BMYCURE4MM_SAFETY_BASE_REF"] = base
        if head is not None:
            env["BMYCURE4MM_SAFETY_HEAD_REF"] = head
        return subprocess.run(
            ["bash", str(SAFETY_SCRIPT)],
            cwd=self.repo,
            env=env,
            check=False,
            capture_output=True,
            text=True,
        )

    def _commit(self, path: str, content: str) -> str:
        (self.repo / path).write_text(content, encoding="utf-8")
        self._git("add", path)
        self._git("commit", "-qm", f"add {path}")
        return self._git("rev-parse", "HEAD").stdout.strip()

    def test_ci_mode_rejects_prohibited_tracked_change_without_echoing_content(self) -> None:
        marker = "AK" + "IA" + "0123456789ABCDEF"
        head = self._commit("synthetic_fixture.txt", f"credential={marker}\n")

        result = self._run_safety(base=self.base, head=head)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("synthetic_fixture.txt", result.stderr)
        self.assertNotIn(marker, result.stdout + result.stderr)

    def test_ci_mode_accepts_clean_synthetic_change(self) -> None:
        head = self._commit("synthetic_fixture.txt", "synthetic=true\n")

        result = self._run_safety(base=self.base, head=head)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("ci mode", result.stdout)

    def test_ci_mode_fails_closed_for_invalid_ref(self) -> None:
        result = self._run_safety(base="missing-base-ref", head=self.base)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid safety base ref", result.stderr)

    def test_local_mode_still_scans_the_staged_diff(self) -> None:
        marker = "Data" + " di nascita"
        (self.repo / "synthetic_fixture.txt").write_text(f"{marker}: synthetic\n", encoding="utf-8")
        self._git("add", "synthetic_fixture.txt")

        result = self._run_safety()

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("synthetic_fixture.txt", result.stderr)
        self.assertNotIn(marker, result.stdout + result.stderr)


class RequiredCheckWrapperTests(TestCase):
    def test_required_wrapper_propagates_uv_failure(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bmycure4mm-wrapper-test-") as temp_dir:
            fake_uv = Path(temp_dir) / "uv"
            fake_uv.write_text("#!/usr/bin/env bash\nexit 23\n", encoding="utf-8")
            fake_uv.chmod(fake_uv.stat().st_mode | stat.S_IXUSR)
            env = os.environ.copy()
            env["PATH"] = f"{temp_dir}{os.pathsep}{env['PATH']}"
            result = subprocess.run(
                ["bash", str(REQUIRED_WRAPPER), "quality"],
                cwd=ROOT,
                env=env,
                check=False,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.returncode, 23, result.stdout + result.stderr)

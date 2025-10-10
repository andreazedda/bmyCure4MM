from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

from django.contrib.auth import get_user_model
from django.test import SimpleTestCase, TestCase, override_settings
from django.urls import reverse
from unittest.mock import patch

from chemtools import models, utils


class UtilsTests(SimpleTestCase):
    def test_run_settings_generator_timeout(self) -> None:
        with patch("chemtools.utils.subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="cmd", timeout=1)):
            with self.assertRaises(utils.ToolRunError) as exc:
                utils.run_settings_generator()
            self.assertIn("timed out", str(exc.exception).lower())

    def test_run_paths_generator_non_zero_exit_logs(self) -> None:
        fake_result = SimpleNamespace(returncode=1, stdout="stdout log", stderr="stderr log")
        with patch("chemtools.utils.subprocess.run", return_value=fake_result):
            with self.assertRaises(utils.ToolRunError) as exc:
                utils.run_paths_generator()
            self.assertIn("stdout log", exc.exception.stdout)
            self.assertIn("stderr log", exc.exception.stderr)

    def test_run_drug_parameter_evaluator_collects_latest(self) -> None:
        media_root_path = Path(tempfile.mkdtemp(prefix="chem_media_root_"))
        pipeline_root = Path(self._create_temp_structure())

        html_output = pipeline_root / "data" / "outputs" / "result.html"

        def fake_run(command, **kwargs):
            html_output.write_text("<html></html>")
            return SimpleNamespace(returncode=0, stdout="ok", stderr="")

        with override_settings(MEDIA_ROOT=str(media_root_path)):
            with (
                patch.object(utils, "PIPELINES_DIR", pipeline_root),
                patch("chemtools.utils.subprocess.run", side_effect=fake_run),
            ):
                target_dir = media_root_path / "tmpdrug"
                target_dir.mkdir(parents=True, exist_ok=True)
                result_path, stdout, stderr = utils.run_drug_parameter_evaluator(smiles="CCO", cid=123, outdir=target_dir)

            self.assertTrue(result_path and result_path.exists())
            self.assertIn("drug_123", result_path.name)
            self.assertEqual(stdout, "ok")
            self.assertEqual(stderr, "")
            cache_copy = utils._drug_cache_dir() / result_path.name
            self.assertTrue(cache_copy.exists())

        shutil.rmtree(media_root_path, ignore_errors=True)
        shutil.rmtree(pipeline_root, ignore_errors=True)

    def test_resolve_outdir_prefers_env_override(self) -> None:
        media_root_path = Path(tempfile.mkdtemp(prefix="chem_media_root_"))
        try:
            env_dir = media_root_path / "env_override"
            with override_settings(MEDIA_ROOT=str(media_root_path)):
                with patch.dict(os.environ, {"MM_OUTDIR": str(env_dir)}, clear=False):
                    resolved = utils._resolve_outdir(None, "drug_tmp")
                self.assertEqual(resolved, env_dir)
                self.assertTrue(env_dir.exists())

                explicit_dir = media_root_path / "explicit"
                result = utils._resolve_outdir(explicit_dir, "drug_tmp")
                self.assertEqual(result, explicit_dir)
        finally:
            shutil.rmtree(media_root_path, ignore_errors=True)

    def _create_temp_structure(self) -> str:
        root = Path(tempfile.mkdtemp(prefix="chem_pipe_"))
        (root / "processes").mkdir(parents=True, exist_ok=True)
        (root / "configs").mkdir(parents=True, exist_ok=True)
        (root / "data" / "outputs").mkdir(parents=True, exist_ok=True)
        (root / "processes" / "drug_parameter_evaluator.py").write_text("print('stub')\n")
        (root / "configs" / "drug_parameter_evaluator.yaml").write_text("drug_id: '123'\n")
        return str(root)


class JobStatusViewTests(TestCase):
    def setUp(self) -> None:
        self.media_root = Path(tempfile.mkdtemp(prefix="chem_media_status_"))
        self.addCleanup(shutil.rmtree, self.media_root, ignore_errors=True)
        self.user = get_user_model().objects.create_user("viewer", password="password")

    def test_job_status_payload(self) -> None:
        with override_settings(MEDIA_ROOT=str(self.media_root), MEDIA_URL="/media/"):
            job = models.ChemJob.objects.create(
                kind=models.ChemJob.PARAM,
                input_a="CCO",
                user=self.user,
            )
            job_dir = self.media_root / "chem" / str(job.pk)
            job_dir.mkdir(parents=True, exist_ok=True)
            html_path = job_dir / f"drug_{job.pk}_v1.html"
            html_path.write_text("<html></html>", encoding="utf-8")
            job.out_html.name = html_path.relative_to(self.media_root).as_posix()
            job.save(update_fields=["out_html"])

            self.client.force_login(self.user)
            response = self.client.get(reverse("chemtools:job_status", args=[job.pk]))
            self.assertEqual(response.status_code, 200)

            payload = response.json()
            self.assertEqual(payload["status"], "Completed")
            self.assertEqual(payload["variant"], "success")
            self.assertTrue(payload["has_html"])
            self.assertFalse(payload["has_csv"])
            self.assertTrue(payload["html_url"].endswith(f"drug_{job.pk}_v1.html"))
            self.assertIsNone(payload["thumbnail_url"])

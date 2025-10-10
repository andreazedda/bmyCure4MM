from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace

from django.test import SimpleTestCase, override_settings
from unittest.mock import patch

from chemtools import utils


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

    def _create_temp_structure(self) -> str:
        root = Path(tempfile.mkdtemp(prefix="chem_pipe_"))
        (root / "processes").mkdir(parents=True, exist_ok=True)
        (root / "configs").mkdir(parents=True, exist_ok=True)
        (root / "data" / "outputs").mkdir(parents=True, exist_ok=True)
        (root / "processes" / "drug_parameter_evaluator.py").write_text("print('stub')\n")
        (root / "configs" / "drug_parameter_evaluator.yaml").write_text("drug_id: '123'\n")
        return str(root)

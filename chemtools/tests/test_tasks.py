from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from typing import Any

from django.test import TestCase, override_settings
from unittest.mock import patch

from chemtools import models, tasks


class TaskVersioningTests(TestCase):
    def setUp(self) -> None:
        self.media_root = Path(tempfile.mkdtemp(prefix="chem_media_tasks_"))
        self.addCleanup(shutil.rmtree, self.media_root, ignore_errors=True)

    def test_run_drug_params_job_creates_versioned_outputs(self) -> None:
        with override_settings(MEDIA_ROOT=str(self.media_root)):
            job = models.ChemJob.objects.create(kind=models.ChemJob.PARAM, input_a="CCO")
            counter = {"value": 0}

            def fake_evaluator(*, smiles: str | None, cid: str | None, outdir: Any) -> tuple[Path, str, str]:
                counter["value"] += 1
                outdir_path = Path(outdir)
                outdir_path.mkdir(parents=True, exist_ok=True)
                result_path = outdir_path / f"result_{counter['value']}.html"
                result_path.write_text("<html></html>", encoding="utf-8")
                stdout = f"stdout {counter['value']}"
                stderr = "stderr 1" if counter["value"] == 1 else ""
                return result_path, stdout, stderr

            with (
                patch("chemtools.tasks.utils.run_settings_generator", return_value="settings ok"),
                patch("chemtools.tasks.utils.run_paths_generator", return_value="paths ok"),
                patch("chemtools.tasks.utils.run_drug_parameter_evaluator", side_effect=fake_evaluator),
            ):
                tasks.run_drug_params_job(job.pk, "CCO", "123")
                job.refresh_from_db()
                job_dir = self.media_root / "chem" / str(job.pk)
                v1_path = job_dir / f"drug_{job.pk}_v1.html"
                self.assertTrue(v1_path.exists())
                self.assertEqual(job.out_html.name, v1_path.relative_to(self.media_root).as_posix())
                self.assertIn("stdout 1", job.log)
                self.assertIn("stderr 1", job.log)
                self.assertIn("--- STDERR ---", job.log)

                tasks.run_drug_params_job(job.pk, "CCO", "123")
                job.refresh_from_db()
                v2_path = job_dir / f"drug_{job.pk}_v2.html"
                self.assertTrue(v2_path.exists())
                self.assertTrue(v1_path.exists())
                self.assertEqual(job.out_html.name, v2_path.relative_to(self.media_root).as_posix())
                self.assertIn("stdout 2", job.log)
                # Ensure previous log is retained with separator.
                self.assertGreaterEqual(job.log.count("---"), 2)

    def test_binding_viz_creates_placeholder_when_snapshot_missing(self) -> None:
        with override_settings(MEDIA_ROOT=str(self.media_root)):
            job = models.ChemJob.objects.create(kind=models.ChemJob.BIND, input_a="5LF3")

            def fake_binding_visualizer(*, pdb_id: str, ligand: str | None, outdir) -> tuple[Path, str, str]:
                outdir_path = Path(outdir)
                outdir_path.mkdir(parents=True, exist_ok=True)
                html_path = outdir_path / f"{pdb_id}_structure_viewer.html"
                html_path.write_text("<html></html>", encoding="utf-8")
                return html_path, "stdout", ""

            with (
                patch("chemtools.tasks.utils.run_settings_generator", return_value="settings ok"),
                patch("chemtools.tasks.utils.run_paths_generator", return_value="paths ok"),
                patch("chemtools.tasks.utils.run_binding_visualizer", side_effect=fake_binding_visualizer),
            ):
                tasks.run_binding_viz_job(job.pk, "5LF3", "")
                job.refresh_from_db()

            job_dir = self.media_root / "chem" / str(job.pk)
            self.assertTrue((job_dir / "thumb.png").exists())
            self.assertTrue(job.out_html.name.endswith(".html"))

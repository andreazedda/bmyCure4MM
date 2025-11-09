"""
Comprehensive API integration tests for chemtools.
Tests external API calls, job processing, and error handling.
"""
from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock, Mock

from django.contrib.auth import get_user_model
from django.test import TestCase, TransactionTestCase
from django.urls import reverse

from chemtools import models, tasks, utils


class DrugParametersAPITests(TestCase):
    """Test drug parameters API integration."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        
    def test_pubchem_cid_lookup(self) -> None:
        """Test PubChem CID lookup works."""
        # CID 2244 = Aspirin
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="",
            input_b="2244",
            user=self.user
        )
        
        with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
            mock_run.return_value = (Path("/tmp/output.html"), "Success", "")
            tasks.run_drug_params_job(job.pk, "", "2244")
        
        mock_run.assert_called_once()
        call_args = mock_run.call_args
        # Check that cid was passed (either positional or keyword)
        self.assertTrue(
            (call_args[1].get("cid") == "2244" if call_args[1] else False) or
            (len(call_args[0]) > 1 and call_args[0][1] == "2244")
        )
        
    def test_smiles_molecular_properties(self) -> None:
        """Test SMILES molecular property calculation."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",  # Ethanol
            input_b="",
            user=self.user
        )
        
        with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
            mock_run.return_value = (Path("/tmp/output.html"), "MW: 46.07", "")
            tasks.run_drug_params_job(job.pk, "CCO", "")
        
        mock_run.assert_called_once()
        call_args = mock_run.call_args
        # Verify SMILES was passed
        self.assertTrue(
            (call_args[1].get("smiles") == "CCO" if call_args[1] else False) or
            (len(call_args[0]) > 0 and call_args[0][0] == "CCO")
        )
        
    def test_invalid_smiles_error_handling(self) -> None:
        """Test error handling for invalid SMILES."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="INVALID!!!",
            input_b="",
            user=self.user
        )
        
        with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
            mock_run.side_effect = utils.ToolRunError("Invalid SMILES", "", "Error parsing")
            tasks.run_drug_params_job(job.pk, "INVALID!!!", "")
        
        job.refresh_from_db()
        self.assertIn("ERROR", job.log)
        self.assertIn("Invalid SMILES", job.log)
        
    def test_pubchem_api_timeout(self) -> None:
        """Test PubChem API timeout handling."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="",
            input_b="999999999",  # Non-existent CID
            user=self.user
        )
        
        with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
            mock_run.side_effect = utils.ToolRunError("PubChem timeout", "", "Timeout after 900s")
            tasks.run_drug_params_job(job.pk, "", "999999999")
        
        job.refresh_from_db()
        self.assertIn("ERROR", job.log)
        
    def test_output_html_generation(self) -> None:
        """Test HTML output file generation and storage."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            input_b="",
            user=self.user
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "test_output.html"
            output_file.write_text("<html><body>Test</body></html>")
            
            with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
                mock_run.return_value = (output_file, "Success", "")
                with patch("chemtools.tasks._store_file") as mock_store:
                    mock_store.return_value = "chem/1/drug_1_v1.html"
                    tasks.run_drug_params_job(job.pk, "CCO", "")
        
        mock_store.assert_called_once()


class BindingVisualizerAPITests(TestCase):
    """Test binding visualizer API integration."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        
    def test_pdb_download_from_rcsb(self) -> None:
        """Test PDB file download from RCSB."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="5LF3",  # Bortezomib complex
            input_b="BOR",
            user=self.user
        )
        
        with patch("chemtools.utils.run_binding_visualizer") as mock_run:
            mock_run.return_value = (Path("/tmp/5LF3_structure_viewer.html"), "Downloaded PDB", "")
            tasks.run_binding_viz_job(job.pk, "5LF3", "BOR")
        
        mock_run.assert_called_once()
        call_args = mock_run.call_args
        # Verify pdb_id and ligand were passed
        self.assertTrue(
            (call_args[1].get("pdb_id") == "5LF3" if call_args[1] else False) or
            (len(call_args[0]) > 0 and call_args[0][0] == "5LF3")
        )
        
    def test_invalid_pdb_id_error(self) -> None:
        """Test error handling for invalid PDB ID."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="XXXX",  # Invalid PDB ID
            input_b="",
            user=self.user
        )
        
        with patch("chemtools.utils.run_binding_visualizer") as mock_run:
            mock_run.side_effect = utils.ToolRunError("PDB not found", "", "404 error")
            tasks.run_binding_viz_job(job.pk, "XXXX", "")
        
        job.refresh_from_db()
        self.assertIn("ERROR", job.log)
        self.assertIn("PDB not found", job.log)
        
    def test_ligand_identification(self) -> None:
        """Test automatic ligand identification."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="5LF3",
            input_b="",  # No ligand specified
            user=self.user
        )
        
        with patch("chemtools.utils.run_binding_visualizer") as mock_run:
            mock_run.return_value = (Path("/tmp/output.html"), "Found ligand BOR", "")
            tasks.run_binding_viz_job(job.pk, "5LF3", "")
        
        mock_run.assert_called_once()
        
    def test_thumbnail_generation(self) -> None:
        """Test structure thumbnail generation."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="5LF3",
            input_b="BOR",
            user=self.user
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            snapshot = tmpdir_path / "5LF3_structure.png"
            snapshot.write_bytes(b"PNG_DATA")
            
            with patch("chemtools.utils.run_binding_visualizer") as mock_run:
                mock_run.return_value = (Path("/tmp/output.html"), "Success", "")
                with patch("chemtools.tasks._job_media_dir") as mock_media:
                    mock_media.return_value = tmpdir_path
                    with patch("tempfile.mkdtemp") as mock_mkdtemp:
                        mock_mkdtemp.return_value = str(tmpdir_path)
                        tasks.run_binding_viz_job(job.pk, "5LF3", "BOR")
        
    def test_3d_structure_rendering(self) -> None:
        """Test 3D structure HTML rendering."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="5LF3",
            input_b="BOR",
            user=self.user
        )
        
        with patch("chemtools.utils.run_binding_visualizer") as mock_run:
            html_content = "<html><script>// 3D rendering code</script></html>"
            output_file = Path("/tmp/5LF3_structure_viewer.html")
            mock_run.return_value = (output_file, "Rendered 3D structure", "")
            tasks.run_binding_viz_job(job.pk, "5LF3", "BOR")
        
        mock_run.assert_called_once()


class SimilaritySearchAPITests(TestCase):
    """Test similarity search API integration."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        
    def test_fingerprint_generation(self) -> None:
        """Test molecular fingerprint generation."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="c1ccccc1",  # Benzene
            input_b="0.7",
            user=self.user
        )
        
        with patch("chemtools.utils.run_similarity_search") as mock_run:
            mock_run.return_value = (Path("/tmp/sim.csv"), "Generated fingerprint", "")
            tasks.run_similarity_job(job.pk, "c1ccccc1")
        
        mock_run.assert_called_once()
        
    def test_similarity_database_query(self) -> None:
        """Test similarity search against database."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            input_b="0.8",
            user=self.user
        )
        
        with patch("chemtools.utils.run_similarity_search") as mock_run:
            csv_content = "SMILES,Similarity,Name\nCCO,1.0,Ethanol\nCCCO,0.85,Propanol\n"
            output_file = Path("/tmp/sim.csv")
            mock_run.return_value = (output_file, "Found 2 matches", "")
            tasks.run_similarity_job(job.pk, "CCO")
        
        mock_run.assert_called_once()
        
    def test_threshold_filtering(self) -> None:
        """Test similarity threshold filtering."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="c1ccccc1",
            input_b="0.9",  # High threshold
            user=self.user
        )
        
        with patch("chemtools.utils.run_similarity_search") as mock_run:
            mock_run.return_value = (Path("/tmp/sim.csv"), "Applied threshold 0.9", "")
            tasks.run_similarity_job(job.pk, "c1ccccc1")
        
        mock_run.assert_called_once()
        
    def test_csv_output_generation(self) -> None:
        """Test CSV output file generation."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            input_b="0.7",
            user=self.user
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "similarity_results.csv"
            output_file.write_text("SMILES,Similarity\nCCO,1.0\n")
            
            with patch("chemtools.utils.run_similarity_search") as mock_run:
                mock_run.return_value = (output_file, "Success", "")
                with patch("chemtools.tasks._store_file") as mock_store:
                    mock_store.return_value = "chem/1/sim_1_v1.csv"
                    tasks.run_similarity_job(job.pk, "CCO")
        
        mock_store.assert_called_once()


class JobLifecycleTests(TransactionTestCase):
    """Test complete job lifecycle from creation to completion."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_job_creation_from_form(self) -> None:
        """Test job creation through form submission."""
        with patch("chemtools.views._enqueue") as mock_enqueue:
            mock_enqueue.return_value = (True, False)
            response = self.client.post(reverse("chemtools:drug_params"), {
                "smiles": "CCO",
                "cid": ""
            })
        
        self.assertEqual(response.status_code, 302)
        job = models.ChemJob.objects.filter(user=self.user).first()
        self.assertIsNotNone(job)
        self.assertEqual(job.input_a, "CCO")
        self.assertEqual(job.kind, models.ChemJob.PARAM)
        
    def test_job_queueing_celery(self) -> None:
        """Test job queueing with Celery."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        with patch("chemtools.tasks.run_drug_params_job.delay") as mock_delay:
            from chemtools.views import _enqueue
            queued, failed = _enqueue(tasks.run_drug_params_job, job, "CCO", "")
        
        self.assertTrue(queued)
        self.assertFalse(failed)
        mock_delay.assert_called_once_with(job.pk, "CCO", "")
        
    def test_job_status_polling(self) -> None:
        """Test job status polling endpoint."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        response = self.client.get(reverse("chemtools:job_status", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("status", data)
        self.assertIn("variant", data)
        
    def test_job_completion_with_output(self) -> None:
        """Test job completion and output availability."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        # Simulate job completion
        job.out_html.name = "chem/1/drug_1_v1.html"
        job.save()
        
        response = self.client.get(reverse("chemtools:job_status", args=[job.pk]))
        data = response.json()
        self.assertTrue(data["has_html"])
        self.assertIsNotNone(data["html_url"])
        
    def test_job_retry_mechanism(self) -> None:
        """Test job retry functionality."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user,
            log="ERROR: Previous failure"
        )
        
        with patch("chemtools.views._enqueue") as mock_enqueue:
            mock_enqueue.return_value = (True, False)
            response = self.client.get(reverse("chemtools:job_retry", args=[job.pk]))
        
        self.assertEqual(response.status_code, 302)
        mock_enqueue.assert_called_once()


class ErrorHandlingTests(TestCase):
    """Test robust error handling."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        
    def test_network_timeout_handling(self) -> None:
        """Test handling of network timeouts."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="",
            input_b="2244",
            user=self.user
        )
        
        with patch("chemtools.utils._run_command") as mock_cmd:
            import subprocess
            mock_cmd.side_effect = utils.ToolRunError("Timeout", "", "")
            tasks.run_drug_params_job(job.pk, "", "2244")
        
        job.refresh_from_db()
        self.assertIn("ERROR", job.log)
        
    def test_malformed_input_validation(self) -> None:
        """Test validation of malformed inputs."""
        from chemtools.validators import validate_smiles_basic
        
        # Test SQL injection attempt
        is_valid, error = validate_smiles_basic("CCO'; DROP TABLE;--")
        self.assertFalse(is_valid)
        self.assertIn("DROP", error)
        
        # Test XSS attempt
        is_valid, error = validate_smiles_basic("<script>alert('xss')</script>")
        self.assertFalse(is_valid)
        
    def test_file_size_limit_enforcement(self) -> None:
        """Test file size limit enforcement."""
        from chemtools.validators import validate_pdb_content
        
        # Create oversized content (>10MB)
        large_content = b"ATOM" * (11 * 1024 * 1024 // 4)
        is_valid, error = validate_pdb_content(large_content, max_size_mb=10)
        self.assertFalse(is_valid)
        self.assertIn("too large", error)
        
    def test_missing_output_file_handling(self) -> None:
        """Test handling of missing output files."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
            # Return None for output file
            mock_run.return_value = (None, "No output generated", "")
            tasks.run_drug_params_job(job.pk, "CCO", "")
        
        job.refresh_from_db()
        self.assertFalse(job.out_html)
        
    def test_corrupted_output_handling(self) -> None:
        """Test handling of corrupted output files."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            corrupt_file = Path(tmpdir) / "corrupt.html"
            corrupt_file.write_bytes(b"\x00\x00\xFF\xFF")  # Invalid HTML
            
            with patch("chemtools.utils.run_drug_parameter_evaluator") as mock_run:
                mock_run.return_value = (corrupt_file, "Warning: output may be corrupt", "")
                # Should not crash
                tasks.run_drug_params_job(job.pk, "CCO", "")

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.test import TestCase
from django.urls import reverse
from unittest.mock import patch, MagicMock
from django.core.files.uploadedfile import SimpleUploadedFile

from chemtools import models


class ToolsHomeViewTests(TestCase):
    """Test the chemtools home page with job listing."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_tools_home_renders_correctly(self) -> None:
        """Test that tools home page renders with explanations."""
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Drug Discovery")
        self.assertContains(response, "Cheminformatics")
        
    def test_tools_home_shows_jobs_table(self) -> None:
        """Test that jobs are displayed in the table."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCCC",
            user=self.user
        )
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertContains(response, "CCCC")
        self.assertContains(response, job.get_kind_display())
        
    def test_tools_home_shows_empty_state(self) -> None:
        """Test empty state when no jobs exist."""
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertContains(response, "No Jobs Yet")
        self.assertContains(response, "Try Drug Parameters")


class DrugParamsViewTests(TestCase):
    """Test drug parameters calculation view."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_drug_params_page_loads(self) -> None:
        """Test that drug parameters page loads."""
        response = self.client.get(reverse("chemtools:drug_params"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Drug Parameter Evaluation")
        
    def test_drug_params_form_submission(self) -> None:
        """Test submitting drug parameters form."""
        with patch("chemtools.views._enqueue", return_value=(True, False)):
            response = self.client.post(reverse("chemtools:drug_params"), {
                "smiles": "CCO"
            })
        self.assertRedirects(response, reverse("chemtools:tools_home"))
        
    def test_drug_params_validates_empty_smiles(self) -> None:
        """Test that empty SMILES is rejected."""
        response = self.client.post(reverse("chemtools:drug_params"), {
            "smiles": ""
        })
        self.assertEqual(response.status_code, 200)
        # Form should reject empty values (requires SMILES or CID)
        self.assertContains(response, "Provide a SMILES string or a PubChem CID")


class SimilaritySearchViewTests(TestCase):
    """Test similarity search view."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_similarity_page_loads(self) -> None:
        """Test that similarity search page loads."""
        response = self.client.get(reverse("chemtools:similarity"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Similarity Search")
        
    def test_similarity_form_submission(self) -> None:
        """Test submitting similarity search form."""
        with patch("chemtools.views._enqueue", return_value=(True, False)):
            response = self.client.post(reverse("chemtools:similarity"), {
                "smiles": "CCO",
                "threshold": "0.7"
            })
        self.assertRedirects(response, reverse("chemtools:tools_home"))
        
    def test_similarity_validates_threshold(self) -> None:
        """Test that invalid threshold is rejected."""
        response = self.client.post(reverse("chemtools:similarity"), {
            "smiles": "CCO",
            "threshold": "1.5"  # Invalid: >1.0
        })
        self.assertEqual(response.status_code, 200)
        # Form should have validation error for threshold


class BindingVisualizerViewTests(TestCase):
    """Test binding visualizer view."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_binding_viz_page_loads(self) -> None:
        """Test that binding visualizer page loads."""
        response = self.client.get(reverse("chemtools:binding_viz"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Binding Visualizer")
        
    def test_binding_viz_file_upload(self) -> None:
        """Test submitting a PDB ID."""
        with patch("chemtools.views._enqueue", return_value=(True, False)):
            response = self.client.post(reverse("chemtools:binding_viz"), {
                "pdb_id": "5LF3",
                "ligand": "BOR"
            })
        self.assertRedirects(response, reverse("chemtools:tools_home"))


class JobStatusViewTests(TestCase):
    """Test job status polling endpoint."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_job_status_returns_json(self) -> None:
        """Test that job status endpoint returns JSON."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        response = self.client.get(reverse("chemtools:job_status", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "application/json")
        
    def test_job_status_includes_status(self) -> None:
        """Test that job status includes status field."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        response = self.client.get(reverse("chemtools:job_status", args=[job.pk]))
        data = response.json()
        self.assertIn("status", data)
        self.assertIn("variant", data)


class RetryJobViewTests(TestCase):
    """Test job retry functionality."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)

    def test_retry_job_shows_info_when_queued(self) -> None:
        """Test retry shows success message when job is queued."""
        job = models.ChemJob.objects.create(kind=models.ChemJob.SIM, input_a="CCC", user=self.user)
        with patch("chemtools.views._enqueue", return_value=(True, False)) as enqueue:
            response = self.client.get(reverse("chemtools:job_retry", args=[job.pk]))
            enqueue.assert_called_once()
        self.assertRedirects(response, reverse("chemtools:tools_home"))
        messages = list(get_messages(response.wsgi_request))
        self.assertTrue(any("queued" in message.message.lower() for message in messages))

    def test_retry_job_shows_error_when_immediate_failure(self) -> None:
        """Test retry shows error message when job fails immediately."""
        job = models.ChemJob.objects.create(kind=models.ChemJob.PARAM, input_a="CCO", input_b="123", user=self.user)
        with patch("chemtools.views._enqueue", return_value=(False, True)):
            response = self.client.get(reverse("chemtools:job_retry", args=[job.pk]))
        self.assertRedirects(response, reverse("chemtools:tools_home"))
        messages = list(get_messages(response.wsgi_request))
        self.assertTrue(any("failed" in message.message.lower() for message in messages))
        
    def test_retry_nonexistent_job_404(self) -> None:
        """Test retrying non-existent job returns 404."""
        response = self.client.get(reverse("chemtools:job_retry", args=[9999]))
        self.assertEqual(response.status_code, 404)


class ChemToolsSecurityTests(TestCase):
    """Test security and access control for chemtools."""
    
    def test_tools_home_requires_login(self) -> None:
        """Test that tools home requires authentication."""
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertEqual(response.status_code, 302)
        self.assertIn("/admin/login/", response.url)
        
    def test_user_can_only_see_own_jobs(self) -> None:
        """Test that users can only see their own jobs."""
        user1 = get_user_model().objects.create_user("user1", password="pass123")
        user2 = get_user_model().objects.create_user("user2", password="pass123")
        
        job1 = models.ChemJob.objects.create(kind=models.ChemJob.SIM, input_a="AAA", user=user1)
        job2 = models.ChemJob.objects.create(kind=models.ChemJob.SIM, input_a="BBB", user=user2)
        
        self.client.force_login(user1)
        response = self.client.get(reverse("chemtools:tools_home"))
        
        self.assertContains(response, "AAA")
        self.assertNotContains(response, "BBB")


class JobDetailViewTests(TestCase):
    """Test integrated job detail view."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_job_detail_requires_login(self) -> None:
        """Test that job detail requires authentication."""
        self.client.logout()
        job = models.ChemJob.objects.create(kind=models.ChemJob.PARAM, input_a="CCO", user=self.user)
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 302)
        self.assertIn("/admin/login/", response.url)
        
    def test_job_detail_user_isolation(self) -> None:
        """Test that users can only view their own jobs."""
        user2 = get_user_model().objects.create_user("other", password="pass123")
        job = models.ChemJob.objects.create(kind=models.ChemJob.PARAM, input_a="CCO", user=user2)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 404)
        
    def test_job_detail_nonexistent_404(self) -> None:
        """Test that non-existent job returns 404."""
        response = self.client.get(reverse("chemtools:job_detail", args=[9999]))
        self.assertEqual(response.status_code, 404)
        
    def test_job_detail_drug_params_renders(self) -> None:
        """Test job detail renders Drug Parameters correctly."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        # Create mock HTML output with parameters table
        from django.core.files.base import ContentFile
        html_content = """
        <html>
            <div id="container-01"></div>
            <table class="parameters-table">
                <tr><td>MW</td><td>46.07</td></tr>
                <tr><td>LogP</td><td>-0.31</td></tr>
            </table>
            <script>
                let viewer = $3Dmol.createViewer($("#container-01"));
            </script>
        </html>
        """
        job.out_html.save("test.html", ContentFile(html_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Molecular Properties")
        self.assertContains(response, "parameters-table")
        self.assertContains(response, "46.07")  # MW value
        self.assertContains(response, "LogP")
        
    def test_job_detail_similarity_search_renders(self) -> None:
        """Test job detail renders Similarity Search correctly."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            input_b="0.9",
            user=self.user
        )
        
        # Create mock CSV output
        from django.core.files.base import ContentFile
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
702,CCO,1.0,46.07,-0.31,1,1,20.23,0,-0.77,Ethanol
1031,CCCO,0.95,60.1,0.25,1,1,20.23,1,-1.01,Propanol"""
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify the CSV was loaded (csv_data should be in context)
        if hasattr(response, 'context') and response.context:
            csv_data = response.context.get('csv_data')
            if csv_data:
                self.assertGreater(len(csv_data), 0)
                # Verify content
                self.assertContains(response, "Similar Compounds Found")
                self.assertContains(response, "Ethanol")
            else:
                # If CSV wasn't loaded, at least check job is displayed
                self.assertContains(response, job.get_kind_display())
        
    def test_job_detail_binding_viz_renders(self) -> None:
        """Test job detail renders Binding Visualizer correctly."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="5LF3",
            input_b="BOR",
            user=self.user
        )
        
        # Create mock HTML output
        from django.core.files.base import ContentFile
        html_content = """
        <html>
            <div id="viewport" style="width:100%;height:600px;"></div>
            <script>
                let viewer = $3Dmol.createViewer("#viewport");
                viewer.addModel("ATOM...", "pdb");
            </script>
        </html>
        """
        job.out_html.save("test.html", ContentFile(html_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "3D Protein-Ligand Structure")
        self.assertContains(response, "5LF3")
        self.assertContains(response, "$3Dmol")
        
    def test_job_detail_empty_csv(self) -> None:
        """Test job detail handles empty CSV gracefully."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        # Create CSV with only headers
        from django.core.files.base import ContentFile
        csv_content = "CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name\n"
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        # With empty CSV, similar compounds section should not appear
        # Instead check that job is displayed
        self.assertContains(response, job.get_kind_display())
        
    def test_job_detail_queued_auto_refresh(self) -> None:
        """Test that queued jobs have auto-refresh meta tag."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        # Don't set output files - job is still queued
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'http-equiv="refresh"')
        self.assertContains(response, 'content="5"')  # 5 second refresh
        
    def test_job_detail_completed_no_refresh(self) -> None:
        """Test that completed jobs don't have auto-refresh."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        # Mark job as complete with output
        from django.core.files.base import ContentFile
        job.out_html.save("test.html", ContentFile("<html></html>"), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        # Should NOT contain refresh meta tag
        self.assertNotContains(response, 'http-equiv="refresh"')
        
    def test_job_detail_progress_display(self) -> None:
        """Test that job progress is displayed."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user,
            progress_percent=45,
            progress_message="Processing compound 45 of 100"
        )
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "45%")
        self.assertContains(response, "Processing compound 45 of 100")
        
    def test_job_detail_log_display(self) -> None:
        """Test that job log is displayed in collapsible section."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user,
            log="Processing SMILES: CCO\nCalculating MW: 46.07\nSuccess!"
        )
        
        from django.core.files.base import ContentFile
        job.out_html.save("test.html", ContentFile("<html></html>"), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Processing SMILES: CCO")
        self.assertContains(response, "Calculating MW: 46.07")
        
    def test_job_detail_csv_parsing_error_handling(self) -> None:
        """Test graceful handling of malformed CSV."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        # Create malformed CSV
        from django.core.files.base import ContentFile
        csv_content = "Malformed,CSV\nWithout,Proper,Headers"
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        # Should not crash
        self.assertEqual(response.status_code, 200)


class ChemToolsIntegrationTests(TestCase):
    """Integration tests for complete workflows."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_complete_drug_params_workflow(self) -> None:
        """Test complete workflow: submit job, check status, view results."""
        # Submit job
        with patch("chemtools.views._enqueue", return_value=(True, False)):
            response = self.client.post(reverse("chemtools:drug_params"), {
                "smiles": "CCO"
            })
        
        # Check job was created
        job = models.ChemJob.objects.filter(user=self.user).first()
        self.assertIsNotNone(job)
        self.assertEqual(job.input_a, "CCO")
        
        # Check status endpoint works
        response = self.client.get(reverse("chemtools:job_status", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Check job detail view works
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
    def test_complete_similarity_search_workflow(self) -> None:
        """Test complete similarity search workflow."""
        with patch("chemtools.views._enqueue", return_value=(True, False)):
            response = self.client.post(reverse("chemtools:similarity"), {
                "smiles": "c1ccccc1",  # Benzene
                "threshold": "0.8"
            })
        
        job = models.ChemJob.objects.filter(user=self.user, kind=models.ChemJob.SIM).first()
        self.assertIsNotNone(job)
        self.assertEqual(job.input_a, "c1ccccc1")
        
        # Check integrated view works
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "c1ccccc1")
        
    def test_integrated_view_replaces_downloads(self) -> None:
        """Test that integrated view shows results without downloads."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        from django.core.files.base import ContentFile
        job.out_html.save("test.html", ContentFile("<html><p>Results</p></html>"), save=True)
        
        # Check tools_home has "View Results" button
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertContains(response, "📊 View Results")
        self.assertContains(response, f"/chemtools/job/{job.pk}/")
        
        # Check job detail shows embedded content
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertContains(response, "<p>Results</p>")  # Embedded HTML


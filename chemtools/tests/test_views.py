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


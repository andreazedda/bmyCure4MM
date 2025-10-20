"""
Tests for virtual cohort simulation module.
"""
from django.contrib.auth import get_user_model
from django.test import TestCase, Client
from django.urls import reverse

from simulator import cohort, models

User = get_user_model()


class CohortSamplingTestCase(TestCase):
    """Test patient parameter sampling."""
    
    def test_sample_patient_params_generates_n_patients(self):
        """Test that sampling generates exactly n patients."""
        n = 50
        patients = cohort.sample_patient_params(seed=42, n=n)
        
        self.assertEqual(len(patients), n)
        
        # Check structure
        for patient in patients:
            self.assertIn("patient_id", patient)
            self.assertIn("baseline_tumor_cells", patient)
            self.assertIn("baseline_healthy_cells", patient)
            self.assertIn("tumor_growth_rate", patient)
            self.assertIn("healthy_growth_rate", patient)
    
    def test_sample_patient_params_is_reproducible(self):
        """Test that same seed produces same results."""
        patients1 = cohort.sample_patient_params(seed=42, n=10)
        patients2 = cohort.sample_patient_params(seed=42, n=10)
        
        for p1, p2 in zip(patients1, patients2):
            self.assertAlmostEqual(
                p1["baseline_tumor_cells"],
                p2["baseline_tumor_cells"],
                places=6
            )
    
    def test_sample_patient_params_variability(self):
        """Test that different patients have different parameters."""
        patients = cohort.sample_patient_params(seed=42, n=100)
        
        tumor_baselines = [p["baseline_tumor_cells"] for p in patients]
        
        # Check that values are not all identical
        self.assertGreater(len(set(tumor_baselines)), 1)
        
        # Check that values are positive
        self.assertTrue(all(t > 0 for t in tumor_baselines))


class CohortSimulationTestCase(TestCase):
    """Test cohort simulation runner."""
    
    def setUp(self):
        """Create test user and scenario."""
        self.user = User.objects.create_user(
            username="researcher",
            email="research@example.com",
            password="testpass123"
        )
        
        self.scenario = models.Scenario.objects.create(
            title="Test Cohort Scenario",
            summary="Test patient for cohort simulation",
            clinical_stage="newly_diagnosed",
        )
    
    def test_run_cohort_returns_aggregates(self):
        """Test that cohort simulation returns aggregate statistics."""
        regimen_params = {
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "time_horizon": 168,
            "interaction_strength": 0.1,
            "len_on_days": 21,
            "bor_weekly": 0,
            "dara_interval": 7,
        }
        
        result = cohort.run_cohort(
            scenario=self.scenario,
            n=20,  # Small for fast testing
            regimen_params=regimen_params,
            user_id=self.user.id,
            seed=42,
        )
        
        # Check structure
        self.assertIn("cohort_id", result)
        self.assertIn("n", result)
        self.assertIn("summaries", result)
        self.assertIn("aggregates", result)
        self.assertEqual(result["n"], 20)
        self.assertEqual(len(result["summaries"]), 20)
        
        # Check aggregates
        agg = result["aggregates"]
        self.assertIn("efficacy_mean", agg)
        self.assertIn("efficacy_p95", agg)
        self.assertIn("toxicity_mean", agg)
        self.assertIn("recurrence_rate", agg)
        
        # Check that aggregates are numeric
        self.assertIsInstance(agg["efficacy_mean"], float)
        self.assertIsInstance(agg["toxicity_mean"], float)


class CohortViewTestCase(TestCase):
    """Test cohort views and endpoints."""
    
    def setUp(self):
        """Create test user and scenario."""
        self.user = User.objects.create_user(
            username="tester",
            email="test@example.com",
            password="testpass123"
        )
        self.user.is_staff = True
        self.user.save()
        
        self.scenario = models.Scenario.objects.create(
            title="View Test Scenario",
            summary="Test",
            clinical_stage="newly_diagnosed",
            active=True,
        )
        
        self.client = Client()
        self.client.login(username="tester", password="testpass123")
    
    def test_cohort_view_requires_login(self):
        """Test that cohort view requires authentication."""
        self.client.logout()
        response = self.client.get(reverse("simulator:cohort"))
        self.assertEqual(response.status_code, 302)  # Redirect to login
    
    def test_cohort_view_displays_form(self):
        """Test that cohort view displays configuration form."""
        response = self.client.get(reverse("simulator:cohort"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Virtual Cohort")
        self.assertContains(response, "n_patients")
    
    def test_cohort_run_with_valid_params(self):
        """Test running cohort with valid parameters."""
        response = self.client.post(
            reverse("simulator:cohort_run"),
            {
                "scenario_id": self.scenario.pk,
                "n_patients": 10,
                "seed": 42,
                "lenalidomide_dose": 25.0,
                "bortezomib_dose": 1.3,
                "daratumumab_dose": 16.0,
                "time_horizon": 168,
                "interaction_strength": 0.1,
                "len_on_days": 21,
                "bor_weekly": 0,
                "dara_interval": 7,
            },
        )
        
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Cohort Results")
        self.assertContains(response, "Efficacy")
    
    def test_cohort_run_rejects_invalid_n_patients(self):
        """Test that invalid n_patients is rejected."""
        response = self.client.post(
            reverse("simulator:cohort_run"),
            {
                "scenario_id": self.scenario.pk,
                "n_patients": 2000,  # Too large
                "seed": 42,
                "lenalidomide_dose": 25.0,
                "bortezomib_dose": 1.3,
                "daratumumab_dose": 16.0,
                "time_horizon": 168,
                "interaction_strength": 0.1,
                "len_on_days": 21,
                "bor_weekly": 0,
                "dara_interval": 7,
            },
        )
        
        self.assertEqual(response.status_code, 400)
        self.assertIn(b"n_patients must be between 10 and 1000", response.content)
    
    def test_cohort_export_returns_csv(self):
        """Test CSV export endpoint."""
        # First run a cohort
        result = cohort.run_cohort(
            scenario=self.scenario,
            n=10,
            regimen_params={
                "lenalidomide_dose": 25.0,
                "bortezomib_dose": 1.3,
                "daratumumab_dose": 16.0,
                "time_horizon": 168,
                "interaction_strength": 0.1,
                "len_on_days": 21,
                "bor_weekly": 0,
                "dara_interval": 7,
            },
            user_id=self.user.id,
            seed=42,
        )
        
        # Cache result manually (since we're bypassing the view)
        from django.core.cache import cache
        cache.set(f"cohort_result_{result['cohort_id']}", result, 1800)
        
        # Export
        response = self.client.get(
            reverse("simulator:cohort_export", args=[result["cohort_id"]])
        )
        
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "text/csv")
        self.assertIn(b"patient_id", response.content)
        self.assertIn(b"tumor_reduction", response.content)
        self.assertIn(b"SUMMARY STATISTICS", response.content)

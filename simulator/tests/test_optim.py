"""
Tests for simulator optimization module (optim.py).

Validates multi-objective optimization with Optuna MOTPE sampler.
"""
from django.test import TestCase
from django.contrib.auth import get_user_model

from simulator import models, optim

User = get_user_model()


class OptimizationTestCase(TestCase):
    """Test Optuna multi-objective optimization."""
    
    def setUp(self):
        """Create test user and scenario."""
        self.user = User.objects.create_user(
            username="researcher",
            email="research@example.com",
            password="testpass123"
        )
        
        # Create simple scenario - model only has title, summary, clinical_stage
        self.scenario = models.Scenario.objects.create(
            title="Test Scenario for Optimization",
            summary="Test patient with MM for optimization experiments",
            clinical_stage="newly_diagnosed",
        )
    
    def test_run_study_returns_pareto_frontier(self):
        """
        Test that run_study() executes and returns Pareto frontier.
        
        Validates:
        - Study completes without errors
        - Result contains 'pareto' key
        - Pareto solutions have expected structure
        """
        result = optim.run_study(
            scenario=self.scenario,
            user_id=self.user.id,
            n_trials=20,  # Small for fast testing
            seed=42
        )
        
        # Check result structure
        self.assertIn("pareto", result)
        self.assertIn("study", result)
        self.assertIsInstance(result["pareto"], list)
        
        # Check Pareto solutions structure (if any found)
        if result["pareto"]:
            solution = result["pareto"][0]
            self.assertIn("efficacy", solution)
            self.assertIn("safety", solution)
            self.assertIn("exposure", solution)
            self.assertIn("params", solution)
            
            # Validate that key optimizer parameters are present
            params = solution["params"]
            expected_keys = {
                "lenalidomide_dose",
                "bortezomib_dose",
                "daratumumab_dose",
                "time_horizon",
            }
            self.assertTrue(expected_keys.issubset(params.keys()))
    
    def test_run_study_respects_constraints(self):
        """
        Test that optimization respects toxicity constraints.
        
        All Pareto solutions should satisfy:
        - healthy_loss <= 0.25 (25% max toxicity)
        """
        result = optim.run_study(
            scenario=self.scenario,
            user_id=self.user.id,
            n_trials=30,
            seed=42
        )
        
        # Check that all Pareto solutions satisfy toxicity constraint
        for solution in result["pareto"]:
            # Safety = 1 - healthy_loss, so safety >= 0.75 means healthy_loss <= 0.25
            self.assertGreaterEqual(
                solution["safety"],
                0.75,
                msg=f"Solution violates toxicity constraint: {solution}"
            )
    
    def test_run_study_with_different_seeds(self):
        """
        Test that different seeds produce different results.
        
        Validates randomness in optimization process.
        """
        result1 = optim.run_study(
            scenario=self.scenario,
            user_id=self.user.id,
            n_trials=15,
            seed=42
        )
        
        result2 = optim.run_study(
            scenario=self.scenario,
            user_id=self.user.id,
            n_trials=15,
            seed=99
        )
        
        # Different seeds should explore different parts of parameter space
        # (not guaranteed, but highly likely with 15 trials)
        if result1["pareto"] and result2["pareto"]:
            # Check if at least one parameter differs in first solution
            params1 = result1["pareto"][0]["params"]
            params2 = result2["pareto"][0]["params"]
            
            # Allow some tolerance for floating point comparison
            has_difference = any(
                abs(params1.get(k, 0) - params2.get(k, 0)) > 1e-6
                for k in params1.keys()
            )
            self.assertTrue(
                has_difference,
                "Different seeds should explore different parameter space"
            )

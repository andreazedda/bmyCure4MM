from __future__ import annotations

from types import SimpleNamespace

from django.test import SimpleTestCase

from twin_engine.robustness import compute_robust_scenario_ranking


class RobustnessDiagnosticsTests(SimpleTestCase):
    def test_robust_ranking_uses_probability_best_and_rank_stability(self) -> None:
        run_a = SimpleNamespace(
            id=1,
            intervention_definition={"label": "Scenario A"},
            comparison_metrics={
                "research_utility_v2": 1.2,
                "uncertainty": {
                    "status": "completed",
                    "metric_summaries": {"research_utility_v2": {"mean": 1.15, "p05": 1.0, "p95": 1.3, "point_estimate": 1.2}},
                    "samples": [
                        {"sample_index": 0, "metrics": {"research_utility_v2": 1.2}},
                        {"sample_index": 1, "metrics": {"research_utility_v2": 1.1}},
                        {"sample_index": 2, "metrics": {"research_utility_v2": 1.3}},
                    ],
                },
            },
        )
        run_b = SimpleNamespace(
            id=2,
            intervention_definition={"label": "Scenario B"},
            comparison_metrics={
                "research_utility_v2": 1.0,
                "uncertainty": {
                    "status": "completed",
                    "metric_summaries": {"research_utility_v2": {"mean": 0.98, "p05": 0.9, "p95": 1.1, "point_estimate": 1.0}},
                    "samples": [
                        {"sample_index": 0, "metrics": {"research_utility_v2": 0.9}},
                        {"sample_index": 1, "metrics": {"research_utility_v2": 1.05}},
                        {"sample_index": 2, "metrics": {"research_utility_v2": 0.95}},
                    ],
                },
            },
        )
        result = compute_robust_scenario_ranking([run_a, run_b])
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["rows"][0]["scenario_label"], "Scenario A")
        self.assertGreater(result["rows"][0]["probability_best"], 0.5)
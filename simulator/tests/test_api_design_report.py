from __future__ import annotations

from django.urls import reverse


def test_design_report_endpoint_smoke(client):
    url = reverse("api_design_report")
    resp = client.get(url)
    assert resp.status_code == 200
    data = resp.json()

    # DesignReport schema (see simulator/design/types.py)
    assert set(
        [
            "inputs",
            "tradeoff_points",
            "pareto_front",
            "dynamics",
            "relapse_forecast",
            "toxicity_decomposition",
            "rationale",
        ]
    ).issubset(set(data.keys()))

    assert isinstance(data["tradeoff_points"], list)
    assert isinstance(data["dynamics"], list)
    assert isinstance(data["relapse_forecast"], dict)
    assert isinstance(data["toxicity_decomposition"], dict)
    assert isinstance(data["rationale"], dict)


def test_design_report_endpoint_accepts_params(client):
    url = reverse("api_design_report")
    resp = client.get(url, {"seed": 123, "steps": 3, "dt_days": 2.0})
    assert resp.status_code == 200
    data = resp.json()
    assert data["inputs"]["seed"] == 123
    assert data["inputs"]["steps"] == 3
    assert data["inputs"]["dt_days"] == 2.0
    assert len(data["dynamics"]) == 3

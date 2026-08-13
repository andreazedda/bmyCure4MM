from __future__ import annotations

from typing import Any


COMPLETED_STATUS = "completed"
INSUFFICIENT_UNCERTAINTY_STATUS = "insufficient_uncertainty"


def compute_robust_scenario_ranking(runs) -> dict[str, Any]:
    scenario_payloads = [_extract_scenario_payload(run) for run in runs]
    scenario_payloads = [item for item in scenario_payloads if item is not None]
    if len(scenario_payloads) < 2:
        return {
            "status": INSUFFICIENT_UNCERTAINTY_STATUS,
            "message": "At least two scenarios with completed uncertainty summaries are required for robust ranking.",
            "rows": [],
            "pairwise_dominance": {},
        }

    shared_indices = sorted(set.intersection(*(set(item["sample_utilities"]) for item in scenario_payloads)))
    if not shared_indices:
        return {
            "status": INSUFFICIENT_UNCERTAINTY_STATUS,
            "message": "No aligned uncertainty samples were available across scenarios.",
            "rows": [],
            "pairwise_dominance": {},
        }

    wins = {item["label"]: 0.0 for item in scenario_payloads}
    for sample_index in shared_indices:
        sample_values = [(item["label"], item["sample_utilities"][sample_index]) for item in scenario_payloads]
        best_value = max(value for _label, value in sample_values)
        winners = [label for label, value in sample_values if abs(float(value) - float(best_value)) <= 1.0e-9]
        share = 1.0 / max(len(winners), 1)
        for label in winners:
            wins[label] += share

    point_rank_order = sorted(scenario_payloads, key=lambda item: (-item["point_utility_v2"], item["label"]))
    robust_rank_order = sorted(
        scenario_payloads,
        key=lambda item: (
            -(wins[item["label"]] / len(shared_indices)),
            -item["mean_utility_v2"],
            item["label"],
        ),
    )
    point_ranks = {item["label"]: index for index, item in enumerate(point_rank_order, start=1)}
    robust_ranks = {item["label"]: index for index, item in enumerate(robust_rank_order, start=1)}

    pairwise_dominance = {}
    for left in scenario_payloads:
        left_payload = {}
        for right in scenario_payloads:
            if left["label"] == right["label"]:
                continue
            wins_left = 0
            for sample_index in shared_indices:
                if left["sample_utilities"][sample_index] > right["sample_utilities"][sample_index]:
                    wins_left += 1
            left_payload[right["label"]] = wins_left / len(shared_indices)
        pairwise_dominance[left["label"]] = left_payload

    rows = []
    for item in robust_rank_order:
        next_item = robust_rank_order[robust_ranks[item["label"]]] if robust_ranks[item["label"]] < len(robust_rank_order) else None
        overlap_with_next = _interval_overlap(item, next_item) if next_item else 0.0
        probability_best = wins[item["label"]] / len(shared_indices)
        rows.append(
            {
                "scenario_label": item["label"],
                "point_rank": point_ranks[item["label"]],
                "robust_rank": robust_ranks[item["label"]],
                "point_utility_v2": item["point_utility_v2"],
                "mean_utility_v2": item["mean_utility_v2"],
                "p05": item["p05"],
                "p95": item["p95"],
                "probability_best": probability_best,
                "overlap_with_next": overlap_with_next,
                "robustness_classification": _classify_robustness(
                    probability_best=probability_best,
                    point_rank=point_ranks[item["label"]],
                    robust_rank=robust_ranks[item["label"]],
                    overlap_with_next=overlap_with_next,
                ),
            }
        )

    return {
        "status": COMPLETED_STATUS,
        "n_scenarios": len(scenario_payloads),
        "n_aligned_samples": len(shared_indices),
        "rows": rows,
        "pairwise_dominance": pairwise_dominance,
        "unstable_rank_flag": bool(rows and rows[0]["point_rank"] != rows[0]["robust_rank"]),
        "limitations": [
            "Robust ranking depends on the stored uncertainty sample traces and remains an exploratory scenario-ordering diagnostic.",
            "Probability-best summaries do not imply clinical superiority or causal effect estimation.",
        ],
    }


def _extract_scenario_payload(run) -> dict[str, Any] | None:
    comparison_metrics = getattr(run, "comparison_metrics", {}) or {}
    uncertainty = _latest_uncertainty(run) or comparison_metrics.get("uncertainty") or {}
    if uncertainty.get("status") != COMPLETED_STATUS:
        return None
    metric_summary = (uncertainty.get("metric_summaries") or {}).get("research_utility_v2") or {}
    sample_utilities = {}
    for sample in uncertainty.get("samples") or []:
        metrics = sample.get("metrics") or {}
        if metrics.get("research_utility_v2") is not None:
            sample_utilities[int(sample["sample_index"])] = float(metrics["research_utility_v2"])
    if not sample_utilities:
        return None
    label = str((getattr(run, "intervention_definition", {}) or {}).get("label") or f"Run {getattr(run, 'id', 'unknown')}")
    point_utility = comparison_metrics.get("research_utility_v2")
    if point_utility is None:
        point_utility = metric_summary.get("point_estimate")
    return {
        "label": label,
        "point_utility_v2": float(point_utility or 0.0),
        "mean_utility_v2": float(metric_summary.get("mean") or 0.0),
        "p05": metric_summary.get("p05"),
        "p95": metric_summary.get("p95"),
        "sample_utilities": sample_utilities,
    }


def _latest_uncertainty(run) -> dict[str, Any]:
    metadata_records = getattr(run, "metadata_records", None)
    if metadata_records is None:
        return {}
    record = metadata_records.filter(solver_name="counterfactual_uncertainty").order_by("-created_at").first()
    if record is None:
        return {}
    return dict((record.solver_parameters or {}).get("diagnostic_summary") or {})


def _interval_overlap(current: dict[str, Any], next_item: dict[str, Any] | None) -> float:
    if not next_item:
        return 0.0
    left_low = current.get("p05")
    left_high = current.get("p95")
    right_low = next_item.get("p05")
    right_high = next_item.get("p95")
    if None in {left_low, left_high, right_low, right_high}:
        return 0.0
    overlap = max(0.0, min(float(left_high), float(right_high)) - max(float(left_low), float(right_low)))
    width = max(float(left_high) - float(left_low), float(right_high) - float(right_low), 1.0e-9)
    return overlap / width


def _classify_robustness(*, probability_best: float, point_rank: int, robust_rank: int, overlap_with_next: float) -> str:
    if point_rank != robust_rank:
        return "fragile"
    if probability_best >= 0.70 and overlap_with_next <= 0.25:
        return "robust"
    if probability_best >= 0.40:
        return "contested"
    return "fragile"

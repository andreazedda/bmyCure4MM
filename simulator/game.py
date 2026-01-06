from __future__ import annotations


def compute_game_metrics(summary: dict | None) -> dict | None:
    """Compute a lightweight game score from existing KPIs.

    This does not change the scientific model; it only turns the already-computed
    KPIs into an easy-to-play feedback loop.

    Score (0-100) = 70% efficacy + 30% safety
    - efficacy: tumor_reduction (higher is better)
    - safety: healthy_loss (lower is better; 0 points if >= 30% loss)
    """
    if not summary:
        return None

    try:
        tumor_reduction = float(summary.get("tumor_reduction") or 0.0)
    except (TypeError, ValueError):
        tumor_reduction = 0.0

    try:
        healthy_loss = float(summary.get("healthy_loss") or 0.0)
    except (TypeError, ValueError):
        healthy_loss = 0.0

    tumor_reduction = max(0.0, min(1.0, tumor_reduction))
    healthy_loss = max(0.0, min(1.0, healthy_loss))

    efficacy_points = 70.0 * tumor_reduction
    safety_points = 30.0 * max(0.0, min(1.0, 1.0 - (healthy_loss / 0.3)))
    score = int(round(efficacy_points + safety_points))

    win = bool(tumor_reduction >= 0.90 and healthy_loss <= 0.20)

    return {
        "score": score,
        "win": win,
        "tumor_reduction": tumor_reduction,
        "healthy_loss": healthy_loss,
        "mission": {
            "tumor_reduction_target": 0.90,
            "healthy_loss_max": 0.20,
        },
    }

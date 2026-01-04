from __future__ import annotations

import random

from .logic import activation_score
from .types import AntigenProfile, LogicGate


def binding_probability(
    profile: AntigenProfile,
    *,
    target_antigen: str,
    affinity: float,
    half_max_density: float = 50.0,
) -> float:
    """Probability-like binding proxy in [0,1].

    We model binding as a saturating function of antigen density.
    `affinity` scales the curve (higher = binds more strongly).
    """

    density = max(0.0, profile.get(target_antigen))
    k = max(1e-6, half_max_density) / max(1e-6, affinity)
    # Michaelis-Menten like saturation
    return float(density / (density + k))


def activation_probability(
    profile: AntigenProfile,
    *,
    logic: LogicGate | None,
    target_antigen: str | None = None,
    affinity: float = 1.0,
    rng: random.Random | None = None,
) -> float:
    """Activation probability for a modality that requires targeting.

    If `logic` is provided, it drives activation.
    Otherwise we fall back to single-target binding probability.
    """

    if logic is not None:
        return activation_score(profile, logic, rng=rng)

    if not target_antigen:
        return 0.0

    return binding_probability(profile, target_antigen=target_antigen, affinity=affinity)


def estimate_tp_fp(
    tumor_profiles: list[AntigenProfile],
    normal_profiles: list[AntigenProfile],
    *,
    logic: LogicGate | None,
    target_antigen: str | None = None,
    affinity: float = 1.0,
) -> dict[str, float]:
    """Estimate coverage/selectivity as TP/FP rates.

    TP = fraction of tumor profiles activated
    FP = fraction of normal profiles activated
    """

    rng = random.Random(0)

    def rate(profiles: list[AntigenProfile]) -> float:
        if not profiles:
            return 0.0
        hits = 0
        for p in profiles:
            if activation_probability(p, logic=logic, target_antigen=target_antigen, affinity=affinity, rng=rng) >= 0.5:
                hits += 1
        return hits / float(len(profiles))

    tp = rate(tumor_profiles)
    fp = rate(normal_profiles)
    return {"tp": float(tp), "fp": float(fp)}

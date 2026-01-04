from __future__ import annotations

import math
import random

from .types import AntigenProfile, LogicGate


def activation_score(profile: AntigenProfile, gate: LogicGate, *, rng: random.Random | None = None) -> float:
    """Compute an activation score in [0, 1] for an antigen logic gate.

    This is intentionally simple but *quantitative*:
    - Each antigen contributes via a sigmoid around its threshold.
    - AND uses product, OR uses noisy-OR.
    - NOT acts as an inhibitory multiplier.

    `gate.noise` adds small Gaussian noise on the final score.
    """

    rng = rng or random.Random()

    def sigmoid(x: float) -> float:
        # Stable-ish sigmoid
        if x >= 30:
            return 1.0
        if x <= -30:
            return 0.0
        return 1.0 / (1.0 + math.exp(-x))

    def antigen_term(antigen: str) -> float:
        density = profile.get(antigen)
        thr = float(gate.thresholds.get(antigen, 0.0))
        # scale makes threshold transition reasonably sharp
        scale = max(1e-6, 0.15 * max(thr, 1.0))
        z = (density - thr) / scale
        return float(sigmoid(z))

    pos = [antigen_term(a) for a in (gate.positive_antigens or [])]
    if not pos:
        pos_score = 0.0
    elif gate.kind == "AND":
        pos_score = 1.0
        for p in pos:
            pos_score *= p
    elif gate.kind == "OR":
        # noisy-OR: 1 - Π(1-p)
        miss = 1.0
        for p in pos:
            miss *= (1.0 - p)
        pos_score = 1.0 - miss
    elif gate.kind == "NOT":
        # NOT without positives is undefined; treat as 0 unless overridden by positives
        pos_score = 0.0
    else:
        pos_score = 0.0

    inhib = [antigen_term(a) for a in (gate.negative_antigens or [])]
    inhib_score = 1.0
    for p in inhib:
        # if inhibitory antigen is high -> strong inhibition
        inhib_score *= (1.0 - p)

    score = pos_score * inhib_score

    if gate.noise:
        score += rng.gauss(0.0, float(gate.noise))

    return float(min(1.0, max(0.0, score)))

from __future__ import annotations

import random

from .types import Clone, TumorState


def _normalize(fracs: dict[str, float]) -> dict[str, float]:
    total = sum(max(0.0, v) for v in fracs.values())
    if total <= 0:
        # fallback to uniform
        n = max(1, len(fracs))
        return {k: 1.0 / n for k in fracs}
    return {k: max(0.0, v) / total for k, v in fracs.items()}


def apply_selection(
    state: TumorState,
    *,
    clones: dict[str, Clone],
    survival_advantage: dict[str, float],
    drift: float = 0.01,
    rng: random.Random | None = None,
) -> TumorState:
    """Update clone fractions under selection + drift.

    `survival_advantage[clone]` > 1 means better survival (resistance).
    """

    rng = rng or random.Random()

    new_fracs: dict[str, dict[str, float]] = {}
    for comp, fracs in state.fractions.items():
        updated = {}
        for name, f in fracs.items():
            adv = float(survival_advantage.get(name, 1.0))
            noise = rng.uniform(-drift, drift)
            updated[name] = max(0.0, f * adv * (1.0 + noise))
        new_fracs[comp] = _normalize(updated)

    return TumorState(time_days=state.time_days, loads=dict(state.loads), fractions=new_fracs)


def antigen_escape(
    clones: dict[str, Clone],
    *,
    target_antigen: str,
    escape_rate: float,
    rng: random.Random | None = None,
) -> dict[str, Clone]:
    """Stochastic antigen-loss (escape) that reduces target antigen density in some clones."""

    rng = rng or random.Random()
    mutated: dict[str, Clone] = {}
    for name, clone in clones.items():
        if rng.random() < escape_rate:
            new_densities = dict(clone.antigen_profile.densities)
            new_densities[target_antigen] = 0.2 * float(new_densities.get(target_antigen, 0.0))
            mutated[name] = Clone(
                name=clone.name,
                antigen_profile=clone.antigen_profile.__class__(new_densities),
                sensitivities=dict(clone.sensitivities),
            )
        else:
            mutated[name] = clone
    return mutated

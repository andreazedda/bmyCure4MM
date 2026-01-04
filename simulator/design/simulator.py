from __future__ import annotations

import math
import random

from .evolution import antigen_escape, apply_selection
from .targeting import activation_probability
from .toxicity import estimate_toxicity
from .types import ADCModel, CARTModel, Clone, SmallMoleculeModel, StepOutcome, Therapy, TumorState


def _clamp01(x: float) -> float:
    return float(min(1.0, max(0.0, x)))


def _kill_strength(
    therapy: Therapy,
    *,
    clone: Clone,
    compartment: str,
    rng: random.Random,
) -> float:
    """Kill fraction per step in [0,1] with bulk/reservoir distinction."""

    base = 0.0

    if therapy.modality == "small_molecule":
        model = SmallMoleculeModel(**therapy.params)
        base = _clamp01(model.exposure * model.on_target_effect * clone.sensitivity("small_molecule"))

    elif therapy.modality == "adc":
        model = ADCModel(**therapy.params)
        act = activation_probability(
            clone.antigen_profile,
            logic=therapy.logic,
            target_antigen=model.target_antigen,
            affinity=model.affinity,
            rng=rng,
        )
        # Kill scales with payload potency + internalization
        base = _clamp01(model.payload_potency * model.internalization_rate * act * clone.sensitivity("adc"))

    elif therapy.modality == "car_t":
        model = CARTModel(**therapy.params)
        density = clone.antigen_profile.get(model.target_antigen)
        act = 1.0 if density >= model.activation_threshold else density / max(1e-6, model.activation_threshold)
        eff = _clamp01(model.expansion_rate * (1.0 - model.exhaustion_rate))
        base = _clamp01(eff * _clamp01(act) * clone.sensitivity("car_t"))

    # Reservoir is harder to clear (quiescence + antigen differences)
    if compartment == "reservoir":
        base *= 0.45

    return _clamp01(base)


def step(
    state: TumorState,
    *,
    clones: dict[str, Clone],
    therapy: Therapy,
    normal_compartments: list,
    dt_days: float,
    rng: random.Random | None = None,
) -> tuple[TumorState, StepOutcome, dict[str, Clone]]:
    rng = rng or random.Random()

    # Kill tumor loads with compartment-specific aggregate kill
    new_loads = dict(state.loads)

    for comp in ("bulk", "reservoir"):
        fracs = state.fractions.get(comp, {})
        comp_kill = 0.0
        for name, f in fracs.items():
            clone = clones.get(name)
            if not clone:
                continue
            k = _kill_strength(therapy, clone=clone, compartment=comp, rng=rng)
            comp_kill += f * k
        comp_kill = _clamp01(comp_kill)
        new_loads[comp] = max(0.0, new_loads.get(comp, 0.0) * (1.0 - comp_kill))

    # Evolution: selection towards more resistant clones
    survival_advantage = {}
    for name, clone in clones.items():
        # Resistance proxy: 1 + (1 - sensitivity)
        sens = clone.sensitivity(therapy.modality)
        survival_advantage[name] = 1.0 + 0.6 * (1.0 - sens)

    evolved_state = apply_selection(state, clones=clones, survival_advantage=survival_advantage, drift=0.02, rng=rng)
    evolved_state.loads = new_loads
    evolved_state.time_days = state.time_days + dt_days

    # Antigen escape under targeted pressures
    evolved_clones = clones
    if therapy.modality in ("adc", "car_t"):
        try:
            target = str(therapy.params.get("target_antigen") or therapy.target or "")
        except Exception:
            target = ""
        if target:
            evolved_clones = antigen_escape(evolved_clones, target_antigen=target, escape_rate=0.05, rng=rng)

    tox_scores, tox_why = estimate_toxicity(therapy, normal_compartments=normal_compartments)

    # Relapse risk proxy: reservoir remaining + diversity
    bulk = float(new_loads.get("bulk", 0.0))
    res = float(new_loads.get("reservoir", 0.0))
    diversity = 0.0
    fr = evolved_state.fractions.get("reservoir", {})
    for p in fr.values():
        if p > 1e-12:
            diversity -= p * math.log(p)
    relapse_risk = _clamp01(0.65 * (res / max(1e-9, bulk + res)) + 0.2 * _clamp01(diversity / 2.0) + 0.15 * tox_scores["marrow_reserve_depletion"])

    outcome = StepOutcome(
        time_days=float(evolved_state.time_days),
        tumor_load_bulk=float(bulk),
        tumor_load_reservoir=float(res),
        clone_fractions_bulk=dict(evolved_state.fractions.get("bulk", {})),
        clone_fractions_reservoir=dict(evolved_state.fractions.get("reservoir", {})),
        toxicity=tox_scores,
        toxicity_explanation=tox_why,
        relapse_risk=float(relapse_risk),
    )

    return evolved_state, outcome, evolved_clones

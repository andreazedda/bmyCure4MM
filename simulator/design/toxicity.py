from __future__ import annotations

from .targeting import activation_probability
from .types import (
    ADCModel,
    CARTModel,
    NormalCompartment,
    SmallMoleculeModel,
    Therapy,
    ToxicityCategory,
)


def _clamp01(x: float) -> float:
    return float(min(1.0, max(0.0, x)))


def estimate_toxicity(
    therapy: Therapy,
    *,
    normal_compartments: list[NormalCompartment],
) -> tuple[dict[ToxicityCategory, float], dict[ToxicityCategory, list[str]]]:
    """Mechanistic toxicity decomposition.

    Output is (scores, explanations). Each score is in [0,1].

    IMPORTANT: This intentionally distinguishes categories causally.
    It does not produce vague statements.
    """

    scores: dict[ToxicityCategory, float] = {
        "on_target_off_tumor": 0.0,
        "off_target_polypharmacology": 0.0,
        "payload_mediated": 0.0,
        "immune_mediated": 0.0,
        "marrow_reserve_depletion": 0.0,
    }
    why: dict[ToxicityCategory, list[str]] = {k: [] for k in scores}

    if therapy.modality == "small_molecule":
        model = SmallMoleculeModel(**therapy.params)
        scores["off_target_polypharmacology"] = _clamp01(model.off_target_polypharmacology)
        why["off_target_polypharmacology"].append(
            f"Polypharmacology baseline={model.off_target_polypharmacology:.2f} (diffusive systemic exposure)."
        )
        # marrow reserve depletion proxy: exposure harms dividing cells
        marrow = [c for c in normal_compartments if "marrow" in c.name.lower()]
        if marrow:
            depletion = _clamp01(0.4 * model.exposure)
            scores["marrow_reserve_depletion"] = max(scores["marrow_reserve_depletion"], depletion)
            why["marrow_reserve_depletion"].append(
                f"Marrow depletion scales with exposure: {model.exposure:.2f} → {depletion:.2f}."
            )

    elif therapy.modality == "adc":
        model = ADCModel(**therapy.params)
        # On-target/off-tumor: normal compartments with target antigen
        max_on_target = 0.0
        for comp in normal_compartments:
            bind = activation_probability(
                comp.antigen_profile,
                logic=therapy.logic,
                target_antigen=model.target_antigen,
                affinity=model.affinity,
            )
            comp_risk = comp.vulnerability * bind
            max_on_target = max(max_on_target, comp_risk)
            if comp_risk > 0.05:
                why["on_target_off_tumor"].append(
                    f"{comp.name} expresses {model.target_antigen}: binding≈{bind:.2f}, vulnerability={comp.vulnerability:.2f}."
                )
        scores["on_target_off_tumor"] = _clamp01(max_on_target)

        # Payload-mediated: bystander + systemic leakage + linker/payload potency
        payload = model.payload_potency * (0.5 * model.bystander_diffusion + 0.5 * model.systemic_leakage)
        payload *= (0.5 + 0.5 * model.linker_cleavage)
        scores["payload_mediated"] = _clamp01(payload)
        why["payload_mediated"].extend(
            [
                f"Bystander diffusion={model.bystander_diffusion:.2f} enables payload spillover to nearby non-target cells.",
                f"Systemic leakage={model.systemic_leakage:.2f} models FcR uptake/catabolism → circulating payload exposure.",
                f"Linker cleavage={model.linker_cleavage:.2f} controls payload release kinetics.",
                f"Payload potency={model.payload_potency:.2f} amplifies exposure-to-damage mapping.",
            ]
        )

        # Marrow reserve depletion: explicit chain for neutropenia
        marrow = [c for c in normal_compartments if "marrow" in c.name.lower() or "neut" in c.name.lower()]
        if marrow:
            marrow_bind = max(
                activation_probability(
                    c.antigen_profile,
                    logic=therapy.logic,
                    target_antigen=model.target_antigen,
                    affinity=model.affinity,
                )
                for c in marrow
            )
            depletion = _clamp01(0.6 * scores["payload_mediated"] + 0.4 * marrow_bind)
            scores["marrow_reserve_depletion"] = max(scores["marrow_reserve_depletion"], depletion)
            why["marrow_reserve_depletion"].extend(
                [
                    "Neutropenia mechanism chain:",
                    "ADC → (binding/internalization + systemic leakage/bystander) → payload exposure in marrow → damage to myeloid precursors → neutropenia.",
                    f"Computed marrow binding proxy={marrow_bind:.2f}; payload-mediated score={scores['payload_mediated']:.2f}.",
                ]
            )

    elif therapy.modality == "car_t":
        model = CARTModel(**therapy.params)
        # On-target/off-tumor depends on normal expression & activation threshold
        max_risk = 0.0
        for comp in normal_compartments:
            density = comp.antigen_profile.get(model.target_antigen)
            bind = 1.0 if density >= model.activation_threshold else density / max(1e-6, model.activation_threshold)
            comp_risk = comp.vulnerability * _clamp01(bind)
            max_risk = max(max_risk, comp_risk)
            if comp_risk > 0.05:
                why["on_target_off_tumor"].append(
                    f"{comp.name}: {model.target_antigen} density={density:.1f} vs threshold={model.activation_threshold:.1f} → activation≈{_clamp01(bind):.2f}."
                )
        scores["on_target_off_tumor"] = _clamp01(max_risk)

        # Immune-mediated toxicity (CRS/neurotox proxy): expansion * cytokine output
        immune = _clamp01(model.expansion_rate * model.cytokine_per_kill)
        scores["immune_mediated"] = immune
        why["immune_mediated"].extend(
            [
                "Immune-mediated toxicity is modeled as cytokine output from an expanding living agent.",
                f"Expansion rate={model.expansion_rate:.2f}; cytokine_per_kill={model.cytokine_per_kill:.2f} → immune score={immune:.2f}.",
                "This is distinct from chemical off-target toxicity.",
            ]
        )

    return scores, why

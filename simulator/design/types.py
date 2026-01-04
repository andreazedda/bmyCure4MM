from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal


CompartmentName = Literal["bulk", "reservoir"]
ToxicityCategory = Literal[
    "on_target_off_tumor",
    "off_target_polypharmacology",
    "payload_mediated",
    "immune_mediated",
    "marrow_reserve_depletion",
]
ModalityKind = Literal["small_molecule", "adc", "car_t"]
LogicKind = Literal["AND", "OR", "NOT"]


@dataclass(frozen=True)
class AntigenProfile:
    """Antigen densities for a cell (tumor clone or normal compartment).

    Values are abstracted densities (e.g., molecules/cell or normalized units).
    """

    densities: dict[str, float]

    def get(self, antigen: str, default: float = 0.0) -> float:
        try:
            return float(self.densities.get(antigen, default))
        except (TypeError, ValueError):
            return default


@dataclass(frozen=True)
class Clone:
    """A tumor subclone with an antigen profile and sensitivities.

    `sensitivities` is a modality-specific scalar in [0, 1] (higher = more sensitive).
    """

    name: str
    antigen_profile: AntigenProfile
    sensitivities: dict[ModalityKind, float] = field(default_factory=dict)

    def sensitivity(self, modality: ModalityKind) -> float:
        return float(self.sensitivities.get(modality, 0.5))


@dataclass
class TumorState:
    """State of the disease with two tumor compartments.

    `loads` are arbitrary cell counts or normalized burdens.
    `fractions` map clone name -> fraction within each compartment.
    """

    time_days: float
    loads: dict[CompartmentName, float]
    fractions: dict[CompartmentName, dict[str, float]]


@dataclass(frozen=True)
class NormalCompartment:
    """A normal compartment that can be harmed (e.g., marrow progenitors)."""

    name: str
    antigen_profile: AntigenProfile
    vulnerability: float = 1.0  # relative susceptibility to damage


@dataclass(frozen=True)
class LogicGate:
    """Boolean logic gate over antigen densities with thresholds and noise.

    - AND/OR operate over `positive_antigens`
    - NOT uses `negative_antigens` as inhibitory conditions

    Activation is computed as a continuous score in [0, 1].
    """

    kind: LogicKind
    positive_antigens: list[str]
    negative_antigens: list[str] = field(default_factory=list)
    thresholds: dict[str, float] = field(default_factory=dict)
    noise: float = 0.0


@dataclass(frozen=True)
class SmallMoleculeModel:
    """Abstract but mechanistic small molecule PK/PD envelope."""

    name: str
    exposure: float  # normalized AUC-like
    on_target_effect: float  # kill strength on tumor given sensitivity
    off_target_polypharmacology: float  # baseline off-target toxicity


@dataclass(frozen=True)
class ADCModel:
    """ADC delivery model (binding/internalization/linker/payload/bystander/leakage)."""

    name: str
    target_antigen: str
    affinity: float
    internalization_rate: float
    linker_cleavage: float
    payload_potency: float
    bystander_diffusion: float
    systemic_leakage: float  # FcR uptake / catabolism -> payload elsewhere


@dataclass(frozen=True)
class CARTModel:
    """CAR-T dynamics proxy."""

    name: str
    target_antigen: str
    activation_threshold: float
    expansion_rate: float
    exhaustion_rate: float
    cytokine_per_kill: float


@dataclass(frozen=True)
class Therapy:
    """A therapy strategy: modality + target logic (for logic-gated constructs)."""

    name: str
    modality: ModalityKind
    target: str | None = None  # e.g. BCMA (target), not the modality
    logic: LogicGate | None = None
    params: dict[str, Any] = field(default_factory=dict)


@dataclass
class StepOutcome:
    """One simulation step outcome including artifacts and explanations."""

    time_days: float
    tumor_load_bulk: float
    tumor_load_reservoir: float
    clone_fractions_bulk: dict[str, float]
    clone_fractions_reservoir: dict[str, float]
    toxicity: dict[ToxicityCategory, float]
    toxicity_explanation: dict[ToxicityCategory, list[str]]
    relapse_risk: float


@dataclass
class DesignReport:
    """JSON-ready report with required artifacts."""

    inputs: dict[str, Any]
    tradeoff_points: list[dict[str, Any]]
    pareto_front: list[dict[str, Any]]
    dynamics: list[dict[str, Any]]
    relapse_forecast: dict[str, Any]
    toxicity_decomposition: dict[str, Any]
    rationale: dict[str, Any]

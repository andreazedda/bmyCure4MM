from __future__ import annotations

from enum import Enum
from typing import Any

INTENDED_USE_LEVEL = "E1_research_prototype"
APP_VERSION = "0.1.0"
API_VERSION = "research-api-v1"
CLINICAL_DECISION_SUPPORT = False
PATIENT_SPECIFIC_PREDICTION_VALIDATED = False
CAUSAL_EFFECT_IDENTIFIED = False
CURRENT_RESEARCH_MODEL_VERSION = "research-twin-v1"
CLAIMS_POLICY_VERSION = "claims-policy-v1"
MODEL_CARD_VERSION = "model-card-v1"
EVIDENCE_GRAPH_VERSION = "evidence-graph-v1"
VALIDATION_PROTOCOL_VERSION = "validation-protocol-v1"
REPORT_TEMPLATE_VERSION = "research-report-v1"


class EpistemicLabel(str, Enum):
    UNKNOWN = "UNKNOWN"
    NEEDS_EVIDENCE = "NEEDS_EVIDENCE"
    OBSERVED = "OBSERVED"
    DERIVED = "DERIVED"
    USER_PROVIDED = "USER_PROVIDED"
    SIMULATED = "SIMULATED"
    HEURISTIC = "HEURISTIC"
    LITERATURE_INFORMED = "LITERATURE_INFORMED"
    HYPOTHETICAL = "HYPOTHETICAL"
    VALIDATED_EXTERNAL = "VALIDATED_EXTERNAL"


FORBIDDEN_CLAIM_CODES = (
    "recommended_or_best_treatment",
    "patient_specific_dose_instruction",
    "predicted_patient_benefit",
    "clinically_superior_regimen",
    "identified_causal_treatment_effect",
    "validated_patient_prognosis_without_governed_evidence",
)


def governance_metadata(
    *,
    epistemic_label: EpistemicLabel | str,
    model_version: str = CURRENT_RESEARCH_MODEL_VERSION,
    output_kind: str,
) -> dict[str, Any]:
    label = epistemic_label.value if isinstance(epistemic_label, EpistemicLabel) else str(epistemic_label)
    return {
        "claims_policy_version": CLAIMS_POLICY_VERSION,
        "intended_use_level": INTENDED_USE_LEVEL,
        "clinical_decision_support": CLINICAL_DECISION_SUPPORT,
        "patient_specific_prediction_validated": (PATIENT_SPECIFIC_PREDICTION_VALIDATED),
        "causal_effect_identified": CAUSAL_EFFECT_IDENTIFIED,
        "epistemic_label": label,
        "model_version": model_version,
        "output_kind": output_kind,
    }


def template_governance_context() -> dict[str, Any]:
    return {
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.HYPOTHETICAL,
            output_kind="global_interface_context",
        ),
        "epistemic_labels": [label.value for label in EpistemicLabel],
    }


def _as_float(payload: dict[str, object], name: str) -> float | None:
    value = payload.get(name)
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def build_model_relative_diagnostics(
    summary: dict[str, object] | None,
) -> dict[str, Any] | None:
    """Classify simulated outputs without proposing a clinical action.

    Thresholds are navigation heuristics over model outputs. They do not
    estimate patient benefit, establish tolerability, or identify a causal
    treatment effect.
    """
    if not summary:
        return None

    tumor_reduction = _as_float(summary, "tumor_reduction")
    healthy_loss = _as_float(summary, "healthy_loss")
    time_to_recurrence = _as_float(summary, "time_to_recurrence")
    signals: list[dict[str, str]] = []

    def add_signal(
        code: str,
        title_en: str,
        title_it: str,
        detail_en: str,
        detail_it: str,
    ) -> None:
        signals.append(
            {
                "code": code,
                "title_en": title_en,
                "title_it": title_it,
                "detail_en": detail_en,
                "detail_it": detail_it,
                "epistemic_label": EpistemicLabel.HEURISTIC.value,
            }
        )

    if healthy_loss is None:
        add_signal(
            "insufficient_evidence",
            "Healthy-cell constraint unavailable",
            "Vincolo sulle cellule sane non disponibile",
            "The run does not contain the modeled healthy-cell-loss metric.",
            "La simulazione non contiene la metrica di perdita delle cellule sane.",
        )
    elif healthy_loss >= 0.30:
        add_signal(
            "constraint_flag",
            "Modeled healthy-cell-loss constraint flag",
            "Segnale di vincolo sulla perdita di cellule sane",
            "The simulated healthy-cell-loss metric is at or above the 0.30 heuristic threshold.",
            "La metrica simulata di perdita delle cellule sane è pari o superiore alla soglia euristica 0,30.",
        )
    elif healthy_loss >= 0.20:
        add_signal(
            "constraint_watch_zone",
            "Modeled healthy-cell-loss watch zone",
            "Zona di attenzione per la perdita di cellule sane",
            "The simulated metric is between the 0.20 and 0.30 heuristic thresholds.",
            "La metrica simulata è compresa tra le soglie euristiche 0,20 e 0,30.",
        )

    if tumor_reduction is None:
        add_signal(
            "insufficient_evidence",
            "Tumor-response signal unavailable",
            "Segnale di risposta tumorale non disponibile",
            "The run does not contain the modeled tumor-reduction metric.",
            "La simulazione non contiene la metrica di riduzione tumorale.",
        )
    elif tumor_reduction < 0:
        add_signal(
            "model_regrowth_signal",
            "Model regrowth signal",
            "Segnale di ricrescita del modello",
            "The simulated endpoint exceeds its modeled starting tumor burden.",
            "L'endpoint simulato supera il carico tumorale iniziale del modello.",
        )
    elif tumor_reduction < 0.30:
        add_signal(
            "simulated_low_impact_zone",
            "Simulated low-impact zone",
            "Zona simulata a basso impatto",
            "The model-relative reduction is below the 0.30 heuristic threshold.",
            "La riduzione relativa al modello è inferiore alla soglia euristica 0,30.",
        )
    elif tumor_reduction >= 0.50:
        add_signal(
            "simulated_high_impact_zone",
            "Simulated high-impact zone",
            "Zona simulata ad alto impatto",
            "The model-relative reduction is at or above the 0.50 heuristic threshold.",
            "La riduzione relativa al modello è pari o superiore alla soglia euristica 0,50.",
        )

    if time_to_recurrence is None:
        add_signal(
            "recurrence_not_observed_in_horizon",
            "No modeled recurrence time in the reported horizon",
            "Nessun tempo di recidiva modellato nell'orizzonte riportato",
            "This is a simulated horizon result, not a patient-specific prognosis.",
            "È un risultato dell'orizzonte simulato, non una prognosi specifica del paziente.",
        )
    elif time_to_recurrence < 90:
        add_signal(
            "short_model_durability_signal",
            "Short model-durability signal",
            "Segnale di breve durabilità del modello",
            "The simulated recurrence time is below the 90-day heuristic threshold.",
            "Il tempo di recidiva simulato è inferiore alla soglia euristica di 90 giorni.",
        )

    if not signals:
        add_signal(
            "no_threshold_flag",
            "No configured diagnostic threshold flag",
            "Nessun segnale dalle soglie diagnostiche configurate",
            "No current heuristic threshold was activated by the reported model outputs.",
            "Nessuna soglia euristica attuale è stata attivata dagli output del modello.",
        )

    return {
        "status": (
            "constraint_flag"
            if any(signal["code"] == "constraint_flag" for signal in signals)
            else "model_relative_diagnostics"
        ),
        "signals": signals,
        "has_signals": bool(signals),
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.SIMULATED,
            output_kind="model_relative_diagnostics",
        ),
        "limitations": (
            "Threshold flags describe this mechanistic simulation only; "
            "they do not authorize a treatment or dose change."
        ),
    }

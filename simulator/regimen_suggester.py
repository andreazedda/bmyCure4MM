"""Heuristic regimen-catalog grouping for exploratory research interfaces.

This module preserves legacy transport keys such as ``preferred`` so existing
callers remain compatible. Those keys describe display buckets only. They are
not treatment recommendations, estimates of patient benefit, or safety
determinations. The catalog's historical evidence metadata has not yet passed
source-level verification, so it fails closed as ``NEEDS_EVIDENCE``.
"""

from __future__ import annotations

from dataclasses import dataclass

from mmportal.governance import EpistemicLabel, governance_metadata


@dataclass(frozen=True)
class RegimenSuggestion:
    """Minimal non-prescriptive catalog record.

    Some field names are retained for response-schema compatibility. In this
    E1 release, ``strength`` is an epistemic placeholder rather than a clinical
    recommendation grade, and unverified fields fail closed.
    """

    name: str
    full_name: str
    components: list[str]
    indication: str
    strength: str = "not_assessed"
    evidence_level: str = "NEEDS_EVIDENCE"
    key_trials: tuple[str, ...] = ()
    considerations: tuple[str, ...] = ()
    contraindications: tuple[str, ...] = ()
    expected_response_rate: str = "NEEDS_EVIDENCE"

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "full_name": self.full_name,
            "components": list(self.components),
            "indication": self.indication,
            "strength": self.strength,
            "evidence_level": self.evidence_level,
            "key_trials": list(self.key_trials),
            "considerations": list(self.considerations),
            "contraindications": list(self.contraindications),
            "expected_response_rate": self.expected_response_rate,
            "epistemic_label": "NEEDS_EVIDENCE",
        }


def _entry(name: str, full_name: str, components: list[str], context: str) -> RegimenSuggestion:
    return RegimenSuggestion(
        name=name,
        full_name=full_name,
        components=components,
        indication=f"Historical catalog context: {context}; source verification pending",
    )


FRONTLINE_REGIMENS = {
    "VRd": _entry(
        "VRd",
        "Bortezomib, Lenalidomide, Dexamethasone",
        ["Bortezomib", "Lenalidomide", "Dexamethasone"],
        "newly diagnosed",
    ),
    "Dara-VRd": _entry(
        "Dara-VRd",
        "Daratumumab + Bortezomib, Lenalidomide, Dexamethasone",
        ["Daratumumab", "Bortezomib", "Lenalidomide", "Dexamethasone"],
        "newly diagnosed",
    ),
    "KRd": _entry(
        "KRd",
        "Carfilzomib, Lenalidomide, Dexamethasone",
        ["Carfilzomib", "Lenalidomide", "Dexamethasone"],
        "newly diagnosed",
    ),
    "DRd": _entry(
        "DRd",
        "Daratumumab, Lenalidomide, Dexamethasone",
        ["Daratumumab", "Lenalidomide", "Dexamethasone"],
        "newly diagnosed",
    ),
    "VMP": _entry(
        "VMP",
        "Bortezomib, Melphalan, Prednisone",
        ["Bortezomib", "Melphalan", "Prednisone"],
        "newly diagnosed",
    ),
    "Dara-VMP": _entry(
        "Dara-VMP",
        "Daratumumab + Bortezomib, Melphalan, Prednisone",
        ["Daratumumab", "Bortezomib", "Melphalan", "Prednisone"],
        "newly diagnosed",
    ),
    "Rd": _entry(
        "Rd",
        "Lenalidomide, Dexamethasone",
        ["Lenalidomide", "Dexamethasone"],
        "newly diagnosed",
    ),
}

RELAPSED_REGIMENS = {
    "DPd": _entry(
        "DPd",
        "Daratumumab, Pomalidomide, Dexamethasone",
        ["Daratumumab", "Pomalidomide", "Dexamethasone"],
        "previously treated",
    ),
    "DVd": _entry(
        "DVd",
        "Daratumumab, Bortezomib, Dexamethasone",
        ["Daratumumab", "Bortezomib", "Dexamethasone"],
        "previously treated",
    ),
    "KPd": _entry(
        "KPd",
        "Carfilzomib, Pomalidomide, Dexamethasone",
        ["Carfilzomib", "Pomalidomide", "Dexamethasone"],
        "previously treated",
    ),
    "IsaPd": _entry(
        "Isa-Pd",
        "Isatuximab, Pomalidomide, Dexamethasone",
        ["Isatuximab", "Pomalidomide", "Dexamethasone"],
        "previously treated",
    ),
    "Selinexor-Vd": _entry(
        "SVd",
        "Selinexor, Bortezomib, Dexamethasone",
        ["Selinexor", "Bortezomib", "Dexamethasone"],
        "previously treated",
    ),
    "Elotuzumab-Pd": _entry(
        "Elo-Pd",
        "Elotuzumab, Pomalidomide, Dexamethasone",
        ["Elotuzumab", "Pomalidomide", "Dexamethasone"],
        "previously treated",
    ),
}

NOVEL_THERAPIES = {
    "Ide-cel": _entry(
        "Ide-cel",
        "Idecabtagene vicleucel (CAR-T)",
        ["Anti-BCMA CAR-T cells"],
        "advanced-line catalog",
    ),
    "Cilta-cel": _entry(
        "Cilta-cel",
        "Ciltacabtagene autoleucel (CAR-T)",
        ["Anti-BCMA CAR-T cells"],
        "advanced-line catalog",
    ),
    "Teclistamab": _entry(
        "Teclistamab",
        "Teclistamab (BCMA bispecific)",
        ["Teclistamab"],
        "advanced-line catalog",
    ),
    "Elranatamab": _entry(
        "Elranatamab",
        "Elranatamab (BCMA bispecific)",
        ["Elranatamab"],
        "advanced-line catalog",
    ),
    "Talquetamab": _entry(
        "Talquetamab",
        "Talquetamab (GPRC5D bispecific)",
        ["Talquetamab"],
        "advanced-line catalog",
    ),
}


def _catalog_card(entry: RegimenSuggestion, rule_code: str) -> dict[str, object]:
    return {
        **entry.to_dict(),
        "rationale": (
            f"Activated heuristic catalog rule {rule_code}; this rule does not "
            "estimate patient benefit or determine a clinical action."
        ),
    }


def suggest_regimens(
    age: int | None = None,
    transplant_eligible: bool | None = None,
    ecog: int | None = None,
    r_iss: str | None = None,
    high_risk_cytogenetics: bool = False,
    line_of_therapy: int = 1,
    prior_therapies: list[str] | None = None,
    refractory_to: list[str] | None = None,
    comorbidities: list[str] | None = None,
    renal_function: str | None = None,
    neuropathy_grade: int = 0,
    cardiac_issues: bool = False,
) -> dict[str, object]:
    """Group catalog entries with transparent heuristic rules.

    The returned ``preferred`` and ``alternative`` keys are legacy display
    bucket names. They carry no clinical preference, benefit, or safety claim.
    """
    del age, r_iss, comorbidities
    prior_therapies = prior_therapies or []
    refractory_lower = [item.lower() for item in (refractory_to or [])]

    catalog_governance = governance_metadata(
        epistemic_label=EpistemicLabel.HEURISTIC,
        output_kind="regimen_catalog_grouping",
    )
    catalog_governance["evidence_status"] = EpistemicLabel.NEEDS_EVIDENCE.value

    result: dict[str, object] = {
        "preferred": [],
        "alternative": [],
        "consider_in_clinical_trial": [],
        "avoid": [],
        "patient_factors": [],
        "disclaimer": (
            "Exploratory heuristic catalog grouping only. No treatment, dose, "
            "patient benefit, clinical safety, or causal effect is determined."
        ),
        "governance": catalog_governance,
    }
    higher = result["preferred"]
    alternative = result["alternative"]
    later_line = result["consider_in_clinical_trial"]
    constraints = result["avoid"]
    factors = result["patient_factors"]
    assert isinstance(higher, list)
    assert isinstance(alternative, list)
    assert isinstance(later_line, list)
    assert isinstance(constraints, list)
    assert isinstance(factors, list)

    if transplant_eligible is None:
        factors.append("Transplant-status input is unresolved.")
    if high_risk_cytogenetics:
        factors.append("A structured high-risk-cytogenetics flag activated a catalog branch.")
    if neuropathy_grade >= 2:
        factors.append("A neuropathy input activated a model constraint branch.")
    if cardiac_issues:
        factors.append("A cardiac-comorbidity input activated a model constraint branch.")
    if renal_function == "dialysis":
        factors.append("A dialysis input activated a renal model constraint branch.")

    if line_of_therapy == 1:
        if transplant_eligible:
            key = "Dara-VRd" if high_risk_cytogenetics else "VRd"
            higher.append(_catalog_card(FRONTLINE_REGIMENS[key], "F01"))
            alternative.append(_catalog_card(FRONTLINE_REGIMENS["KRd"], "F02"))
            if not high_risk_cytogenetics:
                alternative.append(_catalog_card(FRONTLINE_REGIMENS["Dara-VRd"], "F03"))
        else:
            higher.append(_catalog_card(FRONTLINE_REGIMENS["DRd"], "F04"))
            if high_risk_cytogenetics:
                alternative.append(_catalog_card(FRONTLINE_REGIMENS["Dara-VMP"], "F05"))
            if ecog is not None and ecog >= 2:
                alternative.append(_catalog_card(FRONTLINE_REGIMENS["Rd"], "F06"))
            if renal_function == "dialysis":
                alternative.append(_catalog_card(FRONTLINE_REGIMENS["VMP"], "F07"))
    elif line_of_therapy >= 2:
        len_refractory = any("len" in value for value in refractory_lower)
        bort_refractory = any("bort" in value for value in refractory_lower)
        dara_refractory = any("dara" in value for value in refractory_lower)
        if line_of_therapy == 2:
            if len_refractory and not dara_refractory:
                higher.extend(
                    [
                        _catalog_card(RELAPSED_REGIMENS["DPd"], "R01"),
                        _catalog_card(RELAPSED_REGIMENS["DVd"], "R02"),
                    ]
                )
            elif not dara_refractory:
                higher.extend(
                    [
                        _catalog_card(RELAPSED_REGIMENS["DVd"], "R03"),
                        _catalog_card(RELAPSED_REGIMENS["DPd"], "R04"),
                    ]
                )
            if bort_refractory and neuropathy_grade < 2:
                alternative.append(_catalog_card(RELAPSED_REGIMENS["KPd"], "R05"))
            alternative.append(_catalog_card(RELAPSED_REGIMENS["IsaPd"], "R06"))
        else:
            if not dara_refractory:
                higher.append(_catalog_card(RELAPSED_REGIMENS["DPd"], "R07"))
            alternative.extend(
                [
                    _catalog_card(RELAPSED_REGIMENS["Selinexor-Vd"], "R08"),
                    _catalog_card(RELAPSED_REGIMENS["Elotuzumab-Pd"], "R09"),
                ]
            )
        if line_of_therapy >= 4:
            later_line.extend(
                [
                    _catalog_card(NOVEL_THERAPIES["Cilta-cel"], "L01"),
                    _catalog_card(NOVEL_THERAPIES["Teclistamab"], "L02"),
                ]
            )
            bcma_exposed = any(
                "bcma" in item.lower() or "car-t" in item.lower()
                for item in prior_therapies
            )
            if bcma_exposed:
                later_line.append(_catalog_card(NOVEL_THERAPIES["Talquetamab"], "L03"))

    if neuropathy_grade >= 3:
        constraints.append({"agent": "Bortezomib", "reason": "Constraint rule C01 activated."})
    if cardiac_issues:
        constraints.append({"agent": "Carfilzomib", "reason": "Constraint rule C02 activated."})
    if renal_function == "dialysis":
        constraints.append({"agent": "Lenalidomide", "reason": "Constraint rule C03 activated."})

    return result


def get_regimen_details(regimen_name: str) -> dict[str, object] | None:
    """Return a catalog record without adding a clinical interpretation."""
    all_regimens = {**FRONTLINE_REGIMENS, **RELAPSED_REGIMENS, **NOVEL_THERAPIES}
    if regimen_name in all_regimens:
        return all_regimens[regimen_name].to_dict()
    for name, regimen in all_regimens.items():
        if name.lower() == regimen_name.lower():
            return regimen.to_dict()
    return None


def get_all_regimens() -> dict[str, list[dict[str, object]]]:
    """Return the fail-closed catalog grouped by historical setting."""
    return {
        "frontline": [item.to_dict() for item in FRONTLINE_REGIMENS.values()],
        "relapsed": [item.to_dict() for item in RELAPSED_REGIMENS.values()],
        "novel": [item.to_dict() for item in NOVEL_THERAPIES.values()],
    }

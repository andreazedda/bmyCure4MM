from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass
class MathematicalModel:
    """Coupled PK/PD and logistic growth model for MM simulation."""

    baseline_tumor_cells: float
    baseline_healthy_cells: float
    drug_doses: Dict[str, float]
    pk_params: Dict[str, Dict[str, float]]
    pd_params: Dict[str, Dict[str, float]]
    growth_rates: Dict[str, float]
    interaction_matrix: np.ndarray
    time_span: Tuple[float, float]
    carrying_capacity_factor_tumor: float = 10.0
    carrying_capacity_factor_healthy: float = 1.2
    evaluation_points: int = 200

    def simulate(self) -> pd.DataFrame:
        """Integrate the coupled system and return trajectories."""
        start, end = self.time_span
        if end <= start:
            raise ValueError("time_span end must be greater than start.")

        drug_names = list(self.drug_doses.keys())
        n_drugs = len(drug_names)
        if n_drugs == 0:
            raise ValueError("At least one drug dose must be provided.")

        if self.interaction_matrix.size == 0:
            interaction = np.zeros((n_drugs, n_drugs))
        else:
            interaction = np.array(self.interaction_matrix, dtype=float)
            if interaction.shape != (n_drugs, n_drugs):
                interaction = np.zeros((n_drugs, n_drugs))

        time_horizon = end - start
        dose_rates = {
            drug: max(dose, 0.0) / max(time_horizon, 1e-6) for drug, dose in self.drug_doses.items()
        }

        def initial_concentrations() -> Iterable[float]:
            concentrations: list[float] = []
            for drug in drug_names:
                pk = self.pk_params.get(drug, {})
                vd = pk.get("Vd", 1.0)
                concentrations.append(self.drug_doses.get(drug, 0.0) / max(vd, 1e-6))
            return concentrations

        y0 = np.concatenate(
            (
                np.array(
                    [
                        self.baseline_tumor_cells,
                        self.baseline_healthy_cells,
                    ],
                    dtype=float,
                ),
                np.array(list(initial_concentrations()), dtype=float),
            )
        )

        carrying_t = self.baseline_tumor_cells * self.carrying_capacity_factor_tumor
        carrying_h = self.baseline_healthy_cells * self.carrying_capacity_factor_healthy

        def pkpd_effects(concentrations: np.ndarray) -> np.ndarray:
            effects = np.zeros_like(concentrations, dtype=float)
            for idx, drug in enumerate(drug_names):
                pd = self.pd_params.get(drug, {})
                emax = pd.get("Emax", 0.0)
                ec50 = pd.get("EC50", 1.0)
                conc = max(concentrations[idx], 0.0)
                effects[idx] = emax * conc / (ec50 + conc + 1e-9)
            interaction_term = interaction @ effects
            return np.clip(effects + interaction_term, 0.0, 1.0)

        def rhs(_t: float, y: np.ndarray) -> np.ndarray:
            tumor = max(y[0], 0.0)
            healthy = max(y[1], 0.0)
            concentrations = np.maximum(y[2:], 0.0)
            effects = pkpd_effects(concentrations)
            total_effect = effects.sum()
            toxicity_effect = effects.mean() if effects.size else 0.0
            d_tumor = (
                self.growth_rates.get("tumor", 0.0)
                * tumor
                * (1.0 - tumor / max(carrying_t, 1e-6))
                - total_effect * tumor
            )
            d_healthy = (
                self.growth_rates.get("healthy", 0.0)
                * healthy
                * (1.0 - healthy / max(carrying_h, 1e-6))
                - toxicity_effect * healthy
            )
            d_concentrations = np.zeros_like(concentrations)
            for idx, drug in enumerate(drug_names):
                pk = self.pk_params.get(drug, {})
                half_life_hours = pk.get("half_life", 24.0)
                half_life_days = half_life_hours / 24.0
                k_elim = pk.get("k_elim", np.log(2) / max(half_life_days, 1e-6))
                input_rate = dose_rates[drug]
                d_concentrations[idx] = -k_elim * concentrations[idx] + input_rate
            return np.concatenate(([d_tumor, d_healthy], d_concentrations))

        t_eval = np.linspace(start, end, self.evaluation_points)
        solution = solve_ivp(
            rhs,
            t_span=self.time_span,
            y0=y0,
            t_eval=t_eval,
            vectorized=False,
            rtol=1e-6,
            atol=1e-8,
        )
        if not solution.success:
            raise RuntimeError(f"ODE solver failed: {solution.message}")

        data = {
            "time": solution.t,
            "tumor_cells": solution.y[0],
            "healthy_cells": solution.y[1],
        }
        for idx, drug in enumerate(drug_names):
            data[f"{drug}_concentration"] = solution.y[2 + idx]
        return pd.DataFrame(data)

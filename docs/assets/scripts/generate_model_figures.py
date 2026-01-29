#!/usr/bin/env python3
"""
Generate static SVG figures used by the documentation.

The docs are static (GitHub Pages), so we commit the generated SVGs under:
  docs/assets/images/models/
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "assets" / "images" / "models"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _save(fig: plt.Figure, name: str) -> None:
    path = OUT_DIR / name
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def figure_logistic_growth() -> None:
    t = np.linspace(0, 180, 400)  # days
    k = 1e12
    t0 = 1e9
    rates = [0.01, 0.02, 0.04]

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    for r in rates:
        # Logistic closed form
        tt = k / (1.0 + ((k - t0) / t0) * np.exp(-r * t))
        ax.plot(t, tt, label=f"r={r:.2f}/day")
    ax.set_title("Logistic growth: T(t) with different growth rates")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Tumor cells T(t)")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right")
    _save(fig, "logistic_growth.svg")


def figure_emax_curve() -> None:
    c = np.linspace(0, 200, 400)  # mg/L (illustrative)
    emax = 0.9
    ec50s = [5, 25, 75]

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    for ec50 in ec50s:
        e = emax * c / (ec50 + c + 1e-9)
        ax.plot(c, e, label=f"EC50={ec50:g}")
    ax.set_title("Emax PD curve: effect vs concentration")
    ax.set_xlabel("Concentration C (a.u.)")
    ax.set_ylabel("Effect E(C) (0..1)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    _save(fig, "emax_curve.svg")


def figure_pk_profiles() -> None:
    # Simple 1-compartment with first-order elimination and two input schedules.
    t = np.linspace(0, 56, 800)  # days
    half_life_days = 2.0
    k = np.log(2) / half_life_days

    def simulate(u_fn):
        c = np.zeros_like(t)
        dt = float(t[1] - t[0])
        for i in range(1, len(t)):
            c[i] = c[i - 1] + dt * (-k * c[i - 1] + u_fn(t[i - 1]))
        return c

    dose = 100.0
    horizon = 56.0
    continuous_u = lambda _tt: dose / horizon

    days = {0, 7, 14, 21, 28, 35, 42, 49}
    window = 0.35  # days

    def weekly_u(tt: float) -> float:
        d = int(np.floor(tt))
        frac = tt - d
        if d in days and frac < window:
            return dose / window
        return 0.0

    c_cont = simulate(continuous_u)
    c_week = simulate(weekly_u)

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.plot(t, c_cont, label="Continuous infusion-like (u(t)=Dose/H)")
    ax.plot(t, c_week, label="Weekly pulses (windowed)")
    ax.set_title("PK toy model: effect of dosing schedule on concentration")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Concentration C(t) (a.u.)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    _save(fig, "pk_profiles.svg")


def figure_coupled_dynamics() -> None:
    # Illustrative coupled tumor/healthy dynamics with a fixed drug effect.
    t = np.linspace(0, 180, 800)
    dt = float(t[1] - t[0])

    r_t = 0.02
    r_h = 0.015
    k_t = 1e12
    k_h = 5e11
    immune = 1.0

    def simulate(total_effect: float, tox_effect: float):
        T = np.zeros_like(t)
        H = np.zeros_like(t)
        T[0] = 1e9
        H[0] = 5e11
        for i in range(1, len(t)):
            dT = r_t * T[i - 1] * (1.0 - T[i - 1] / k_t) - total_effect * T[i - 1]
            dH = r_h * H[i - 1] * (1.0 - H[i - 1] / k_h) - immune * tox_effect * H[i - 1]
            T[i] = max(0.0, T[i - 1] + dt * dT)
            H[i] = max(0.0, H[i - 1] + dt * dH)
        return T, H

    T0, H0 = simulate(total_effect=0.0, tox_effect=0.0)
    T1, H1 = simulate(total_effect=0.03, tox_effect=0.01)

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.plot(t, T0, label="Tumor (no drug)", color="#d62728", alpha=0.6)
    ax.plot(t, T1, label="Tumor (with drug)", color="#d62728")
    ax.plot(t, H0, label="Healthy (no drug)", color="#1f77b4", alpha=0.6)
    ax.plot(t, H1, label="Healthy (with drug)", color="#1f77b4")
    ax.set_title("Coupled dynamics: tumor kill vs healthy toxicity (illustrative)")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Cells (log scale)")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left")
    _save(fig, "coupled_dynamics.svg")

def figure_uncertainty_bands() -> None:
    # Illustrative uncertainty bands for tumor reduction over time (toy).
    rng = np.random.default_rng(7)
    t = np.linspace(0, 180, 181)
    baseline = np.exp(-t / 80.0)  # decreasing proxy

    samples = []
    for _ in range(200):
        noise = rng.normal(0.0, 0.06, size=t.shape)
        trend = baseline * np.exp(noise)
        samples.append(trend)
    arr = np.asarray(samples)
    p05 = np.percentile(arr, 5, axis=0)
    p50 = np.percentile(arr, 50, axis=0)
    p95 = np.percentile(arr, 95, axis=0)

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.fill_between(t, p05, p95, alpha=0.25, label="p05–p95 band")
    ax.plot(t, p50, label="median (p50)")
    ax.set_title("Uncertainty bands (toy): trajectory variability")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Relative tumor burden (a.u.)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    _save(fig, "uncertainty_bands.svg")


def figure_pareto_front() -> None:
    # Illustrative Pareto front scatter in 2D (efficacy vs toxicity).
    rng = np.random.default_rng(13)
    n = 120
    efficacy = rng.uniform(0.2, 0.98, size=n)
    toxicity = rng.uniform(0.05, 0.45, size=n)
    # Create a few "good" points
    efficacy[:15] = rng.uniform(0.75, 0.98, size=15)
    toxicity[:15] = rng.uniform(0.06, 0.22, size=15)

    # Pareto: maximize efficacy, minimize toxicity
    points = np.column_stack([efficacy, toxicity])
    pareto = np.ones(n, dtype=bool)
    for i in range(n):
        if not pareto[i]:
            continue
        for j in range(n):
            if i == j:
                continue
            if (points[j, 0] >= points[i, 0]) and (points[j, 1] <= points[i, 1]) and (
                (points[j, 0] > points[i, 0]) or (points[j, 1] < points[i, 1])
            ):
                pareto[i] = False
                break

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.scatter(efficacy, toxicity, s=18, alpha=0.6, label="Trials")
    ax.scatter(efficacy[pareto], toxicity[pareto], s=28, label="Pareto front")
    ax.set_title("Pareto front (toy): efficacy vs toxicity")
    ax.set_xlabel("Efficacy (tumor_reduction)")
    ax.set_ylabel("Toxicity proxy (healthy_loss)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    _save(fig, "pareto_front.svg")


def main() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "font.size": 10,
        }
    )
    figure_logistic_growth()
    figure_emax_curve()
    figure_pk_profiles()
    figure_coupled_dynamics()
    figure_uncertainty_bands()
    figure_pareto_front()
    print(f"Generated figures in {OUT_DIR}")


if __name__ == "__main__":
    main()

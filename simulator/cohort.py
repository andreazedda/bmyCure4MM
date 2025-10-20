"""
Virtual cohort simulation for in-silico clinical trials.

Generates synthetic patient populations with parameter variability
and runs regimen simulations across the cohort to assess aggregate outcomes.
"""
from __future__ import annotations

import random
from typing import Dict, List

import numpy as np
from django.db import transaction
from django.utils import timezone

from . import models


def sample_patient_params(seed: int, n: int) -> List[Dict]:
    """
    Generate n synthetic patients with parameter variability.
    
    Uses log-normal distributions for biological parameters to ensure
    positive values with realistic spread.
    
    Parameters:
        seed: Random seed for reproducibility
        n: Number of patients to generate
        
    Returns:
        List of parameter dicts with baseline cells, growth rates, PK/PD noise
    """
    rng = np.random.default_rng(seed)
    
    patients = []
    for i in range(n):
        # Baseline tumor burden: log-normal around 1e9 cells (±50% CV)
        tumor_baseline = rng.lognormal(
            mean=np.log(1.0e9),
            sigma=0.5,  # CV ≈ 50%
        )
        
        # Baseline healthy cells: log-normal around 5e11 (±30% CV)
        healthy_baseline = rng.lognormal(
            mean=np.log(5.0e11),
            sigma=0.3,  # CV ≈ 30%
        )
        
        # Tumor growth rate: log-normal around 0.023/day (±40% CV)
        tumor_growth = rng.lognormal(
            mean=np.log(0.023),
            sigma=0.4,
        )
        
        # Healthy growth rate: log-normal around 0.015/day (±30% CV)
        healthy_growth = rng.lognormal(
            mean=np.log(0.015),
            sigma=0.3,
        )
        
        # PK/PD variability: clearance rates (±25% CV)
        len_clearance = rng.lognormal(mean=np.log(0.10), sigma=0.25)
        bor_clearance = rng.lognormal(mean=np.log(0.15), sigma=0.25)
        dara_clearance = rng.lognormal(mean=np.log(0.05), sigma=0.25)
        
        patients.append({
            "patient_id": i,
            "baseline_tumor_cells": float(tumor_baseline),
            "baseline_healthy_cells": float(healthy_baseline),
            "tumor_growth_rate": float(tumor_growth),
            "healthy_growth_rate": float(healthy_growth),
            "lenalidomide_clearance_rate": float(len_clearance),
            "bortezomib_clearance_rate": float(bor_clearance),
            "daratumumab_clearance_rate": float(dara_clearance),
        })
    
    return patients


def run_cohort(
    *,
    scenario: models.Scenario,
    n: int,
    regimen_params: Dict,
    user_id: int,
    seed: int = 42,
) -> Dict:
    """
    Run virtual cohort simulation.
    
    Generates n synthetic patients and simulates regimen outcomes
    for each, collecting aggregate statistics.
    
    Parameters:
        scenario: Clinical scenario context
        n: Number of patients in cohort
        regimen_params: Regimen parameters (doses, schedule, horizon)
        user_id: User ID for simulation attempts
        seed: Random seed for reproducibility
        
    Returns:
        Dict with keys:
            - n: cohort size
            - summaries: list of individual patient summaries
            - aggregates: mean/p95 efficacy, toxicity, recurrence rate
            - cohort_id: unique identifier for export
    """
    import uuid
    
    cohort_id = str(uuid.uuid4())
    patients = sample_patient_params(seed=seed, n=n)
    summaries = []
    
    for patient in patients:
        # Merge patient-specific params with regimen params
        merged_params = {**regimen_params, **patient}
        
        # Create simulation attempt
        attempt = models.SimulationAttempt.objects.create(
            scenario=scenario,
            user_id=user_id,
            parameters=merged_params,
            submitted=timezone.now(),
        )
        
        # Run simulation
        summary = attempt.run_model()
        
        # Augment with patient ID for tracking
        summary["patient_id"] = patient["patient_id"]
        summary["attempt_id"] = attempt.id
        summaries.append(summary)
    
    # Compute aggregates
    efficacies = [s.get("tumor_reduction", 0.0) for s in summaries]
    toxicities = [s.get("healthy_loss", 0.0) for s in summaries]
    recurrences = [1 if s.get("time_to_recurrence") is not None else 0 for s in summaries]
    
    # AUC totals
    auc_totals = []
    for s in summaries:
        auc = s.get("auc", {})
        total = sum([
            auc.get("lenalidomide", 0.0),
            auc.get("bortezomib", 0.0),
            auc.get("daratumumab", 0.0),
        ])
        auc_totals.append(total)
    
    aggregates = {
        "efficacy_mean": float(np.mean(efficacies)),
        "efficacy_p95": float(np.percentile(efficacies, 95)),
        "efficacy_p05": float(np.percentile(efficacies, 5)),
        "toxicity_mean": float(np.mean(toxicities)),
        "toxicity_p95": float(np.percentile(toxicities, 95)),
        "toxicity_p05": float(np.percentile(toxicities, 5)),
        "recurrence_rate": float(np.mean(recurrences)),
        "auc_mean": float(np.mean(auc_totals)),
        "auc_p95": float(np.percentile(auc_totals, 95)),
    }
    
    return {
        "cohort_id": cohort_id,
        "n": n,
        "summaries": summaries,
        "aggregates": aggregates,
        "seed": seed,
    }

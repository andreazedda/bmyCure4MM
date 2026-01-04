"""Design & simulation engine for rational MM therapy strategies.

This package is intentionally *modality-aware* and *evolution-aware*.
It is not a drug catalog and it is not a clinical recommender.

Core goals (MVP):
- Multi-compartment tumor model: bulk vs reservoir (stem-like/progenitor)
- Multi-clone representation with antigen profiles and drug sensitivities
- Targeting/activation engine with Boolean logic gates (AND/OR/NOT) + thresholds
- Mechanistic toxicity decomposition (on-target/off-tumor, off-target, payload, immune, marrow)
- Simple evolutionary dynamics + adaptive planning loop
- Report artifacts suitable for UI/plots (returned as JSON-ready dicts)
"""

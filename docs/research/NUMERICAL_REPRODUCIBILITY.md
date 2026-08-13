# Numerical Reproducibility Policy

## Governed baseline

Before dependency resolution changed, deterministic ODE and seeded design
outputs were characterized at Git `d99ff81dd73dec0245e3d6971cee9291c1c5572b`
with NumPy 2.0.2 and SciPy 1.13.1. The locked Python 3.11 environment reproduced
every value and the seeded report hash exactly. The machine-readable reference
is `tests/fixtures/numerical_baseline_v1.json`.

Run the fail-closed comparison with:

```bash
uv run python -m scripts.check_numerical_baseline
```

ODE scalar comparisons use relative tolerance `1e-7` and absolute tolerance
`1e-10`. The pure-Python seeded design report requires an exact canonical JSON
SHA-256. A zero process exit means only that these declared regression checks
passed; it is not evidence of calibration, clinical validity, or causal effect.

## Randomness

- Scientific stochastic entry points must accept an explicit seed.
- Use a local `random.Random(seed)` or NumPy `Generator`, not ambient global
  state.
- Monte Carlo sample count, seed, algorithm, and dependency identity must be
  retained with the result.
- Repeated runs with the same seed and version vector must be identical within
  the declared floating-point policy.
- Every persisted scientific execution records a seed in its immutable run
  manifest. Deterministic paths use the explicit default seed `0`.

## Floating point and parallelism

CI fixes `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and
`VECLIB_MAXIMUM_THREADS` to `1` for numerical gates. Parallel reduction order,
BLAS implementation, CPU architecture, compiler, NumPy, and SciPy can change
rounding. Comparisons outside the governed platform must retain the full
reproducibility report and may not claim bitwise identity.

The current baseline checks output tolerances, not solver correctness or
identifiability. Step-refinement properties and broader scientific invariants
belong to issue #23.

## Drift and invalidation

Any dependency, solver, tolerance, seed, equation, or parameter change must
produce an output-diff report. Unexpected drift blocks the change. Explained
drift requires scientific review, an explicit model-version decision, updated
regression evidence, and identification of prior runs that are no longer
comparable. Never overwrite this baseline merely to make CI pass.

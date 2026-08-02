# Solver C — complex-scaled spectral solver

**Status: built, validated, and independent.** This closes the gate that has
blocked the project for two sessions.

## Independence

Solver C uses **none** of: the Frobenius coefficient recurrence, the Gaussian
reduction, the Leaver continued fraction, the Hill determinant, Wynn
acceleration, or Solver A/B's spectral residual. Its root condition is the
smallest singular value of an independently assembled collocation matrix. This
is enforced by a test that greps the module source for the forbidden symbols
(`test_solver_c_uses_no_recurrence_machinery`).

## Why complex scaling, and why the earlier attempt failed

`docs/FAILED_APPROACHES.md` #3: boundary-factored collocation on the **real**
axis fails because after peeling `exp(iΩr)` the ingoing solution goes as
`exp(−2iΩr)`, which for `Im ω < 0` *decays*. Boundedness therefore does not
exclude it, and the matrix is singular everywhere.

Rotating the contour repairs exactly that. On

```
r(ρ) = r₊ + ρ e^{iθ},   ρ ∈ [0, L]
```

the ρ-dependent part of the outgoing exponent has real part
`ρ(−Ω_I cosθ − Ω_R sinθ)`, so the outgoing solution **decays** once

```
tan θ > −Ω_I / Ω_R
```

while the ingoing solution, carrying the opposite sign, **grows**. It is the
complex scaling — not boundedness — that enforces the radiation condition.

### The criterion predicts its own failure

For the static overtone (`ω = 0.3719 − 1.3226i`) the criterion requires
`θ > 74.3°`. At `θ = 60°` Solver C converges to a *wrong* root
(`0.5433 − 0.9092i`); at `θ = 82°` it returns `0.371859583 − 1.322609277i`,
matching the published value to `1.3e-8`. That failure is pinned as a
regression test.

## Contour safety

`r = 0` and `r = ±r₊`, `±r₋` all lie on the real axis, i.e. at `θ = 0` or `π`.
For `0 < θ < π` the contour never meets a singular point of the radial equation
and `Δ ≠ 0` along it. The branch of `Ω = √(ω²−μ²)` is fixed to `Re Ω ≥ 0`,
matching Solvers A/B. The angular eigenvalue is recomputed at every `ω` and
never frozen.

## Validation

| Test | ω (Solver C) | reference | difference |
| --- | --- | --- | --- |
| C1 static fundamental | `0.533834335 − 0.383373593i` | published | `2.2e-6` |
| C2 static overtone (θ=82°) | `0.371859583 − 1.322609277i` | published | `1.3e-8` |
| C3 static `l=1` | `1.016016912 − 0.362328023i` | published | `1.1e-9` |
| C4 two-spin `(0.2,0.3,μ=0.1)` | `1.681109025 − 0.347226588i` | HH CFM | `2.9e-5` |
| C5 two-spin massive `(0.4,0.2,μ=0.9)` | `1.822163662 − 0.311199991i` | HH CFM | `7.4e-5` |
| C6 two-spin `k=1 (0.3,0.1)` | `2.634399881 − 0.349098759i` | HH CFM | `1.8e-6` |

### Invariance (the test that separates resonances from continuum artifacts)

Static fundamental, deviation from the published value:

| varied | range | max deviation |
| --- | --- | --- |
| scaling angle θ | 45°–75° | `1.9e-6` |
| contour length L | 40–130 | `7.9e-8` |
| resolution N | 150–380 | `2.2e-6` |

A physical QNM stays put; a rotated-continuum eigenvalue would move with θ.

## Three-solver validation

`results/three_solver_validation.json` — **10/10 points agree below `1e-5`**:

| point | r₂ | Δ_AB | Δ_AC | Δ_BC |
| --- | --- | --- | --- | --- |
| static fundamental | 0 | `1.8e-12` | `7.1e-10` | `7.1e-10` |
| static overtone | 0 | `1.7e-07` | `1.6e-08` | `1.6e-07` |
| static `l=1` | 0 | `5.0e-15` | `3.1e-09` | `3.1e-09` |
| singly rotating massive | 0 | `8.1e-11` | `2.7e-09` | `2.7e-09` |
| two-spin r1 | 0.064 | `5.7e-12` | `2.4e-10` | `2.4e-10` |
| two-spin massive | 0.090 | `4.7e-12` | `8.8e-10` | `8.8e-10` |
| two-spin `k=1` | 0.032 | `2.6e-11` | `3.8e-10` | `3.8e-10` |
| equal spin | 0.067 | `9.4e-12` | `1.5e-10` | `1.5e-10` |
| **large r₂ = 0.140** | 0.140 | `2.5e-07` | `4.7e-11` | `2.5e-07` |
| **large r₂ = 0.190** | 0.190 | `3.3e-09` | `1.6e-09` | `4.9e-09` |

The last two matter most: they sit beyond the `r₂ ≳ 0.1` limit Huang and Huang
state for their own continued-fraction method, and a solver with no recurrence
in it now confirms the recurrence family there.

## Limitations

* `θ` must be chosen above the criterion for the mode being sought. It is not
  yet chosen automatically from a seed, so a caller hunting an unknown deep
  overtone must set it.
* Double precision only; no arbitrary-precision port.
* `r₂ ≳ 0.26` (the slowly convergent recurrence region) has **not** been tested.
* Not yet applied to long-lived massive branches.

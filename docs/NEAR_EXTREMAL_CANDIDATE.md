# The near-extremal candidate — adjudicated and refuted

**Classification: numerical artifact of the recurrence in its non-convergent
regime. Not an exceptional point.**

## The candidate

Newton on the augmented EP conditions reached

```
ω = 1.97826525 − 0.45044647i
s = 0.49777, δ = −0.11359, μ = 1.9   →   a = 0.38418, b = 0.61136
r₊ = 0.53410, r₂ = 0.43976, extremality = 0.00890, κ = 0.17202
|F| = 1.5×10⁻¹¹,  |dF/dω| = 4.3×10⁻⁶
```

`|dF/dω|` seven orders below its value at a simple root looked like a
near-double root. It is not.

## Three independent tests, all negative

**1. Solver C sees no root there.** At the candidate frequency Solver C's
`σ_min = 3.87×10⁻⁷`, while at a generic nearby frequency it is
`5.17×10⁻⁸` — *smaller*. Contrast 0.13. At a genuine QNM Solver C shows a
contrast of order `10¹¹`. There is no spectral condition being satisfied at ω_c.

**2. Solver A does not converge at these parameters.** Tracking an ordinary mode
at the same `(a, b, μ)`:

| depth | ω | movement |
| --- | --- | --- |
| 200 | `2.1404411 − 0.4903574i` | — |
| 400 | `2.1403019 − 0.4902030i` | `2.1e-4` |
| 800 | `2.1404805 − 0.4901583i` | `1.8e-4` |
| 1600 | `2.1398633 − 0.4903641i` | `6.5e-4` |
| 3200 | `2.1395202 − 0.4900188i` | `4.9e-4` |

The root **drifts and does not settle**, while the CF residual sits at `~10⁻¹³`
throughout. That is the failure mode already documented in
`docs/PRECISION_VALIDATION.md`: the residual measures the *truncated* equation,
not the truncation error.

**3. The three solvers agree only to 3.5×10⁻³ here.**

| pair | difference |
| --- | --- |
| A–B | `1.2×10⁻³` |
| A–C | `4.1×10⁻³` |
| B–C | `3.5×10⁻³` |

Solver C is internally stable (contour length 20–60, resolution 280–650 all give
`2.1423 ± 1×10⁻³`), so the spread is not Solver C's error alone.

**A repeated-root diagnostic of `4.3×10⁻⁶` is three orders of magnitude below
the frequency uncertainty of `3.5×10⁻³`.** It carries no information.

## Mechanism

Near extremality the spurious singularity `x(−r₊) = 2r₊/(r₊+r₋)` approaches the
unit circle (`r₊ − r₋ = 0.094` here), so the Frobenius series' radius of
convergence degenerates — the mechanism already identified in
`docs/ROTATING_VALIDATION.md` as the reason Huang–Huang's own continued fraction
degrades past `r₂ ≳ 0.1`. The continued fraction then returns depth-dependent
roots whose residual is small but whose value drifts, and Newton on
`(F, dF/dω)` happily finds a spurious stationary point in that drift.

Ruled out as the cause: `|σ| = 1.69` at the candidate, only ~2× its ordinary
value, so the horizon exponent is not blowing up; and Solver C's resolution,
since its answer is stable under refinement.

## What is and is not established

* **Established:** this candidate is not evidence of an exceptional point, and
  the region `r₂ ≈ 0.44` is outside the validated domain of *all three* solvers
  (best achievable agreement `3.5×10⁻³`).
* **Not established:** that no EP exists near extremality. The region is simply
  unresolvable at present accuracy. Settling it needs a method that resolves the
  near-horizon scale `r₊ − r₋ = 0.094` directly — multidomain spectral matching
  or near-horizon matched asymptotics. Neither was built.

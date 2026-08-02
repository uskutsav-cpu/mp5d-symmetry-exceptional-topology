# Final solver domains

`z₋` is the **inner** horizon radius squared — the near-extremality difficulty
parameter. Extremality is `M − (|a|+|b|)²`. Both are properties of the
background, not of the mode.

| solver | module | formulation | recurrence-free | validated domain | known failure |
| --- | --- | --- | --- | --- | --- |
| **A** | `radial/qnm.py::solve_qnm_cf` | Leaver continued fraction, Gaussian-reduced | no | `z₋ ≤ 0.20`, extremality `≥ 0.05` | drifts `2–6e-4` per depth doubling at `z₋ = 0.44` |
| **B** | `radial/qnm.py::solve_qnm_hill` | raw recurrence + Hill determinant + Wynn | no | as A | shares A's recurrence family |
| **C** | `radial/solver_c.py` | exterior complex scaling, Chebyshev collocation | **yes** | ordinary modes to `z₋ = 0.19`; `θ ∈ [45°,75°]`, `L ∈ [40,130]`, `N ∈ [150,380]` | no root found at `z₋ = 0.44` (contrast 0.13) |
| **D** | `radial/multidomain_near_horizon.py` | multidomain near-horizon + exterior scaling | **yes** | ordinary modes | does not converge at `z₋ = 0.44`; wanders in a ball of radius `~2e-3` |
| **HP** | `radial/highprec.py` | full mpmath chain (coefficients included) | no | arbitrary precision, same parameter domain as A | slow |

## Independence

A and B share the Frobenius recurrence and are **not** independent of each
other. C and D are the recurrence-free pair; C's independence is enforced
structurally by a test that greps the module for `reduce_to_three_term`,
`continued_fraction`, `hill_determinant`, `wynn_epsilon` and `recurrence_row`
and fails if any appears.

**Rule:** any major claim requires agreement between at least one recurrence
solver and at least one recurrence-free solver. Two wrappers around the same
discretization do not count.

## Measured agreement

Three solvers agree below `1e-5` on 10/10 benchmark points, including
`z₋ = 0.140` and `z₋ = 0.190` — beyond the `r₂ ≳ 0.1` limit Huang–Huang state
for their own continued-fraction method.

## The near-extremal region, honestly

At `z₋ = 0.44` (extremality `0.0089`) **no solver converges** and the best
achievable three-solver spread is `3.5e-3`. Two candidate explanations were
measured and refuted:

* **near-horizon scale** — refuted: Solver D resolves `r₊ − r₋` by construction
  and still fails;
* **double-precision conditioning** — refuted: `cond(M) = 5.3e11` at the hard
  point versus `3.2e11` at an easy one that converges to `1e-9`.

No claim of any kind is made in this region. It is outside every declared
domain, and the one candidate found there is recorded as falsified (C34).

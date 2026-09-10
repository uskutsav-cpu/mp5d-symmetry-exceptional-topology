# The EP2 verification gate

A candidate is **not** an exceptional point until every applicable condition
below passes. Frequency proximity is not one of the conditions — it is only a
reason to open an investigation.

Each condition names the routine that implements it, and each routine is
validated against closed-form synthetic systems in
`tests/unit/test_ep_machinery_synthetic.py` — including negative controls, so a
detector that could never fire is caught.

## Gate conditions

| # | condition | routine | synthetic validation |
| --- | --- | --- | --- |
| G1 | scale-free root separation `\|a₁/a₂\| → 0` under refinement | `diagnostics.root_separation` | recovers the true root distance to `1e-6`; invariant under rescaling by `1e7 e^{3z}` |
| G2 | argument-principle zero count is **2** inside a contour that shrinks with the candidate | `diagnostics.zero_count` | counts zeros and poles correctly; raises when unresolved |
| G3 | both roots located without a starting guess | `verification.roots_in_disc` | finds both roots of a pair separated by `1e-3` |
| G4 | monodromy: one loop transposes the branches | `verification.monodromy` | swaps for `λ²−t`; **does not** swap for the avoided crossing `λ²−(t²+g²)` |
| G5 | two loops restore the assignment | `verification.track_roots_on_loop` | verified |
| G6 | Puiseux: splitting `∼ t^{1/2}`, preferred over `t^1` | `verification.fit_puiseux` | exponent `0.5 ± 1e-6`; rejects sqrt for linear splitting |
| G7 | Jordan chain solvable, geometric multiplicity < algebraic | `verification.jordan_chain_matrix` | fires on a Jordan block; **silent** on semisimple and on block-diagonal coincidence |
| G8 | stable under resolution, depth and precision refinement | solver-level | — |
| G9 | confirmed by a **recurrence-free** solver (C or D) | cross-solver | — |
| G10 | inside the validated solver domain, or independently validated there | domain guard | — |

## Why G9 and G10 are not optional

The one previous candidate that looked convincing (`|dF/dω| = 4.3e-6`) failed
exactly these two. It sat at `z₋ = 0.44` with extremality `0.0089`, outside
every solver's validated domain, and the recurrence-free solvers showed no root
there at all. It is recorded as **falsified** (C34), not as "promising".

## Why G1 replaces the old `|dF/dω|` criterion

`|dF/dω|` at a root is not invariant: the spectral condition may be multiplied
by any nonvanishing analytic function without changing its zeros, and that
scales `|dF/dω|` arbitrarily. `|a₁/a₂|` is the ratio of the first two Taylor
coefficients and the prefactor cancels, so it estimates the actual distance to
the partner root. See
`test_root_separation_is_invariant_under_rescaling`.

## Higher order

An EP3 additionally requires: three coalescing branches, zero count 3,
`a₁ = a₂ = 0` with `a₃ ≠ 0`, cube-root splitting, a **3-cycle** monodromy with
no fixed point, and Jordan-chain length 3. A three-mode cluster on a plot is not
an EP3.

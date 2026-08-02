# Exceptional-point search — status

**No exceptional point has been verified. None is claimed.**

## What was built

`src/mp5d/exceptional/ep_solver.py` — the augmented EP2 conditions

```
F(ω, p) = 0,      dF/dω = 0
```

solved by Newton in `(ω, p₁, p₂)`. Both conditions are complex, so the system is
4 real equations in 4 real unknowns (ω plus **two** of `s, δ, μ`). EP2 is
therefore codimension 2, and in the three-parameter space the exceptional set is
generically a one-dimensional locus.

**The total derivative is exact.** `dF/dω` includes `∂F/∂Λ · dΛ/dω` because
`cf_value` recomputes `Λ(ω)` from the angular solver at every evaluation. Λ is
never frozen — freezing it is the standard way to manufacture spurious EPs.

Verified behaviour: at a known QNM, `|F| = 4×10⁻⁸` and `|dF/dω| = 4.6`
(simple root); at a generic frequency `|F| = 0.84`.

## Where we looked

Sector `(m₁,m₂) = (1,1)`, branches `l = 2` (k=0) and `l = 4` (k=1), plus the
radial overtone ladder `N = 0,1,2`, continued in `μ ∈ [0.1, 3.0]` at
`a = 0.2, b = 0.3`.

* `l = 2` vs `l = 4`: separation never below **0.84**. No collision.
* `N = 0` vs `N = 1`: separation shrinks with μ (0.72 → 0.60) but never closes.
* Near `μ ≈ 2.0` the `l = 2` branch crosses toward `ω² = μ²` and jumps to the
  quasi-bound regime with `−Im ω ≈ 0.003`. That is the long-lived behaviour
  Huang–Huang report, and it is a **branch point, not an EP**.

## One candidate, rejected

Newton reached `|F| = 1.5×10⁻¹¹` with `|dF/dω| = 4.3×10⁻⁶` — seven orders below
the value at a simple root, the signature of a near-double root — at

```
s = 0.4978, δ = −0.1136, μ = 1.9   →   a = 0.384, b = 0.611
```

**Rejected**, because that point has `r₂ = 0.44` against a validated recurrence
domain of `r₂ ≲ 0.20`, and extremality `M − (|a|+|b|)² = 0.0089`, i.e. essentially
extremal — exactly where the asymptotic tail breaks down (`c → 0`) and where
Solver C has never been tested. Newton also stagnated near `10⁻⁶` instead of
converging. Recorded in `results/rejected_collisions.json` with what would settle
it.

This is neither evidence for nor against an EP. It is a candidate in an
unvalidated region.

## Bounded negative result

Scanning the **safe** domain (`r₂ ≤ 0.20`, extremality `≥ 0.05`) over
`s ∈ {0.15, 0.25, 0.35}`, `δ ∈ {−0.15 … 0.15}`, `μ ∈ {0.1 … 1.7}` — 74
converged points in sector `(1,1)`, `l = 2`, `N = 0`:

| quantity | value |
| --- | --- |
| min `|dF/dω|` | **1.02×10⁻²** |
| median | 3.82 |
| max | 7.31 |

An EP requires `|dF/dω| → 0`. The minimum over the box is bounded away from
zero by two orders of magnitude, so **there is no EP2 in this box** for this
sector, branch and overtone. `results/negative_regions.json`.

This is a genuine bounded exclusion, not "nothing showed up in a coarse scan":
the diagnostic is the actual EP condition, evaluated at converged roots.

## What this does not establish

* Only **one sector**, one `l`, one overtone, one parameter box.
* Nothing about `r₂ > 0.20` or near-extremal parameters — where the one
  interesting candidate appeared and where the solvers are unvalidated.
* No Jordan chain, Puiseux, or monodromy machinery was built, because no
  candidate survived to need it.

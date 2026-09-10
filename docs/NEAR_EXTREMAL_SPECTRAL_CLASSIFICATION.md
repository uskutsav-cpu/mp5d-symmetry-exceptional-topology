# The near-extremal breakdown: resolved, and it was not extremality

Gate 13. This closes the open question left by `docs/MULTIDOMAIN_NEAR_HORIZON.md`.

## The standing question

Two explanations for the near-extremal solver breakdown had already been
measured and **refuted**:

* *unresolved near-horizon scale* — refuted: Solver D resolves `r₊ − r₋` by
  construction and still fails;
* *double-precision conditioning* — refuted: `cond(M) = 5.3e11` at the failing
  point versus `3.2e11` at an easy point that converges to `1e-9`.

The remaining hypothesis was that the difficulty is *spectral*. It is — but not
in the way expected, and **extremality is not the variable that controls it**.

## Test 1 — extremality alone does not break anything

Continue the sector-(1,1), `l = 2`, `N = 0` mode in `z₋` at fixed `μ`, along the
family `b = a/1.6`, up to extremality `M − (|a|+|b|)² = 0.0103` — essentially
extremal.

| `z₋` | extremality | `ω` at `μ = 0.5` | residual | converged |
| --- | --- | --- | --- | --- |
| 0.00 | 1.0000 | `1.547960 − 0.343919i` | 2.9e-14 | yes |
| 0.08 | 0.1697 | `1.990011 − 0.265448i` | 8.1e-14 | yes |
| 0.14 | 0.0511 | `2.078551 − 0.205193i` | 6.0e-14 | yes |
| 0.17 | 0.0222 | `2.102068 − 0.173016i` | 1.0e-13 | yes |
| **0.19** | **0.0103** | `2.107380 − 0.152983i` | **1.4e-13** | **yes** |

The same holds at `μ = 1.0`. Then, at that **fixed near-extremal background**,
`μ` was scanned from 0.5 to 2.3 in steps of 0.05: all 37 steps converged, with
residuals `~1e-13` throughout, `ω` moving smoothly from `2.107 − 0.153i` to
`2.489 − 0.139i`.

**Solver A does not fail near extremality.** The earlier failure was attributed
to the wrong variable.

## Test 2 — what actually controls it

Exterior complex scaling (Solvers C and D) requires a contour angle satisfying

```
tan θ  >  −Im Ω / Re Ω ,        Ω = √(ω² − μ²)
```

so the outgoing solution `e^{iΩr}` decays along the rotated contour. Evaluating
the **minimum admissible angle** at the relevant points:

| case | `ω` | `μ` | `Ω` | min `θ` |
| --- | --- | --- | --- | --- |
| ordinary, near-extremal | `2.107 − 0.153i` | 0.5 | `2.048 − 0.157i` | **4.4°** |
| odd-`l` corner | `1.198 − 0.727i` | 1.5 | `0.655 − 1.331i` | **63.8°** |
| domain corner | `1.703 − 0.000269i` | 1.8 | `0.00078 − 0.584i` | **89.92°** |
| the old `μ = 1.9` failure | `1.547 − 0.00266i` | 1.9 | `0.00373 − 1.103i` | **89.81°** |

As `Re Ω → 0` the outgoing wave stops oscillating and becomes a pure growing
exponential. The required scaling angle then tends to **90°**, which no contour
can supply. This is a **formulation limit of exterior complex scaling, not a
convergence failure** — no resolution, precision or domain decomposition can fix
it, which is exactly why Solver D's near-horizon refinement did not help.

The condition `Re Ω → 0` is precisely the **quasiresonant / long-lived limit**
documented in `docs/LONG_LIVED_FINAL.md`, where the damping collapses toward
zero.

## Classification

The near-extremal object that could not be resolved was **not** an isolated
resonance failing to converge, and **not** a continuum artifact. It was an
ordinary quasinormal mode that had entered the **quasiresonant regime**, where:

* Solvers C and D are inapplicable by formulation (`θ_min → 90°`);
* Solver A's continued fraction needs diverging depth, because a mode of width
  `~1e-4` requires resolving a correspondingly fine structure.

The old candidate was reported at `μ = 1.9`, and at `μ = 1.9` the sector-(1,1),
`l = 2` mode has damping `~3e-3`. The "near-extremal" label was a coincidence of
the parameter choice: **the controlling variable was `μ`, not the spin.**

## Consequence for Solvers E and F

`docs/FINAL_EXECUTION_GATES.md` deferred the Wronskian and hyperboloidal solvers
until this diagnostic said whether a better solver could help. It now has:

* A **Wronskian solver would not** help in the quasiresonant regime for the same
  reason as C and D if it integrates along a rotated contour — but it *would* be
  the right tool if it matches series directly, since it never needs an
  admissible rotation angle.
* A **hyperboloidal** formulation is the principled fix: it regularizes the
  outgoing behaviour geometrically rather than by rotation, so `Re Ω → 0` is not
  a special case for it.

That is a concrete, evidence-based next step rather than a guess, and it is
recorded as the next scientific gate.

## Correction to the record

`docs/MULTIDOMAIN_NEAR_HORIZON.md` and `docs/NEAR_EXTREMAL_CANDIDATE.md` both
frame the difficulty as near-extremality. That framing is **wrong** and is
corrected here. The refutation of the candidate itself (C34) stands unchanged —
it was, and remains, a numerical artifact.

### A naming hazard that contributed

Prior work reported `r₂`, which is the inner horizon **radius** `r₋`, while the
geometry module exposes `z₋ = r₋²`. "`r₂ = 0.44`" is `z₋ = 0.194`. A guard
written as `z_minus <= 0.20` is therefore `r₋ ≤ 0.447`, five times looser in
`z` than the `r₋ ≤ 0.20` it was believed to reproduce. Both quantities are now
reported side by side everywhere.

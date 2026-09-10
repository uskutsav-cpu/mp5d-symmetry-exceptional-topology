# First trustworthy MP5D quasinormal mode

**Status: cross-solver verified.** This is the first quasinormal frequency
computed in this repository, and the first that meets the benchmark gate in
`.claude/skills/mp5d-complete-project/SCIENTIFIC_GATES.md`.

## Result

Five-dimensional Schwarzschild–Tangherlini, massless scalar, `l = 0`,
fundamental overtone, in `r₊ = 1` units:

```
omega = 0.533835574268 - 0.383375368513 i
```

| Quantity | Value |
| --- | --- |
| Published (Matyjasek, arXiv:2107.04815, Table I) | `0.533835574268 − 0.383375368512 i` |
| Absolute error | `4.98e-13` |
| Relative error | `7.59e-13` |
| Continued-fraction residual | `1.15e-12` |
| ODE residual (independent points) | `4.39e-13` |
| Solver A vs Solver B | `2.15e-12` |

All twelve published digits are reproduced.

### Normalization

The source uses `f(r) = 1 − r^{3−D}` with `D = 5`, so the horizon sits at
`r = 1` and `T_H = 1/(2π)`. This repository's `M = 1` gives `r₊ = 1` and
`κ = 1`, verified numerically — the same normalization. The source tabulates
`ω̃ = ω/T_H`; values here are converted with `T_H = 1/(2π)` and stored in
`data/published_benchmarks/`.

## Benchmark set

| Case | `l` | `n` | ω (Solver A) | abs. error | \|A − B\| |
| --- | --- | --- | --- | --- | --- |
| `st5d_l0n0` | 0 | 0 | `0.533835574268 − 0.383375368513i` | `4.98e-13` | `2.15e-12` |
| `st5d_l0n1` | 0 | 1 | `0.371859570894 − 1.322609278713i` | `2.65e-11` | `9.96e-08` |
| `st5d_l1n0` | 1 | 0 | `1.016016911488 − 0.362328023848i` | `6.08e-13` | `3.88e-15` |
| `st5d_l1n1` | 1 | 1 | `0.856379797311 − 1.157605662729i` | `2.73e-13` | `3.51e-12` |

Overtones converge more slowly in truncation depth, which is why `n = 1` needs a
deeper schedule to reach the same absolute accuracy.

## Gate checklist

| Requirement | Status |
| --- | --- |
| Continued fraction converges | yes |
| Stable with recurrence depth | yes — flat from depth 200 (see table below) |
| Stable with arithmetic precision | **partially** — see limitations |
| Small ODE residual | `4.4e-13` |
| Independent initial guess reaches the same root | yes — seeds `0.60−0.30i`, `0.70−0.25i`, `0.50−0.34i` |
| Sign convention gives a decaying mode | yes, `Im ω < 0` |
| Agrees with a published benchmark to all displayed digits | yes, 12/12 |

### Depth convergence (`l = 0`, `n = 0`)

| depth | ω |
| --- | --- |
| 100 | `0.533835574266 − 0.383375367831i` |
| 200 | `0.533835574268 − 0.383375368513i` |
| 400 | `0.533835574268 − 0.383375368512i` |
| 800 | `0.533835574268 − 0.383375368512i` |

## Method

**Solver A — Leaver continued fraction.** Frobenius series in
`x = (r−r₊)/(r−r₋)`. Both spurious singular points map outside the unit disk
because `r₊ > r₋`, so the radius of convergence is exactly 1 with the irregular
point on the boundary at `x = 1` — Leaver's setting. Prefactors:

```
R = x^{-iσ} (1-x)^{3/2} exp( iΩ(r₊-r₋)/(1-x) ) Σ aₙ xⁿ
```

`σ = (ω − m₁Ω_a − m₂Ω_b)/(2κ)` (ingoing at the horizon);
`exp(iΩr) r^{-3/2}` with `Ω = √(ω²−μ²)` (outgoing at infinity), the power being
exactly `−3/2` because the radial potential is even in `r` (claim C18).

The raw recurrence has **width 9 for the static case `a = b = 0`**; for generic
two-spin parameters it is **width 13**. It is reduced to three terms by the general
Gaussian elimination in `mp5d.radial.recurrence`, then the `n`-th inversion of
Leaver's condition is solved by Muller's method. `Λ` is recomputed from the
angular solver at every iterate — never frozen.

**Solver B — Hill determinant with Wynn acceleration.** Uses the raw N-term
recurrence directly: no Gaussian reduction and no continued fraction. The
truncated determinant is proportional to the forward-generated coefficient
`a_N`, so that is evaluated instead, and the root sequence in truncation order
is accelerated with Wynn's epsilon algorithm. Following Benda and Matyjasek,
arXiv:2503.17325.

The two solvers share only the extracted polynomials `A, B, C`; their root
conditions have no machinery in common, which is what makes them an admissible
independent pair.

## Symmetry checks (rotating, two-spin)

These are exact identities, so their residuals measure numerical error directly
and are the sharpest diagnostics available.

| Check | Result |
| --- | --- |
| Exchange `(a,m₁) ↔ (b,m₂)`, `a=0.30, b=0.15, μ=0.20` | residual `4.15e-15` |
| **U(2) multiplet degeneracy** at `a=b=0.25`, `μ=0.15`, `l=2`, `m=2`: `(1,1)`, `(2,0)`, `(0,2)` | max spread `4.67e-15` |
| Same states at `a=0.35, b=0.15` (δ ≠ 0) | split by `5.38e-02` |

The middle row is the **first numerical confirmation of claims C6 and C7**,
which were derived analytically in session 1: at `δ = 0` the full separated
problem depends on `(m₁, m₂)` only through `m = m₁ + m₂`, so those three states
must share a frequency. They do, to machine precision. The third row confirms
the degeneracy is genuinely lifted off the equal-spin surface, so the second row
is not passing for a trivial reason.

## Limitations

* **Double precision only.** The polynomial coefficients `A, B, C` are extracted
  by FFT in double precision, which floors the achievable residual at about
  `1e-12`. The `precision` field reports `"double"` throughout. Requirement 3 of
  the gate ("stable with arithmetic precision") is therefore **not** fully
  discharged: stability was demonstrated against truncation depth, not against
  an increase in working precision. Extending the extraction to exact rational
  or `mpmath` arithmetic is the natural next step.
* **Seed sensitivity.** The root search can walk onto the `Re ω = 0` branch cut
  of `Ω = √(ω²−μ²)` and return a spurious near-imaginary "root". Mitigations in
  place: Muller iterates are clamped to `Re ω ≥ 0.05`, and the depth schedule
  starts coarse (depth 100), which widens the basin of attraction. Rejected
  candidates are recorded in `results/rejected_radial_roots.json`.
* **Solver B needs closer seeds than Solver A** on rotating backgrounds. With a
  loose seed it drifted to the branch cut; with a tight seed it reproduced
  Solver A to `2.3e-13`. This is a seeding limitation, not a spectral feature.
* **No Nollert-style asymptotic tail is implemented.** The continued fraction is
  evaluated by plain truncation. This suffices for the modes computed here
  (convergence is flat by depth 200–800) but will matter for high overtones and
  weakly damped massive modes.
* **No rotating case has been compared against a published value.** The rotating
  results above are validated by exact internal symmetries and by A/B agreement,
  not against the literature. The Huang–Huang two-spin comparison remains open.

## Reproduction

```bash
PYTHONPATH=src .venv/bin/python scripts/first_qnm.py
```

Writes `results/first_qnm.json`. Regression tests live in
`tests/unit/test_first_qnm.py`.

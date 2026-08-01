# Numerical validation

## Done

| Check | Method | Result |
| --- | --- | --- |
| `√\|g\|` and `det B` | exact symbolic | passes exactly |
| Separation of the frequency term | exact symbolic (`∂²/∂r∂u ≡ 0`) | passes exactly |
| Separated potentials vs. implementation | exact symbolic | passes exactly |
| Schwarzschild-Tangherlini radial limit | exact symbolic | passes exactly |
| `T_H = 1/(2π r₊)` at `a = b = 0` | closed form | passes |
| Extremal limit `κ → 0` at `M = (\|a\|+\|b\|)²` | closed form | passes |
| Horizon Vieta relations `z₊z₋ = a²b²`, `z₊+z₋ = M−a²−b²` | closed form | passes |
| Exchange symmetry of horizon data | numeric | passes |
| Angular `Â = l(l+2)` at `c² = 0` | numeric, `l ≤ 3`, all multiplet members | `< 1e-10` |
| Angular exchange relation `Â(m₁,m₂,c²) = Â(m₂,m₁,−c²) − c²` | numeric, real and complex `c²` | `< 1e-9` |
| Angular spectral vs. independent finite difference | two discretizations sharing no machinery | `< 5e-3` on the lowest 3 branches at FD `N = 1500` (FD is 2nd-order; this is its accuracy, not the spectral method's) |
| Angular truncation convergence | `N` refinement | machine precision by `N ≈ 10` for `\|c²\| ≲ 100`; monotone geometric decay verified up to `\|c²\| ~ 4×10³` |

Total: 46 tests, all passing.

## Not done

Everything radial. In particular **no QNM frequency has been computed**, so none
of the specification's §8 baseline benchmarks (Schwarzschild-Tangherlini QNM,
singly-rotating points, Huang-Huang two-spin points, long-lived branch,
second-method cross-check) has been attempted yet.

## Declared tolerances (to be applied once the radial solver exists)

* benchmark agreement with published values: relative `1e-4` on `Re ω` and
  `Im ω`, or the published precision, whichever is looser;
* convergence: relative change `< 1e-8` under doubling the radial resolution;
* cross-solver agreement: relative `1e-8`;
* exchange-symmetry check: relative `1e-10` (this is an exact identity, so it
  measures numerical error directly and is the most sensitive diagnostic
  available).

# Reuse map

Rule: wrap or adapt existing tested components; implement only the
Myers-Perry-specific operators and missing interfaces.

Legend for "verified": **yes** = installed/inspected in this environment;
*claimed* = known from documentation but not yet exercised here.

| Capability | Existing source | Reuse directly? | Missing component | Action | Verified |
| --- | --- | --- | --- | --- | --- |
| Interval / ball arithmetic | `python-flint` (Arb) 0.9.0 | yes | none | installed; use for Stage F certification | **yes** |
| Arbitrary-precision scalars, root finding | `mpmath` 1.3.0 | yes | none | use for precision escalation in the radial CF | **yes** |
| Exact symbolic derivation of the separated system | `sympy` 1.14 | yes | none | used; drives `tests/unit/test_separation.py` | **yes** |
| Dense eigensolve, SVD, condition numbers | `numpy`/`scipy` (LAPACK) | yes | none | used for the angular operator; will carry the radial collocation matrices | **yes** |
| Jacobi polynomial recurrence | classical three-term recurrence | yes | none | implemented directly (10 lines); a library dependency is not warranted | **yes** |
| Kerr QNM reference values | `qnm` (Stein, arXiv:1908.10377) | no — Kerr only, `D = 4` | 5D MP has no analogue | not applicable; do **not** wrap | *claimed* |
| Black Hole Perturbation Toolkit | BHPToolkit | no — 4D Kerr/Schwarzschild | 5D MP radial/angular operators | reuse conventions and method descriptions only | *claimed* |
| Nonlinear eigenvalue solvers | SLEPc `NEP` / `slepc4py` | candidate | build not yet attempted on this machine | evaluate before writing any bespoke NEP solver | not yet |
| Continuation / bifurcation tracking | `BifurcationKit.jl` | candidate | Julia not installed on this machine | either install Julia or implement pseudo-arclength directly (small systems) | not yet |
| Certified root finding | `HomotopyContinuation.jl`, Krawczyk | candidate | Julia not installed | prefer `python-flint` + Krawczyk in Python to avoid a second toolchain | not yet |
| Pseudospectra | EigTool / Chebfun | no — MATLAB | — | compute resolvent norms directly from SVDs; the matrices are small | *claimed* |
| Spectral collocation | Chebfun / Dedalus / `eigentools` | partial | our operator is a small dense nonlinear eigenproblem, not a PDE stack | implement Chebyshev differentiation matrices directly | not yet |

## Notes

* `uv` is not installed on this machine and `julia` is absent. `requirements.lock`
  (a `pip freeze`) substitutes for `uv.lock`. If any Julia component is adopted,
  the toolchain decision must be recorded here first.
* The angular sector needed **no** external special-function library: the
  equal-spin problem is exactly Jacobi, and the spheroidal correction is exactly
  tridiagonal in that basis. This is a genuine simplification, not a shortcut.

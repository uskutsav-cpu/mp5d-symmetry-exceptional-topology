# Session handoff

## Completed (session 1, 2026-08-01)

Stage A essentially complete; a substantial and unplanned part of Stage C
(equal-spin symmetry analysis) turned out to be reachable analytically and was
completed first, because it sharply narrows the numerical search.

* Repo, venv, dependency lock, CI (pytest + gitleaks), package skeleton.
* **Separation re-derived from the metric symbolically**, not transcribed.
  Conventions in `docs/CONVENTIONS.md` are enforced by CI.
* Geometry module (horizons, `Ω_a`, `Ω_b`, `κ`, `T_H`, `s`/`δ`) + tests.
* Angular sector **complete and cross-verified**: exactly-tridiagonal
  complex-symmetric Jacobi construction + independent finite-difference solver.
* Exact symmetry results C4–C10 in `docs/CLAIM_LEDGER.md`, including a
  structural **no-go** (C8) that removes the specification's primary target
  scenario §2.1 and redirects the search (`docs/SYMMETRY_STRUCTURE.md` §5).

## Radial sector: structure done, solver not

Derived and pinned (C15-C18): `r = 0` is an ordinary point for `ab != 0`; the
problem is confluent-Heun type; the horizon exponent is
`(ω − m₁Ω_a − m₂Ω_b)/(2κ)`; the potential is exactly even in `r` so the
asymptotic exponent is exactly `−3/2` with no Coulomb phase.

That fixes every ingredient a radial method needs *except* the eigenvalue
solver itself.

## Not started

* **Radial eigenvalue solver.** This is the critical path. Nothing downstream (branch
  atlas, degeneracy detection, exceptional sets, certification) can begin
  without it. No QNM frequency has been computed yet (claim C14).
* Baseline benchmark validation against Huang–Huang and the singly-rotating
  literature.
* Author-code recovery (`docs/AUTHOR_CODE_REQUEST.md` drafted, not sent).

## Next actions, in order

1. Implement the radial solver. Recommended: Leaver-type continued fraction in
   `z = r²` (regular singular points at `z = 0, z₊, z₋`; irregular point at
   infinity of half-integer rank, so the expansion variable should be `r`, not
   `z`) **plus** an independent Chebyshev collocation on a compactified
   coordinate. The collocation form is required regardless, because the
   exceptional-point analysis needs left/right eigenvectors and condition
   numbers, which a continued fraction cannot supply.
2. Reproduce the minimum benchmark set (§8 of the specification).
3. Build the equal-spin branch atlas at fixed `(m₁, m₂)` — note that by C6 the
   atlas is indexed by `(l, m, N)` at `δ = 0`, which is much smaller than a
   naive `(n, m₁, m₂, N)` enumeration.
4. Search within diagonal sectors `m₁ = m₂` first (C10).

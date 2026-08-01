# Session handoff

## Session 5 (2026-08-01) — asymptotic tail; Solver C NOT built

### Achieved

* PR #2 merged; work continues on `research/solver-c`.
* **Tail structure derived.** Characteristic equation `A(1/R) = 0`; `x = 1` is a
  multiple root (`A` multiplicity ≥5, `B`,`C` exactly 4), so the `n^{-3/2}`
  balance reads `−u₁A'(1) = 0` with `A'(1) = 0` and half-integer powers survive:
  `R_n = 1 + u₁n^{-1/2} + …`, `u₁ = −√(−2c)`, `c = iΩ(r₊−r₋)`.
* Confirmed two independent ways: the exponent test (`n^{1/2}` scaling
  converges, `n^1` diverges) and the `exp(−k√N)` depth signature.
* **Tail implemented** at leading order, opt-in via `tail_order`; higher orders
  raise rather than返回 a wrong tail. Measured **~1.35× depth reduction**
  (+1.4 digits at depth 200).
* 127 tests passing.

### NOT done — the session's central deliverable

**Solver C was not built.** There is still no independent non-recurrence solver.
Everything in the difficult region continues to rest on one recurrence family
plus Huang–Huang's published matrix values. Gate items 3–10 are therefore open:
static fundamental/overtone via C, two rotating modes via C, large-`r₂` via C,
three-way A/B/C agreement, long-lived branch, and the revised validity map.

### Blocker discovered

The Gaussian reduction degrades past `n ~ 500–800`: the backward ratio recursion
collapses to a spurious, problem-independent `R_n − 1 = 2/n` (identical for two
unrelated parameter sets). Clean window for asymptotics is `n = 200–400`. This
blocks numerical extraction of `u₂` and caps the achievable tail order. It does
**not** affect the CF values or the benchmark agreements.

Untested fixes: raise working precision inside the reduction independently of
the solve precision; or restructure the elimination to shorten the recursive
chain.

### Next, in order

1. Fix the reduction degradation, then derive/extract `u₂`, `u₃`.
2. Build Solver C (complex-contour Wronskian matching) — horizon Frobenius
   series outward, outgoing asymptotic series inward along a rotated contour,
   match logarithmic derivatives at a complex interior point. Must not reuse the
   recurrence, the reduction, the CF, or the Hill determinant.
3. Only then: three-way validation, long-lived branch, revised validity map.

No branch atlas and no EP search until that gate closes.

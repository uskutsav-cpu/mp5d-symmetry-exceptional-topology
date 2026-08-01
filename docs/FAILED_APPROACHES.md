# Failed approaches

Kept deliberately, per specification §25 Outcome D: failures are preserved, not
overwritten.

## 1. Full 5x5 symbolic metric inversion (session 1)

Building the inverse metric with `sympy` on the full 5x5 matrix with symbolic
`a, b, M, r, θ` did not complete in 10 minutes. **Cause:** `sympy.simplify` on
trigonometric rational expressions of this size, not the inversion itself.

**Fix that worked:** exploit the block structure — `(r, θ)` are diagonal, only
the `(t, φ, ψ)` 3x3 block needs inverting — and substitute `u = cos²θ` so every
entry is a *rational* function. Then use `cancel`/`together` only, never
`simplify`. Runtime dropped from >600 s to 2.2 s for the entire test module.

Lesson recorded because it will recur: in this project, always move to `u` and
stay in rational arithmetic.

## 2. Guessing the separated radial potential from memory (session 1)

An initial hand-written candidate for the radial potential disagreed with the
derived expression. Rather than patch it, the potential was extracted directly
from the derivation by partial-fractioning the exact expression in `u`. All
coded potentials now come from that extraction and are pinned by CI.

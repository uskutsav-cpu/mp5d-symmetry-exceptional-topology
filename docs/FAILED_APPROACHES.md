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

## 3. Boundary-factored Chebyshev collocation for the radial QNM problem (session 2)

**Attempted.** `src/mp5d/radial/collocation.py`: peel
`R = e^{iΩr}(r−r₊)^{−iσ}r^{−3/2+iσ} f(x)` with `x = 1 − 2r₊/r`, expand `f` in
Chebyshev polynomials, collocate at interior Chebyshev–Gauss points, and find
`ω` where the matrix is singular.

**Result: fails structurally.** A scan of the smallest singular value over the
`ω` plane shows it decreasing *monotonically* as `Im ω` becomes more negative,
with **no isolated minimum anywhere**. The matrix is numerically singular
throughout the lower half plane.

**Diagnosis (this is the useful part).** After peeling the outgoing factor
`e^{iΩr}`, the *ingoing* solution behaves as `e^{−2iΩr}`. For a quasinormal mode
`Im ω < 0`, so `e^{−2iΩr}` **decays** at infinity rather than diverging.
Requiring the peeled function to be a polynomial — i.e. analytic and bounded —
therefore does **not** exclude the ingoing solution: both solutions are
representable, the boundary condition is never actually imposed, and the
discretized operator has a genuine two-dimensional near-null space at every `ω`.

This is not a resolution, conditioning, or seeding problem, and it cannot be
fixed by tuning. Any method that selects the QNM by "peel the outgoing factor
and demand boundedness" is wrong for `Im ω < 0`. The correct selections are
(a) Leaver's minimal-solution condition on the Frobenius coefficients, or
(b) hyperboloidal slicing, which makes the ingoing solution genuinely singular
at the boundary rather than merely small.

**Not revisited in this form.** The module is retained, clearly marked
non-functional, because the diagnosis is reusable and because its coefficient
assembly is correct and shared with the Leaver construction.

*Aside:* hyperboloidal slicing is the standard fix in the pseudospectrum
literature, but it is not obviously available here — massive fields do not
propagate along null rays, so compactification at null infinity is not the
natural construction for `μ ≠ 0`. This should be checked before assuming
route (b) is open.

## 4. Fully symbolic derivation of the Leaver recurrence (session 2)

`sympy` with all of `a, b, M, ω, μ, m₁, m₂, Λ, Ω, σ, r₊, r₋` symbolic did not
terminate in 10 minutes when transforming the radial ODE to the Leaver variable
and clearing denominators. Same root cause as failure 1: symbolic simplification
of large rational-trigonometric expressions.

**Replaced by:** numerical recovery of the polynomial coefficients by FFT on a
circle (`mp5d.radial.leaver.polynomial_coefficients`). The coefficient functions
are evaluated numerically, multiplied by a clearing factor, sampled on the unit
circle and transformed. This is exact for polynomials and **self-validating**:
if the clearing factor missed a pole the recovered high-order tail does not
vanish and the function raises. Recovered degrees (12, 11, 10) are stable under
changes of the degree bound, and the routine reproduces a known test polynomial
to 2e-15.

Two bugs were found and fixed inside this approach, both worth recording:

* sampling on a circle of radius `ρ = 0.6` requires dividing by `ρⁿ ≈ 1e-57` at
  `n ≈ 250`, which amplifies roundoff catastrophically. Fixed by sampling on the
  **unit** circle with a half-sample phase offset, which also guarantees no
  sample lands on the singular point `u = 1`.
* `numpy.fft.ifft` has the opposite sign convention to the one needed and
  returns the coefficients **index-reversed**; the correct transform is
  `fft(vals)/N`. The symptom was a recovered "tail" of exactly 1.0 — the tail was
  the head.

## 5. Raw Hill determinant from the multi-term Leaver recurrence (session 2, unresolved)

After stripping the common power of `u` that the clearing factor introduces
(without which the first three rows of the recurrence vanish identically and the
matrix is singular for *every* `ω` — observed as a uniform `rel σ_min ≈ 1e-32`),
the Hill matrix has correct structure: no zero rows, `A` valuation 1 degree 9,
`B` 0/8, `C` 0/7.

It still shows **no isolated root**. The smallest singular values
(4.6e-5, 9.5e-6, 9.9e-7) form a decreasing sequence tied to the truncation edge:
the last columns have visibly reduced norm (0.62, 0.22 versus 1.0) because
`a_{K−1}`, `a_{K−2}` appear in fewer equations than interior coefficients. These
edge modes are spurious and mask the genuine singularity condition.

**Status: unresolved, and this is the blocking item.** The fix is the standard
one and is not conceptually hard: reduce the multi-term recurrence to three
terms by Leaver's Gaussian elimination, then impose the *minimal-solution*
condition through the continued fraction with its correct large-`n` tail
estimate, rather than reading a determinant off a hard truncation. That closes
the truncation edge properly. It was not reached before the session's resource
ceiling.

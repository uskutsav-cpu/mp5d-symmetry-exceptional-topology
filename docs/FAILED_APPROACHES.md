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

## 5b. CORRECTION to entry 5 (session 3)

Entry 5 attributed the Hill-truncation failure to an **edge-localized null
vector**. That diagnosis is **wrong** and is corrected here rather than edited
away.

Measured: the null vector's weight in the last five coefficients is `0.19`, not
dominant. The reproduction test
`tests/unit/test_radial_known_failures.py::test_hill_svd_condition_is_omega_insensitive`
now asserts `< 0.5` so the incorrect story cannot be re-asserted.

**Actual cause: dynamic range.** Recurrence rows carry an `n(n-1)` factor, so
Hill-matrix entries span many orders of magnitude. After row equilibration the
smallest singular value is tiny *and nearly independent of* `omega`, which
destroys the root condition. Taking an SVD was the mistake, not the truncation.

**Fix (session 3, working).** Do not form an SVD at all. The truncated
determinant is proportional to the forward-generated coefficient `a_N`, which
carries the same zeros with none of the dynamic range. Evaluating `a_N` and
accelerating the root sequence with Wynn's epsilon algorithm gives Solver B,
which now reproduces Solver A to `1e-12`.

## 6. Muller iterates walking onto the branch cut (session 3, mitigated)

`Om = sqrt(omega^2 - mu^2)` has a branch cut at `Re omega = 0`. Root searches
seeded loosely repeatedly converged to spurious near-imaginary "roots" there
(recorded with their values in `results/rejected_radial_roots.json`).

Mitigations that worked: clamp Muller iterates to `Re omega >= 0.05`, and start
the depth schedule coarse (depth 100), which widens the basin of attraction.

Mitigation that did **not** work, recorded so it is not retried: widening the
initial Muller triangle from `h = 1e-4|x0|` to `2e-2|x0|`. It made matters worse,
turning one failing seed into three. The narrow triangle is correct here.

## 7. Fully symbolic Leaver derivation, second attempt (session 3)

Re-confirmed that a fully symbolic transformation of the radial ODE to the
Leaver variable with all parameters symbolic does not terminate in 10 minutes.
The FFT-based numerical extraction is the route that works. Do not retry the
symbolic path without first fixing the parameters to numbers.

## 8. Reading the large-n ratio from the reduced recurrence (session 5, partial)

Attempting to extract `R_n = a_{n+1}/a_n` at `n` of order `10^3` from the
Gaussian-reduced three-term recurrence gives a **spurious, problem-independent**
answer: `R_n - 1 = 2/n` exactly, identical for two completely unrelated
parameter sets (`a=b=mu=0` and `a=0.2, b=0.3, mu=0.1`). That coincidence is the
giveaway that it is an artifact, not physics.

**Cause.** The reduction is a long recursive elimination; cancellation
accumulates and past `n ~ 500-800` the reduced coefficients lose their leading
structure. Diagnostic: `(alpha_n + beta_n + gamma_n)/alpha_n` should decay
smoothly like `1/n`, but is measured non-monotonic at the `1e-3`-`1e-5` level.

**Consequence.** The clean window for reading the asymptotics is `n = 200-400`,
which is enough to fix the leading coefficient `u1` but not enough to extract
`u2` reliably. This is the current obstacle to a higher-order tail.

**Not fixed.** Candidate fixes, untested: raise working precision inside the
reduction independently of the solve precision; or restructure the elimination
to avoid the long recursive chain.

Note this does *not* invalidate the CF values themselves -- deep terms
contribute little to the backward evaluation, and the benchmark agreements are
unaffected.

## 9. Trusting a small |dF/domega| outside the validated solver domain (session 8)

The augmented EP Newton solve produced |F| = 1.5e-11 and |dF/domega| = 4.3e-6 at
r2 = 0.44, extremality 0.0089 -- the signature of a near-double root. It was
correctly held back as unverified, and has now been **refuted**.

Cause: at r2 = 0.44 the Leaver series radius degenerates (the spurious
singularity x(-r_+) = 2 r_+/(r_+ + r_-) approaches the unit circle), so Solver A
returns depth-dependent roots that drift 2e-4 to 6e-4 per depth doubling while
its CF residual stays at 1e-13. Newton on (F, dF/domega) then locates a spurious
stationary point in that drift.

**Lesson, now enforced by the record:** a small residual does NOT imply a
converged root. Before treating any repeated-root diagnostic as physical, check
that the frequency itself is converged to better than the diagnostic. Here the
diagnostic (4.3e-6) was three orders BELOW the three-solver frequency spread
(3.5e-3), so it could not have carried information.

Ruled out as causes: horizon-exponent blow-up (|sigma| = 1.69, only ~2x
ordinary) and Solver C resolution (stable under contour-length and resolution
refinement).

## 10. Near-horizon multidomain resolution as the fix for r2 ~ 0.44 (session 9)

**Hypothesis:** the near-extremal failure was caused by the inner-horizon
singularity crowding the domain endpoint (Bernstein ellipse collapse), so a
scaled near-horizon coordinate y = (r-r_+)/(r_+-r_-) plus domain decomposition
should fix it.

**Built and validated:** Solver D reproduces ST5D l=0 (5.8e-6), ST5D l=1
(1.4e-6) and HH Table III r1 (2.9e-5).

**Hypothesis REFUTED.** At r2 = 0.44 Solver D still does not converge: the root
wanders in a ball of radius ~2e-3 as n_inner goes 60 -> 280 and n_outer
160 -> 400, with movement NOT decreasing. Resolving the near-horizon scale was
necessary but not sufficient.

**Second hypothesis also refuted:** double-precision conditioning. Measured
cond(M) at r2=0.44 is 5.3e11, comparable to the ST5D case at 3.2e11 which
converges to 1e-9. The hard case is not worse conditioned.

**Still unidentified.** The failing point sits at mu = 1.9, near where the l=2
branch transitions toward omega^2 = mu^2 (observed near mu ~ 2.0). The
difficulty may be spectral (branch point / mode accumulation) rather than
discretization, in which case no spatial resolution helps. UNTESTED.

Value retained: Solver D is a second recurrence-free method, and it agrees with
Solver C (2.1423) against the recurrence family (A 2.1395, B 2.1407),
strengthening the finding that the recurrence is the outlier near extremality.

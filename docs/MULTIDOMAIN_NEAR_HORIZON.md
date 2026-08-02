# Solver D — multidomain near-horizon solver

**Status: built and validated on ordinary modes. Does NOT resolve the
near-extremal region.** The hypothesis that motivated it is refuted.

## Construction

Exterior complex scaling with the near-horizon scale resolved explicitly:

* **inner domain** `y = (r − r₊)/(r₊ − r₋) ∈ [0, Y_m]`, kept real. In this
  scaled coordinate the inner-horizon singularity sits at `y = −1`, a distance
  of order unity from the domain *regardless of how small* `r₊ − r₋` becomes.
* **outer domain** `r = r_m + t e^{iθ}`, complex-scaled so the outgoing solution
  decays (`tan θ > −Ω_I/Ω_R`). Scaling begins only after the real near-horizon
  region — exterior complex scaling proper.
* **matching** at `r_m = r₊ + Y_m(r₊−r₋)`: continuity of value and of `dR/dr`,
  both expressed in `r` so the two parametrizations compare on equal footing.

Recurrence-free: no Frobenius recurrence, no Gaussian reduction, no continued
fraction, no Hill determinant, no Wynn acceleration. `Λ` recomputed at every `ω`.

## Validation on ordinary modes — passes

| case | Solver D | reference | difference |
| --- | --- | --- | --- |
| ST5D `l=0, n=0` | `0.53383120 − 0.38337163i` | published | `5.8e-6` |
| ST5D `l=1` | `1.01601617 − 0.36232687i` | published | `1.4e-6` |
| HH Table III r1 | `1.68110923 − 0.34722673i` | their CFM | `2.9e-5` |

## At r₂ = 0.44 — does not converge

Parameter invariance (`Y_m` 6–12, `n_inner` 90–140, `n_outer` 200–300, `θ`
50–70°, `L` 60–90) gives a spread of `1.2e-3`. Worse, a direct resolution study
shows **no convergence**:

| `n_inner` | ω | movement |
| --- | --- | --- |
| 60 | `2.142323 − 0.493082i` | — |
| 90 | `2.142097 − 0.492944i` | `2.7e-4` |
| 140 | `2.141930 − 0.490136i` | `2.8e-3` |
| 200 | `2.141460 − 0.491443i` | `1.4e-3` |
| 280 | `2.141963 − 0.491898i` | `6.8e-4` |

`n_outer` behaves the same way (`~2e-3`, not decreasing). The value wanders
inside a ball of radius `~2e-3` instead of settling.

## What this rules in and out

**Rules in:** Solver D agrees with Solver C (`2.1423`) and separates from the
recurrence family (A `2.1395`, B `2.1407`). Two recurrence-free methods now sit
together, which supports the earlier finding that the recurrence is the outlier
in this region.

**Rules out — the motivating hypothesis.** The near-horizon scale was *not* the
obstruction. Solver D resolves `r₊ − r₋ = 0.094` by construction and still fails
to converge, so the diagnosis in `docs/NEAR_EXTREMAL_CANDIDATE.md` (Bernstein
ellipse collapse) was at best incomplete.

**Rules out — conditioning.** Measured `cond(M)`:

| case | cond | converges? |
| --- | --- | --- |
| ST5D `l=0` | `3.2e11` | yes, to `1e-9` |
| two-spin `r₂=0.064` | `8.5e10` | yes |
| `r₂=0.19` | `7.1e9` | yes |
| **`r₂=0.44`** | **`5.3e11`** | **no** |

The hard case is no worse conditioned than the easy one. Double-precision
conditioning is therefore **not** the explanation either. This hypothesis is
recorded as refuted rather than quietly dropped.

## Remaining possibility, untested

The failing point sits at `μ = 1.9`, close to where the `l = 2` branch was
observed to transition toward `ω² = μ²` near `μ ≈ 2.0`. The mode may be near a
genuine branch point or in a region of dense mode accumulation, in which case no
amount of spatial resolution helps and the difficulty is spectral, not
discretization. **This has not been tested.** Distinguishing it needs the
remaining fallbacks — Wronskian matching, hyperboloidal slicing, or near-horizon
matched asymptotics — none of which were built.

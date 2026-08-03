# Primary result: a bounded multi-sector exclusion of exceptional points

**Outcome C.** No exceptional point exists in the searched domain, and the
exclusion is quantitative rather than "nothing showed up".

## What was searched

| | |
| --- | --- |
| atlas points | **50 558** |
| independent sectors | **7**: `(0,0) (1,1) (−1,−1) (1,0) (2,0) (2,1) (2,2)` |
| angular branches | 3 per sector: `l = l_min, l_min+2, l_min+4` |
| overtones | `N = 0,1,2,3` |
| parameter points with ≥2 branches | 4 459 |
| `s = (a+b)/2` | `0 → 0.42` |
| `δ = (a−b)/2` | `−0.20 → +0.20` |
| `μ` | `0 → 1.8` |
| inner horizon `r₋` | `0 → 0.2287` (`z₋ = r₋² → 0.0523`) |
| extremality `M−(|a|+|b|)²` | `≥ 0.294` |
| max continued-fraction residual | `3.6e-9` |
| continuation failures | 118 predictor-miss + 33 non-converged, all recorded |

The domain is independently validated: Solver A (continued fraction) versus
Solver C (exterior complex scaling, structurally recurrence-free) agree to
**max `|A−C| = 2.06e-9`** on 9 spanning points up to `r₋ = 0.2287` — beyond the
previously validated `r₋ ≤ 0.19`.

## The bounds

Two scale-invariant quantities, both bounded away from zero:

| quantity | minimum over the domain |
| --- | --- |
| pairwise branch gap `min\|ω_i − ω_j\|` (same sector, same parameters) | **0.5116** |
| scale-free root separation `\|a₁/a₂\|` (Cauchy, sees *unmodelled* partners too) | **0.3107** |
| discriminant `min\|D\| = min\|ω₊−ω₋\|²` across fitted effective models | **0.2618** |
| eigenvalue condition number `κ` at the strongest interaction | **8.73e3** (finite) |

An EP2 requires the branch gap, the root separation and `D` to vanish
*simultaneously*, and `κ` to diverge. None approaches zero: the minimum gap is
`≈ 0.51` where `|ω| ≈ 2–4`, i.e. a relative separation of order 15 %.

### Why these bounds and not `|dF/dω|`

The previous exclusion in this repository used `min |dF/dω| = 1.02e-2`. That
reasoning is **invalid** and has been withdrawn. The spectral condition `F` may
be multiplied by any nonvanishing analytic function without changing its zero
set, so `|dF/dω|` at a root can be made arbitrarily small by rescaling and
bounds nothing. A test exhibits two spectral conditions with identical zeros
whose `|dF/dω|` differ by `>1e5` while `|a₁/a₂|` agrees to 2 %.

Every bound above is invariant under such rescaling: they are built from
eigenvalue *differences*, not from the value of any particular condition
function.

## Verification of the tightest candidates

The ten smallest-gap candidates were put through the gate in
`docs/EP_VERIFICATION_GATE.md`. All ten are **REJECTED**, with the same
signature:

| gate | result |
| --- | --- |
| G2 argument-principle zero count = 2 | **pass** |
| G3 both roots located without a guess | **pass** |
| G1 scale-free separation → 0 | **fail** |
| G4 monodromy transposes the branches | **fail** |
| G6 splitting `∼ t^{1/2}` | **fail** |

G2 and G3 are preconditions — they say only that two distinct roots sit near
each other. The EP-specific evidence (G1, G4, G6) fails uniformly. The
classification is therefore **avoided or ordinary crossing**, not "unresolved".

## Structural explanation, not just absence

Three independent facts explain why exceptional points are hard to make here,
and together they are stronger than the numerics alone:

1. **Cross-sector coalescence is impossible.** Jordan chains of a direct sum are
   the union of the summands' chains, and `(m₁,m₂)` labels permanently decoupled
   blocks because `∂_φ, ∂_ψ` are Killing for arbitrary `a,b`. A coincidence of
   frequencies from different sectors is always semisimple
   (`docs/CROSS_SECTOR_EP_NO_GO.md`, unconditional form). This removes the
   entire `U(2)` multiplet channel — the most natural candidate — by symmetry.

2. **Within a sector, the surviving interactions are between different `(l,N)`,
   and they repel weakly.** The measured gaps decrease slowly and smoothly
   toward large `s` and large `μ` and never approach zero. The condition number
   at the strongest interaction is `8.7e3` — elevated over the ordinary value
   `18.3`, i.e. a genuine avoided crossing, but finite.

3. **The one regime where the spectrum does something dramatic — the
   quasiresonant limit — is a branch point in `μ`, not a coalescence.** There the
   damping collapses by 8.6 orders of magnitude
   (`docs/LONG_LIVED_FINAL.md`) while the neighbouring overtone stays at
   `−Im ω ~ 0.1`, so the branches separate rather than merge.

## Honest limitations

* **The minimum sits on the boundary of the searched box** (`s = 0.42`,
  `δ = −0.20`, `μ = 1.8`), and the gap decreases toward that corner. A probe
  beyond the box shows the gap continuing to fall to `0.472` at `s = 0.45`
  (extremality 0.19) at `μ = 1.8`, but *rising* again at `μ = 3.0` and `μ = 3.6`
  (0.70–1.12). So the trend does not extrapolate to zero — but the region
  `s > 0.45` was not systematically searched and **no claim is made there**.
* **`μ > 1.8` is not covered** by the atlas, only by the targeted probe.
* **The quasiresonant regime is not searched for EPs.** There Solvers C and D
  are inapplicable by formulation (required scaling angle → 90°, see
  `docs/NEAR_EXTREMAL_SPECTRAL_CLASSIFICATION.md`), so no independent
  confirmation is available and the standard for a claim cannot be met.
* **Higher overtones (`N > 3`) and higher `l` were not searched.**
* **This is not a theorem.** It is a bounded numerical exclusion over a declared
  domain, supported by a symmetry argument that rules out one specific channel
  exactly. The two must not be conflated.

## Machine-readable record

`results/negative_regions.json`, `results/all_collision_candidates.json`,
`results/all_rejected_candidates.json`, `results/verified_exceptional_points.json`
(empty, `n_verified = 0`), `results/effective_models.json`,
`results/pseudospectral_results.json`, `results/final_multisolver_validation.json`,
`results/full_branch_atlas.json` (hashed manifest; the 24 MB atlas itself is
regenerated by `scripts/build_branch_atlas.py`).

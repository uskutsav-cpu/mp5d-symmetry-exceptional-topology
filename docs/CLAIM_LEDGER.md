# Claim ledger

Statuses: *established in literature*, *minimally reproduced*, *numerically
observed*, *convergence verified*, *cross-solver verified*, *symmetry derived*,
*interval certified*, *continuum proved*, *conjectured*, *falsified*,
*unresolved*.

No claim may be called new / first / universal / symmetry-protected /
topological / certified until its evidence requirement is met.

| # | Claim | Status | Evidence |
| - | ----- | ------ | -------- |
| C1 | The massive KG equation on 5D MP separates in the stated coordinates with `√\|g\| = r Σ sinθ cosθ`. | established in literature; **independently re-derived** | `tests/unit/test_separation.py` (exact rational arithmetic) |
| C2 | Separated radial/angular potentials as written in `CONVENTIONS.md` §4. | minimally reproduced (re-derived from the metric, not transcribed) | `test_separated_potentials_match_implementation` |
| C3 | `(a,m₁) ↔ (b,m₂)` is an exact symmetry of the separated system. | symmetry derived | `test_symmetrized_radial_potential_is_exchange_invariant` |
| C4 | Spheroidicity is exactly `c² = 4sδ(ω²−μ²)`; it vanishes identically on `δ=0` **and** on `s=0`. | symmetry derived | `CONVENTIONS.md` §5; `test_degeneracy_is_lifted_only_by_the_product_s_delta` |
| C5 | At `δ=0`, `Â = l(l+2)` exactly for all `(n,m₁,m₂)` with `l = 2n+\|m₁\|+\|m₂\|`. | symmetry derived | `test_equal_spin_eigenvalue_is_exactly_l_l_plus_2` |
| C6 | At `δ=0` the full separated problem (radial + horizon condition) depends on `(m₁,m₂)` only through `m = m₁+m₂`. | symmetry derived | `test_equal_spin_u2_structure.py` |
| C7 | The equal-spin degenerate set at fixed `(l,m)` is a single irreducible SU(2) multiplet of dimension `l+1`; summing over `m` gives `(l+1)²`. | symmetry derived | `test_fixed_l_and_m_multiplet_is_a_single_su2_irrep` |
| C8 | **No-go:** the `U(2)` multiplet degeneracy cannot become defective at any `δ`, because its members carry distinct `(m₁,m₂)` and `(m₁,m₂)` labels permanently decoupled blocks. | symmetry derived | `docs/SYMMETRY_STRUCTURE.md` §3 |
| C9 | Any exceptional point in this system must occur within a single `(m₁,m₂)` sector, between different `(l,N)` branches. | symmetry derived (corollary of C8) | `docs/SYMMETRY_STRUCTURE.md` §4 |
| C10 | In diagonal sectors `m₁ = m₂`, all spectral quantities are exactly even in `δ` at fixed `s`. | symmetry derived | `docs/SYMMETRY_STRUCTURE.md` §4b |
| C11 | Angular spectral solver agrees with an independent finite-difference discretization. | cross-solver verified (angular only) | `test_spectral_matches_independent_finite_difference` |
| C12 | EPs in diagonal sectors occur in `±δ` pairs; an on-surface EP is codimension 2. | **conjectured** | — needs radial solver |
| C13 | Exceptional sets come in quadruples under `s→−s`, `δ→−δ`. | **conjectured, probably too strong** | radial operator depends on `s`, `δ` separately |
| C15 | `z = r² = 0` is an **ordinary point** of the radial equation whenever `ab ≠ 0`: the apparent `1/z` singularities cancel exactly because `W(0) = abG` and `P(0) = a²b²`. | symmetry derived | `test_origin_is_an_ordinary_point` |
| C16 | The radial problem has exactly two regular singular points (`z₊`, `z₋`) plus an irregular point at infinity — confluent Heun type. | symmetry derived | `docs/../src/mp5d/radial/structure.py`; C15 + indicial tests |
| C17 | Horizon exponent is exactly `σ₊ = (ω − m₁Ω_a − m₂Ω_b)/(2κ)`; ingoing branch is `(z−z₊)^{−iσ₊}`. | symmetry derived | `test_horizon_indicial_exponent_derivation`, `test_W_at_horizon_factorizes_through_angular_velocities` |
| C18 | The radial potential is exactly even in `r`, so it has **no** `1/r` tail and the asymptotic exponent is exactly `−3/2`: no Coulomb phase, unlike 4D massive Kerr. Follows from `D−3=2`. | symmetry derived | `test_potential_is_even_in_r_so_no_coulomb_tail` |
| C14 | First trustworthy QNM: 5D Schwarzschild-Tangherlini massless scalar `l=0`, `n=0`, `ω = 0.533835574268 − 0.383375368513i` (`r₊=1`). | **cross-solver verified** | matches Matyjasek arXiv:2107.04815 Table I to all 12 digits (`4.98e-13`); Solvers A and B agree to `2.15e-12`; `docs/FIRST_QNM_VALIDATION.md` |
| C19 | Four Schwarzschild-Tangherlini benchmarks (`l=0,1` × `n=0,1`) reproduced to `1e-13`–`1e-11`. | cross-solver verified | `results/first_qnm.json`, `tests/unit/test_first_qnm.py` |
| C20 | The raw Frobenius recurrence has width **9 for the static case a=b=0** and **13 for generic two spins** (corrected: an earlier entry said 9 universally); the general Gaussian reduction to three terms is exact. | numerically observed + validated on synthetic recurrences | `tests/unit/test_recurrence_synthetic.py` |
| C6/C7 numerical | The equal-spin U(2) multiplet degeneracy is **confirmed numerically**: `(1,1)`, `(2,0)`, `(0,2)` at `l=2, m=2`, `a=b=0.25` share one frequency to `4.67e-15`, and split by `5.38e-02` once `δ ≠ 0`. | **numerically observed**, confirming a symmetry-derived prediction | `tests/unit/test_first_qnm.py` |
| C21 | Exchange symmetry `(a,m₁) ↔ (b,m₂)` holds numerically on two-spin backgrounds to `4.15e-15`. | cross-solver verified | `test_exchange_symmetry_is_exact_for_two_spins` |
| C22 | The recurrence's characteristic equation is `A(1/R)=0`; `x=1` is a multiple root (`A` mult ≥5, `B`,`C` mult 4), so `R_n = 1 + u₁n^{-1/2} + …` with half-integer powers. | symbolically derived + numerically confirmed | `docs/RECURRENCE_ASYMPTOTICS.md`; exponent test `n^{1/2}` converges / `n^1` diverges |
| C23 | Leading tail coefficient `u₁ = −√(−2c)`, `c = iΩ(r₊−r₋)`, sign fixed by minimality. | numerically observed, consistent with derivation | measured `(R_n−1)√n` → `−0.523+1.020i` |
| C24 | The leading-order tail measurably reduces required depth (~1.35×, +1.4 digits at depth 200). | numerically observed | `results/tail_validation.json` |
| C25 | Higher tail orders (`u₂`+). | **not derived** | blocked by reduction degradation past `n~500–800` |
| C26 | An independent non-recurrence solver (Solver C). | **not built** | — |
| C27 | Solver C (exterior complex scaling) is recurrence-free and reproduces the static fundamental, static overtone, `l=1`, and all three Huang–Huang Table III two-spin points. | **independently verified** | `docs/SOLVER_C.md`; independence enforced by test |
| C28 | Physical QNMs are invariant under scaling angle (45–75°), contour length (40–130) and resolution (150–380); continuum artifacts are not. | resolution verified | `tests/unit/test_solver_c.py` |
| C29 | Three solvers agree below `1e-5` on 10/10 benchmark points, including `r₂ = 0.140` and `0.190` — beyond the `r₂ ≳ 0.1` limit Huang–Huang state for their own CFM. | **independently verified** | `results/three_solver_validation.json` |
| C30 | Any exceptional point in this system. | **NOT SEARCHED** | no atlas, no EP work |
| C31 | EP2 augmented conditions implemented with an exact total `dF/dω` (Λ recomputed every evaluation, never frozen). | numerically observed | `src/mp5d/exceptional/ep_solver.py` |
| C32 | ~~No EP2 in sector (1,1) … min `|dF/dω| = 1.02e-2`~~ | **WITHDRAWN — invalid reasoning** | `|dF/dω|` is not invariant under rescaling of the spectral condition, so it bounds nothing. Superseded by C37. See `test_root_separation_is_invariant_under_rescaling`. |
| C33 | A near-double root at `r₂=0.44`, extremality `0.0089`. | **REJECTED — outside the verified domain**, not evidence either way | `results/rejected_collisions.json` |
| C34 | The near-extremal candidate (`r₂=0.44`, extremality `0.0089`) is a **numerical artifact** of Solver A's non-convergent drift, not an EP. Solver C shows contrast 0.13 (no root); Solver A drifts 2–6×10⁻⁴ per depth doubling; three-solver spread `3.5×10⁻³` vastly exceeds the claimed `|dF/dω| = 4.3×10⁻⁶`. | **falsified** | `docs/NEAR_EXTREMAL_CANDIDATE.md` |
| C35 | The region `r₂ ≈ 0.44` is outside the validated domain of **all three** solvers; best achievable agreement is `3.5×10⁻³`. | numerically observed | same |
| C36 | No EP exists near extremality. | **NOT ESTABLISHED** — region unresolvable at present accuracy | — |
| C37 | **No EP2 in the searched domain**: 7 sectors, 3 `l` each, `N=0..3`, `s ≤ 0.42`, `|δ| ≤ 0.20`, `μ ≤ 1.8`, `r₋ ≤ 0.2287`, extremality `≥ 0.294`, 50 558 atlas points. Bounds: branch gap `≥ 0.5116`, scale-free root separation `≥ 0.3107`, `min|D| ≥ 0.2618`, `κ ≤ 8.73e3` (finite). | **cross-solver verified (bounded exclusion)** | `docs/FINAL_NEGATIVE_RESULT.md`; A-vs-C max difference `2.06e-9` on 9 spanning points |
| C38 | The ten tightest candidates are **avoided or ordinary crossings**, not EPs: G2/G3 pass, G1/G4/G6 fail uniformly. | **cross-solver verified** | `results/all_rejected_candidates.json` |
| C39 | Full discrete symmetry group is Klein four `{1,E,P,EP}` on (parameters × sector labels). Diagonal sectors even in `δ`; anti-diagonal even in `s`; `(0,0)` both. | symmetry derived | `docs/SYMMETRY_GROUP.md`; verified `<1e-12` |
| C40 | C13 as originally stated (**quadruples under independent `s→−s`, `δ→−δ`**) is **false in general**, true exactly in sector `(0,0)`. | **falsified / corrected** | negative control in `test_symmetry_group.py` |
| C41 | Equal-spin multiplet dimension is exactly `l+1` for every allowed `m`, and `(l+1)²` per level, **for all `l`** (O7 discharged). | symmetry derived (explicit count) | verified exactly to `l = 40` |
| C42 | Cross-sector no-go now holds **unconditionally**: Jordan chains of a direct sum are the union of the summands' chains, so cross-sector coincidence cannot create defectiveness whatever the within-block multiplicity. O6.5–O6.7 discharged. | symmetry derived | `docs/CROSS_SECTOR_EP_NO_GO.md` |
| C43 | Diagonal-sector exceptional lines, **if any exist**, must come in mirror pairs and terminate quadratically on the equal-spin surface (`δ = ±√(cτ)`). Recorded before the search. | symmetry derived (conditional) | `docs/SYMMETRY_GROUP.md` §3 |
| C44 | The near-extremal solver breakdown is **not caused by extremality**. Ordinary modes continue to extremality `0.0103` with residual `1e-13`. The controlling variable is `μ`: as `Re Ω → 0` the required complex-scaling angle → 90°, so Solvers C/D are inapplicable *by formulation*. | **numerically observed + derived**; corrects the earlier framing | `docs/NEAR_EXTREMAL_SPECTRAL_CLASSIFICATION.md` |
| C45 | Massive-scalar branches reach damping `8.65e-10` (8.6 orders below the massless value) in the quasiresonant limit; the transition is a **branch point in `μ`**, not a coalescence. | numerically observed | `docs/LONG_LIVED_FINAL.md` |
| C46 | Effective-model out-of-sample error in `t=δ²` is `1.7e-5` for diagonal sector `(1,1)` versus `1.4e-1` for off-diagonal `(2,0)` — four orders apart, independently confirming `δ`-evenness. | numerically observed | `results/effective_models.json` |
| C47 | Prior work's `r₂` is the inner horizon **radius** `r₋`; the code exposes `z₋ = r₋²`. A guard `z_minus ≤ 0.20` is `r₋ ≤ 0.447`, not `r₋ ≤ 0.20`. | documentation correction | `docs/NEAR_EXTREMAL_SPECTRAL_CLASSIFICATION.md` |

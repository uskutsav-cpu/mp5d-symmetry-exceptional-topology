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

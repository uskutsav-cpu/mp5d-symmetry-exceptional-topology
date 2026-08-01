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
| C14 | Any QNM frequency of this system. | **not yet computed** | — |

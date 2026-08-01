# Proof obligations

## Discharged

* **O1.** Separation of the massive KG equation in the stated conventions —
  discharged exactly and symbolically (`tests/unit/test_separation.py`).
* **O2.** Exchange invariance under `(a,m₁) ↔ (b,m₂)` — discharged exactly.
* **O3.** `c² = 4sδ(ω²−μ²)` — discharged exactly.
* **O4.** Equal-spin dependence on `m₁+m₂` only — discharged exactly.
* **O5.** `(l+1)`-dimensional single-irrep structure of the equal-spin multiplet
  and `(l+1)²` total — discharged (exact counting, `l ≤ 6`; the general
  statement follows from the SO(4) → U(2) branching and should be written out).

## Outstanding

* **O6.** The no-go (claim C8) is currently argued at the level of the
  *separated* problem. To be airtight it needs a statement that the separated
  block decomposition is complete, i.e. that the `(m₁,m₂)` Fourier
  decomposition of the scalar field on this background is a direct-sum
  decomposition of the relevant function space, with no cross-block resolvent
  coupling. This is standard but should be written down rather than assumed.
* **O7.** The `l ≤ 6` multiplet counting should be replaced by a proof for all
  `l`.
* **O8.** Every certification claim must separate
  `E_total ≤ E_discretization + E_truncation + E_arithmetic + E_root`, and must
  not promote a finite-dimensional certificate into a continuum theorem.
* **O9.** Claims C12, C13 are conjectures and carry no evidence yet.

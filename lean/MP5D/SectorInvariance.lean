import Mathlib
namespace MP5D
section
variable {K V : Type*} [Field K] [AddCommGroup V] [Module K V]

/-- The fixed space of a commuting projection is invariant. Idempotence is
not needed for this implication; it is supplied when interpreting P as a projector. -/
theorem projector_fixed_space_invariant (A P : V →ₗ[K] V)
    (hcomm : ∀ x, A (P x) = P (A x)) (v : V) (hv : P v = v) :
    P (A v) = A v := by
  rw [← hcomm v, hv]

/-- Commutation with both symmetry actions preserves a simultaneous eigensector. -/
theorem two_symmetry_sector_invariant (A S T : V →ₗ[K] V)
    (hS : ∀ x, A (S x) = S (A x)) (hT : ∀ x, A (T x) = T (A x))
    (z q : K) (v : V) (hvS : S v = z • v) (hvT : T v = q • v) :
    S (A v) = z • A v ∧ T (A v) = q • A v := by
  constructor
  · rw [← hS v, hvS, map_smul]
  · rw [← hT v, hvT, map_smul]
end
end MP5D

import Mathlib
namespace MP5D

/-- Exchange covariance gives SET parity in a diagonal sector. No labelling
or branch uniqueness is silently added. -/
theorem diagonal_spectral_set_parity {W : Type*}
    (Spec : ℝ → ℤ → ℤ → Set W)
    (hexchange : ∀ d m n, Spec (-d) n m = Spec d m n) (d : ℝ) (m : ℤ) :
    Spec (-d) m m = Spec d m m := hexchange d m m

/-- Branch parity follows only after uniqueness in the common spectral set. -/
theorem unique_branch_parity {W : Type*} (Spec : ℝ → Set W) (branch : ℝ → W)
    (hset : ∀ d, Spec (-d) = Spec d)
    (hmem : ∀ d, branch d ∈ Spec d)
    (hunique : ∀ d, ∀ x ∈ Spec d, x = branch d) (d : ℝ) :
    branch (-d) = branch d := by
  apply hunique d
  rw [← hset d]
  exact hmem (-d)

/-- An even differentiable real function has zero derivative at the origin.
This is not a geometric theorem about exceptional-line termination. -/
theorem even_first_derivative_zero (f : ℝ → ℝ) (d : ℝ)
    (heven : ∀ x, f (-x) = f x) (hf : HasDerivAt f d 0) : d = 0 := by
  have hneg := hf.comp 0 ((hasDerivAt_id (0 : ℝ)).neg)
  have heq : (fun x => f (-x)) = f := funext heven
  simp only [neg_zero] at hneg
  rw [heq] at hneg
  have h := hf.unique hneg
  linarith
end MP5D

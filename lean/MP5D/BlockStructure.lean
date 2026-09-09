import Mathlib

/-! Algebraic product-space statements only. No unbounded operator domain,
resolvent continuation, continuum radial equation, or QNM existence is assumed
or established here. Kernel compilation must succeed before citing this file. -/
namespace MP5D

def iterateOp {V : Type*} (A : V → V) : ℕ → V → V
  | 0, v => v
  | n + 1, v => A (iterateOp A n v)

def block {V W : Type*} (A : V → V) (B : W → W) (v : V × W) : V × W :=
  (A v.1, B v.2)

theorem block_power {V W : Type*} (A : V → V) (B : W → W)
    (n : ℕ) (v : V) (w : W) :
    iterateOp (block A B) n (v, w) = (iterateOp A n v, iterateOp B n w) := by
  induction n with
  | zero => rfl
  | succ n ih => simp [iterateOp, ih, block]

section Linear
variable {K V W : Type*} [Field K]
variable [AddCommGroup V] [Module K V] [AddCommGroup W] [Module K W]

def shifted (A : V →ₗ[K] V) (z : K) (v : V) : V := A v - z • v

/-- A generalized-kernel condition splits into exactly its block conditions.
This is stated for algebraic linear maps, not the paper's continuum QNM pencil. -/
theorem generalized_kernel_product (A : V →ₗ[K] V) (B : W →ₗ[K] W)
    (z : K) (n : ℕ) (v : V) (w : W) :
    iterateOp (block (shifted A z) (shifted B z)) n (v, w) = (0, 0) ↔
      iterateOp (shifted A z) n v = 0 ∧ iterateOp (shifted B z) n w = 0 := by
  rw [block_power]
  simp only [Prod.mk.injEq]

/-- Each link in a generalized chain projects to a link in each block. -/
theorem jordan_link_projects (A : V →ₗ[K] V) (B : W →ₗ[K] W)
    (z : K) (v v' : V) (w w' : W)
    (h : block (shifted A z) (shifted B z) (v', w') = (v, w)) :
    shifted A z v' = v ∧ shifted B z w' = w := by
  exact Prod.mk.inj h

/-- If the two blocks have no length-two generalized-kernel enlargement, neither
can their algebraic direct product. Existing defects inside a block are not denied. -/
theorem no_chain_created_between_blocks (A : V →ₗ[K] V) (B : W →ₗ[K] W)
    (z : K)
    (hA : ∀ v, iterateOp (shifted A z) 2 v = 0 → shifted A z v = 0)
    (hB : ∀ w, iterateOp (shifted B z) 2 w = 0 → shifted B z w = 0)
    (v : V) (w : W)
    (h : iterateOp (block (shifted A z) (shifted B z)) 2 (v, w) = (0, 0)) :
    block (shifted A z) (shifted B z) (v, w) = (0, 0) := by
  have hvw := (generalized_kernel_product A B z 2 v w).mp h
  simp only [block, hA v hvw.1, hB w hvw.2]
end Linear
end MP5D

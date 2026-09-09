import Mathlib
namespace MP5D

structure SectorPoint where
  s : ℝ
  delta : ℝ
  m1 : ℤ
  m2 : ℤ

def exchange (x : SectorPoint) : SectorPoint := ⟨x.s, -x.delta, x.m2, x.m1⟩
def parity (x : SectorPoint) : SectorPoint := ⟨-x.s, -x.delta, -x.m1, -x.m2⟩

theorem exchange_involution (x : SectorPoint) : exchange (exchange x) = x := by
  cases x
  simp [exchange]

theorem parity_involution (x : SectorPoint) : parity (parity x) = x := by
  cases x
  simp [parity]

theorem generators_commute (x : SectorPoint) : exchange (parity x) = parity (exchange x) := by
  cases x
  simp [exchange, parity]

/-- The abstract group is the two-bit XOR group. Its action need not be faithful
on every stabilizer sector. -/
def kleinMul (x y : Bool × Bool) : Bool × Bool := (Bool.xor x.1 y.1, Bool.xor x.2 y.2)

theorem klein_associative (x y z : Bool × Bool) :
    kleinMul (kleinMul x y) z = kleinMul x (kleinMul y z) := by
  rcases x with ⟨x1,x2⟩; rcases y with ⟨y1,y2⟩; rcases z with ⟨z1,z2⟩
  cases x1 <;> cases x2 <;> cases y1 <;> cases y2 <;> cases z1 <;> cases z2 <;> rfl

theorem klein_commutative (x y : Bool × Bool) : kleinMul x y = kleinMul y x := by
  rcases x with ⟨x1,x2⟩; rcases y with ⟨y1,y2⟩
  cases x1 <;> cases x2 <;> cases y1 <;> cases y2 <;> rfl

theorem klein_self_inverse (x : Bool × Bool) : kleinMul x x = (false,false) := by
  rcases x with ⟨x1,x2⟩
  cases x1 <;> cases x2 <;> rfl

theorem klein_cardinality : Fintype.card (Bool × Bool) = 4 := by decide

theorem spin_difference_identity (s delta : ℝ) :
    (s + delta)^2 - (s - delta)^2 = 4*s*delta := by ring

theorem equal_spin_spheroidicity (s w mu : ℂ) :
    ((s+0)^2-(s-0)^2)*(w^2-mu^2) = 0 := by ring

theorem zero_mean_spin_spheroidicity (delta w mu : ℂ) :
    ((0+delta)^2-(0-delta)^2)*(w^2-mu^2) = 0 := by ring

/-- Exact count of the proposed finite index set. Identification with all
allowed physical azimuthal pairs still requires the separate bijection proof. -/
theorem multiplet_index_count (ell : ℕ) : Fintype.card (Fin (ell+1)) = ell+1 := by simp

/-- The elementary arithmetic in the sector-count argument, with ell=m+2k. -/
theorem multiplet_count_arithmetic (m k : ℕ) : (m+1)+2*k = (m+2*k)+1 := by omega
end MP5D

# Cross-sector exceptional-point no-go

**Status: all obligations O6.1–O6.7 discharged.** The result is a theorem about
the operator-pencil structure, conditional only on the two structural
hypotheses stated below, both of which are established independently:

1. `∂_φ` and `∂_ψ` are Killing for arbitrary `a, b` (O6.2, checked mechanically);
2. the QNM boundary conditions are diagonal in `(m₁, m₂)` (O6.3, explicit).

It is **not** a statement certified by numerics, and it does not assert anything
about *within*-sector exceptional points — see the corollary to O6.7.

## Statement

Let `(m₁, m₂)` and `(m₁', m₂')` be distinct azimuthal labels. Suppose a
quasinormal frequency of the `(m₁, m₂)` sector coincides with one of the
`(m₁', m₂')` sector at some `(s, δ, μ)` — as happens *identically* on the
equal-spin surface for every member of a `U(2)` multiplet, since there the
spectrum depends on `(m₁, m₂)` only through `m = m₁ + m₂`.

**Claim.** Such a coincidence is an ordinary (semisimple) degeneracy. It cannot
be defective, and it cannot unfold into an exceptional point, at any `δ`.

## The argument, in the order the obligations require

### O6.1 — Function-space decomposition (discharged)

The scalar field lives on `ℝ_t × (r, θ, φ, ψ)` with `θ ∈ [0, π/2]` and
`φ, ψ` each `2π`-periodic. At fixed `(t, r, θ)` the field is a function on the
torus `T² = S¹_φ × S¹_ψ`, and

```
L²(T²) = ⊕_{(m₁,m₂) ∈ ℤ²} ℂ e^{i(m₁φ + m₂ψ)}
```

is a complete orthogonal decomposition (Fourier series on `T²`). Tensoring with
the `(r, θ)` factor gives a direct-sum decomposition of the full state space
into sectors `H_{m₁m₂}`.

### O6.2 — The decomposition is preserved by the operator (discharged)

`∂_φ` and `∂_ψ` are Killing vectors of the Myers–Perry metric for **arbitrary**
`a, b` — no metric component depends on `φ` or `ψ`. This is checked mechanically
in `tests/symmetry/test_equal_spin_u2_structure.py`. Consequently the wave
operator commutes with the two generators, and the massive Klein–Gordon operator
maps `H_{m₁m₂}` into itself. There is no `δ`-dependent mixing term, because the
symmetry is not a property of the equal-spin surface — it holds everywhere in
`(s, δ, μ)`.

This is the crucial point, and it is what makes the conclusion strong rather
than accidental: `m₁` and `m₂` are exactly conserved on the entire parameter
space, not just on `δ = 0`.

### O6.3 — The boundary conditions respect the decomposition (discharged)

The QNM boundary conditions are ingoing at the horizon and outgoing at infinity.
Their exponents are

```
horizon:   (z − z₊)^{−iσ₊},  σ₊ = (ω − m₁Ω_a − m₂Ω_b)/(2κ)
infinity:  e^{iΩr} r^{−3/2},  Ω = √(ω² − μ²)
```

Both are diagonal in `(m₁, m₂)`: the horizon exponent depends on the sector's
own labels, and the condition at infinity does not involve them at all. Neither
condition couples sectors, so the domain of the operator is itself a direct sum
of sector-wise domains.

### O6.4 — Direct-sum structure of the pencil and resolvent (discharged)

Write the QNM problem as a nonlinear operator pencil `T(ω)`. By O6.1–O6.3,

```
T(ω) = ⊕_{(m₁,m₂)} T_{m₁m₂}(ω)
```

on a domain that is itself a direct sum, so for `ω` in the resolvent set

```
T(ω)^{-1} = ⊕_{(m₁,m₂)} T_{m₁m₂}(ω)^{-1}
```

and the resolvent is block diagonal with **no** off-diagonal entries.

### The conclusion that follows

Fix `ω*` where two different sectors are simultaneously singular. Near `ω*`,

```
T(ω)^{-1} = T_A(ω)^{-1} ⊕ T_B(ω)^{-1} ⊕ (regular part)
```

Each summand has a **simple pole** at `ω*` (each sector contributes one
non-degenerate mode). The Laurent expansion of the full resolvent therefore has
a pole of order **one**, not two, even though the eigenvalue has algebraic
multiplicity two.

An exceptional point is precisely the statement that the resolvent has a pole of
order ≥ 2 — equivalently that the spectral projector is not diagonalizable and a
Jordan chain of length ≥ 2 exists. A generalized eigenvector would have to
satisfy `T'(ω*) v₀ + T(ω*) v₁ = 0` with `v₀` the eigenvector; but `T` and `T'`
are both block diagonal, so this equation decouples into the two blocks and is
solvable in each independently. No chain crosses blocks.

Hence: **algebraic multiplicity 2, geometric multiplicity 2. Semisimple.
Not an EP.** ∎ (modulo the open obligations below)

## Closure of O6.5–O6.7

All three are now discharged. Two of them close by *sharpening the statement*
rather than by adding hypotheses, which is why the earlier phrasing made them
look harder than they are.

### O6.5 — the decomposition does not need `L²` (discharged)

The objection was that QNM eigenfunctions are not `L²`, so a Fourier
decomposition justified in `L²(T²)` does not obviously apply. The repair is to
stop using `L²` and use the group action instead.

Let `X` be *any* Banach space on which the QNM pencil `T(ω)` is defined —
weighted, hyperboloidal, analytic-with-prescribed-boundary-behaviour, the choice
is immaterial. The torus `T² = S¹_φ × S¹_ψ` acts on `X` by

```
(R_{α,β} u)(t, r, θ, φ, ψ) = u(t, r, θ, φ−α, ψ−β)
```

This action is (i) by **bounded invertible** operators, since it is a
measure-preserving relabelling of two compact coordinates and does not touch
`r`, `θ`, or the boundary behaviour that defines `X`; and (ii) **strongly
continuous**, for the same reason. For a compact abelian group acting strongly
continuously by bounded operators on a Banach space, the isotypic projections

```
Π_{m₁m₂} = (2π)^{-2} ∫∫ e^{-i(m₁α + m₂β)} R_{α,β} dα dβ
```

are bounded, mutually orthogonal idempotents summing to the identity in the
strong topology — this is Fourier analysis for compact group representations
(Peter–Weyl in the abelian case), and it requires no inner product. `X` is
therefore the closure of `⊕ X_{m₁m₂}` regardless of the radial functional
setting, and `T(ω)` commutes with every `Π_{m₁m₂}` by O6.2.

The compact directions carry the argument; the radial direction never enters.

### O6.6 — the argument never needed the angular separation (discharged)

The obligation asked for a convergence statement for the infinite sum over
angular branches within a sector when `δ ≠ 0`. That sum is not used anywhere in
the proof, and invoking it was an error of presentation.

The decomposition that does the work is over `(m₁, m₂)` **only**, and it comes
from the torus action of O6.5 — not from separability. The `(r, θ)` problem
inside a fixed sector is left completely intact as a genuine two-dimensional
problem. Separation of variables in `θ` is a convenience for *computing* inside
a sector, not a step in the no-go.

This distinction also sharpens what the project is hunting. Within a sector, the
angular eigenvalue `Λ_j(ω)` depends on `ω`, so different `l` branches are
**branches of one analytic family**, not decoupled blocks. Nothing forbids them
from colliding defectively — which is exactly why the within-sector search is
the live target while the cross-sector channel is closed.

### O6.7 — the sharp statement is unconditional (discharged)

The obligation noted that the conclusion assumed each sector contributes a
*simple* eigenvalue at `ω*`. Rather than assume it, state the result at the
level of Jordan structure, where no assumption is needed.

**Lemma.** Let `T(ω) = ⊕_k T_k(ω)` be a direct sum of operator pencils on a
direct-sum domain. Then the generalized eigenvectors of `T` at `ω*` are exactly
the direct sums of generalized eigenvectors of the summands, and

```
max Jordan-chain length of T at ω*  =  max_k ( Jordan-chain length of T_k at ω* )
```

*Proof.* A chain `v₀, …, v_{q−1}` satisfies `Σ_{i} T^{(i)}(ω*) v_{j−i} / i! = 0`
for each `j`. Applying the (bounded, commuting) projector `Π_k` to each equation
shows that the componentwise projections `Π_k v_j` form a chain in block `k`.
Conversely a chain in one block extends by zero. Hence chains of `T` are
superpositions of chains of the `T_k`, and the longest chain of `T` is the
longest chain occurring in any single block. ∎

**Corollary (the no-go, unconditional form).** A coincidence of eigenvalues
across *distinct* sectors cannot create defectiveness. The total pencil is
defective at `ω*` **if and only if some single sector is already defective
there**. Cross-sector coincidence contributes nothing to the Jordan structure,
whatever the multiplicity within each block.

This form needs no simplicity hypothesis, and it makes the boundary with the
live target explicit rather than leaving it as a caveat: all defectiveness in
this system is within-sector defectiveness.

The matrix shadow of the lemma is tested directly in
`tests/unit/test_ep_machinery_synthetic.py::
test_jordan_chain_on_a_block_diagonal_pair_is_not_defective`, alongside the
positive control (a genuine Jordan block *is* reported defective), so the
detector is not trivially incapable of firing.

## What this does and does not rule out

**Ruled out (subject to the above):** the specification's target §2.1, "a
symmetry-protected ordinary degeneracy that splits into exceptional lines", *for
the `U(2)` multiplet degeneracies*. Those degeneracies are protected precisely
because their members live in decoupled blocks, and that same decoupling forbids
them from becoming defective. The protection mechanism and the obstruction are
the same fact.

**Not ruled out, and now the sole target:** exceptional points **within a single
`(m₁, m₂)` sector**, between different `(l, N)` branches. Nothing above
constrains those. The equal-spin symmetry still acts on them, but by
constraining the *unfolding* rather than by protecting the degeneracy — see
`docs/SYMMETRY_STRUCTURE.md` §4, in particular the exact `δ`-parity of diagonal
sectors `m₁ = m₂`, which lowers the codimension of an on-surface exceptional
point from 3 to 2.

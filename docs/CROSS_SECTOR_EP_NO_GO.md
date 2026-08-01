# Cross-sector exceptional-point no-go

**Status: provisional.** Obligations O6.1–O6.4 below are discharged;
O6.5–O6.7 are **open**. Until they close this is an argument, not a theorem, and
must not be described as one.

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

## Open obligations

* **O6.5 (open).** The completeness of the `(m₁, m₂)` decomposition is stated for
  `L²`. QNMs are *not* `L²` — the outgoing eigenfunctions diverge at infinity.
  The argument therefore needs the decomposition on the space where the QNM
  pencil is actually defined (a weighted or hyperboloidal space, or a space of
  analytic functions with the appropriate boundary behaviour). The Fourier
  decomposition in `φ, ψ` is expected to survive unchanged, because it concerns
  compact directions and is untouched by the radial functional-analytic setup,
  but this has not been written out.
* **O6.6 (open).** Written for the *separated* problem. The step from "each
  separated block has a simple eigenvalue" to "the unseparated operator pencil
  has a semisimple eigenvalue" uses the angular decomposition as well, which for
  `δ ≠ 0` is an infinite sum over angular branches within a sector rather than a
  finite one. Needs a convergence statement.
* **O6.7 (open).** Assumes each sector contributes a *simple* eigenvalue at `ω*`.
  If a single sector were itself defective there — a **within**-sector EP — the
  conclusion about that block does not apply. This is not a gap in the
  cross-sector claim, but it is exactly the case the project is hunting, and the
  two must not be conflated.

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

# The discrete symmetry group and what it forces on exceptional sets

Status: **symmetry derived**, every statement checked numerically to `< 1e-12`
in `tests/symmetry/test_symmetry_group.py`.

## 1. The group

Two generators act on the pair *(parameters, sector labels)*:

| | action on `(a, b)` | action on `(m₁, m₂)` | in `(s, δ)` |
| --- | --- | --- | --- |
| `E` (exchange) | `(a,b) → (b,a)` | `(m₁,m₂) → (m₂,m₁)` | `(s,δ) → (s,−δ)` |
| `P` (parity) | `(a,b) → (−a,−b)` | `(m₁,m₂) → (−m₁,−m₂)` | `(s,δ) → (−s,−δ)` |

`E` is the relabelling of the two rotation planes. `P` is the isometry
`φ → −φ`, `ψ → −ψ` combined with reversing both rotation senses. Both are
exact for arbitrary `a, b, μ`; `E² = P² = 1` and `EP = PE`, so

```
G = {1, E, P, EP} ≅ Z₂ × Z₂     (Klein four-group)
```

`EP` acts as `(s,δ) → (−s, δ)` with `(m₁,m₂) → (−m₂,−m₁)`.

## 2. The essential point: no generator acts on parameters alone

Every element of `G` moves the sector as well as the parameters. A symmetry of
the **parameter space at fixed sector** therefore exists only in sectors that
the group element *stabilizes*:

| stabilized sector | stabilizing element | consequence |
| --- | --- | --- |
| `m₁ = m₂` (diagonal) | `E` | `ω` is exactly **even in `δ`** at fixed `s` |
| `m₁ = −m₂` (anti-diagonal) | `EP` | `ω` is exactly **even in `s`** at fixed `δ` |
| `(0,0)` | all of `G` | even in `δ` **and** in `s` separately |

The anti-diagonal statement was not previously recorded. The `(0,0)` statement
settles claim **C13**.

### C13 as originally stated is false

C13 conjectured that exceptional sets come in quadruples under `s → −s` and
`δ → −δ` *independently*, and was flagged "probably too strong". It is: those
two reflections generate a parameter-space symmetry **only in the sector
`(0,0)`**. In a general sector neither reflection is a symmetry on its own —
each must be accompanied by a change of sector. The negative control
`test_delta_evenness_fails_in_a_non_diagonal_sector` confirms this explicitly
(sector `(2,0)` differs by `> 1e-3` under `δ → −δ`).

## 3. Consequence for exceptional sets in diagonal sectors

This is the structural prediction, and it is falsifiable.

In a diagonal sector `m₁ = m₂ = m`, every branch `ω_j(s, δ, μ)` is even in `δ`.
Assuming analyticity in the parameters, each branch is therefore an analytic
function of

```
t ≡ δ²
```

Let `D(s, t, μ)` be the discriminant of the local two-branch problem (the
product of squared branch differences). `D` inherits the dependence on `δ` only
through `t`. An EP2 is `D = 0`, which is **two real equations**.

* In `(s, t, μ)` space: two equations, three unknowns → the exceptional set is a
  **one-dimensional curve**.
* The physical `δ`-space is the double cover `δ = ±√t`, defined only for `t ≥ 0`.

Two structural consequences follow, neither of which requires knowing whether an
EP exists:

**(a) Exceptional lines come in mirror pairs.** A curve segment at `t > 0` lifts
to two curves at `δ = +√t` and `δ = −√t`, exchanged by `E`.

**(b) Exceptional lines can *terminate* on the equal-spin surface.** Where the
`(s,t,μ)` curve crosses `t = 0` transversally into `t < 0`, the real `δ`-locus
has no continuation: the two mirror branches meet at `δ = 0` and stop. Near such
a point, parametrizing the curve by arclength `τ` with `t ≈ c τ`,

```
δ = ±√(cτ)
```

so the exceptional line meets the equal-spin surface **quadratically** — the
pair of mirror lines forms a smooth parabola-like arc whose vertex lies exactly
on `δ = 0`, approaching it with infinite slope in `δ`. It does **not** cross to
the other side as a single smooth line.

This is a symmetry-controlled termination law. It is the mechanism by which the
enhanced equal-spin symmetry organizes exceptional topology in this system —
not by protecting a degeneracy (that channel is closed by the cross-sector
no-go, `docs/CROSS_SECTOR_EP_NO_GO.md`), but by forcing any diagonal-sector
exceptional line to be mirror-symmetric and to terminate rather than cross.

**Status:** the derivation is unconditional given analyticity and evenness, both
of which are established. Whether any exceptional line exists in a diagonal
sector to exhibit this behaviour is a separate, open, numerical question. The
prediction is recorded before the search so it cannot be retrofitted.

## 4. The equal-spin multiplet count, for all `l` (obligation O7)

At `δ = 0` the spectrum depends on `(m₁,m₂)` only through `m = m₁+m₂` (C6), and
`Λ = l(l+2)` with `l = 2n + |m₁| + |m₂|`, `n ≥ 0` (C5). Fix `l` and `m ≥ 0`
with `|m| ≤ l` and `l ≡ m (mod 2)`. Setting `m₁ = m − m₂`:

* for `m₂ ∈ [0, m]` — that is `m+1` pairs — one has `|m₁|+|m₂| = m`;
* for each `j ≥ 1`, `|m₁|+|m₂| = m + 2j` for **exactly two** pairs, `m₂ = m+j`
  and `m₂ = −j`.

The constraint `|m₁|+|m₂| ≤ l` with `l − |m₁| − |m₂|` even admits `j ≤ (l−m)/2`.
Hence

```
#{(m₁,m₂)} = (m+1) + 2 · (l−m)/2 = l + 1
```

**independently of `m`**. The case `m < 0` follows from `P`. The number of
allowed `m` is `l+1` (namely `m = −l, −l+2, …, l`), so the level total is
`(l+1)²`, the dimension of the SO(4) scalar harmonic level. ∎

This replaces the previous `l ≤ 6` enumeration; O7 is discharged. Verified
exactly to `l = 40`.

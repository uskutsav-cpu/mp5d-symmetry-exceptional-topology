# Conventions

Every statement here is enforced by `tests/unit/test_separation.py`, which
re-derives the separation symbolically in exact rational arithmetic on every CI
run. Nothing in this file is transcribed from a paper without being re-derived.

## 1. Symbols

| Symbol | Meaning |
| --- | --- |
| `M` | Myers–Perry mass parameter, dimension length². Default `M = 1` fixes the length unit. |
| `a`, `b` | Rotation parameters (the two independent spins). |
| `s = (a+b)/2` | Average rotation. |
| `δ = (a−b)/2` | Spin asymmetry. `δ = 0` is the equal-spin surface. |
| `μ` | Scalar field mass. **Never** used for the mass parameter. |
| `ω` | QNM frequency. |
| `m₁`, `m₂` | Azimuthal numbers conjugate to `φ`, `ψ`. Integers. |
| `n` | Angular (polar) node number, `n ≥ 0`. |
| `l = 2n + \|m₁\| + \|m₂\|` | S³ angular momentum. |
| `Λ` | Angular separation constant (see §4 for the exact convention). |
| `N` | Overtone number, `N = 0` fundamental. |

## 2. Metric

Signature `(−++++)`, coordinates `(t, r, θ, φ, ψ)` with `θ ∈ [0, π/2]`,
`φ, ψ ∈ [0, 2π)`:

```
ds² = −dt² + (M/Σ)(dt − a sin²θ dφ − b cos²θ dψ)²
      + (Σ/Δ) dr² + Σ dθ²
      + (r²+a²) sin²θ dφ² + (r²+b²) cos²θ dψ²
```

with

```
Σ = r² + a² cos²θ + b² sin²θ
Π = (r²+a²)(r²+b²)
Δ = (Π − M r²)/r²
```

Verified: `√|g| = r Σ sinθ cosθ` (exact, `test_metric_determinant`).

## 3. Horizons and horizon data

Writing `z = r²` and `P(z) = (z+a²)(z+b²) − M z = (z − z₊)(z − z₋)`:

```
z± = ½ [ (M − a² − b²) ± √( (M − a² − b²)² − 4a²b² ) ]
z₊ z₋ = a²b²          z₊ + z₋ = M − a² − b²
```

* Horizon exists iff `M ≥ (|a| + |b|)²`; extremality at equality.
* `Ω_a = a/(r₊²+a²)`, `Ω_b = b/(r₊²+b²)`.
* `κ = r₊ (r₊² − r₋²) / [(r₊²+a²)(r₊²+b²)]`, `T_H = κ/2π`.
  Checked against `T_H = (D−3)/(4π r₊) = 1/(2π r₊)` at `a = b = 0`.

## 4. Separated massive Klein–Gordon system

Fourier/mode convention:

```
Φ = e^{−iωt + i m₁ φ + i m₂ ψ} R(r) S(θ)
```

`(□ − μ²)Φ = 0` separates exactly. The separation constant `Λ` is fixed by the
requirement `Λ = l(l+2)` at `a = b = 0`.

**Angular equation**

```
(1/(sinθ cosθ)) d/dθ [ sinθ cosθ dS/dθ ]
  + [ Λ + (ω²−μ²)(a² cos²θ + b² sin²θ) − m₁²/sin²θ − m₂²/cos²θ ] S = 0
```

**Radial equation**

```
(1/r) d/dr [ r Δ dR/dr ]
  + [ W²/(r⁴Δ) − G²/r² − (a²+b²)ω² + 2ω(a m₁ + b m₂) − μ² r² − Λ ] R = 0
```

with

```
W(r) = Π ω − m₁ a (r²+b²) − m₂ b (r²+a²)
G    = a b ω − a m₂ − b m₁
```

Sanity limit (`test_schwarzschild_tangherlini_limit`): at `a = b = μ = 0`,
`m₁ = m₂ = 0` the radial bracket reduces to `ω² r⁴/(r²−M) − Λ` with
`Λ = l(l+2)`, the standard 5D Schwarzschild–Tangherlini scalar equation.

### Sign of Im ω

`Φ ∝ e^{−iωt}`, so a **damped** QNM has `Im ω < 0`. Superradiant/unstable modes
have `Im ω > 0`.

## 5. The two derived quantities that organize the project

Substituting `u = cos²θ` and pulling out the constant piece of the potential:

```
Â  ≡ Λ + (ω²−μ²) b²          (shifted angular eigenvalue)
c² ≡ (ω²−μ²)(a²−b²)          (spheroidicity)
```

Because `a² − b² = (a+b)(a−b) = (2s)(2δ)`:

> **c² = 4 s δ (ω² − μ²), exactly.**

Consequences, all exact rather than perturbative:

1. On the equal-spin surface `δ = 0`, `c² ≡ 0` for **every** `s`, `μ`, `ω`.
   The angular operator is then the round S³ Laplacian and
   `Â = l(l+2)`, `l = 2n + |m₁| + |m₂|`, with the full `(l+1)²`-fold SO(4)
   degeneracy.
2. The same holds at `s = 0` (counter-rotating `a = −b`), which is a *second*
   exactly-degenerate surface in the `(s, δ)` plane. The angular degeneracy is
   lifted only by the **product** `s δ`.
3. At `δ = 0` the **radial** equation and the horizon boundary condition depend
   on `(m₁, m₂)` only through the sum `m ≡ m₁ + m₂` (verified symbolically in
   `tests/symmetry/test_equal_spin_u2_structure.py`; at equal spin
   `Ω_a = Ω_b`, so `ω − m₁Ω_a − m₂Ω_b = ω − mΩ`).

Combining 1 and 3: **at `δ = 0` the entire separated QNM problem depends on
`(m₁, m₂)` only via `m`.** The equal-spin degenerate set at fixed `l` is
therefore `{(n, m₁, m₂) : 2n+|m₁|+|m₂| = l, m₁+m₂ = m}`, which has exactly
`l + 1` members for each of the `l + 1` allowed `m` (those with `m ≡ l mod 2`,
`|m| ≤ l`), recovering `(l+1)²`. Each such set is a *single irreducible* SU(2)
multiplet of the `U(2) = (SU(2) × U(1))/ℤ₂` isometry, with `U(1)` charge `m`.

See `docs/SYMMETRY_STRUCTURE.md` for the consequence, which is a structural
no-go on where exceptional points can occur.

## 6. Exchange symmetry

The map `(a, m₁) ↔ (b, m₂)` is an exact symmetry of the full system
(`test_symmetrized_radial_potential_is_exchange_invariant`). At fixed `s` it
acts as `δ → −δ`. In the angular variables it is `u → 1 − u` and gives the
exact relation

```
Â(m₁, m₂, c²) = Â(m₂, m₁, −c²) − c²
```

verified numerically to 1e-9 in `tests/symmetry/test_angular_equal_spin.py`.
This is used throughout as a free consistency check on every continuation in
`δ`: any branch computed at `+δ` must have an exact partner at `−δ` with the
azimuthal labels swapped.

## 7. Angular discretization

With `S = u^{|m₂|/2}(1−u)^{|m₁|/2} f(u)`, `x = 1 − 2u`, the `c² = 0` problem is
Jacobi's equation, so orthonormal Jacobi polynomials `p_n^{(|m₂|,|m₁|)}` are the
*exact* eigenbasis on the equal-spin surface. Multiplication by `u = (1−x)/2` is
tridiagonal, so the angular operator is an exactly tridiagonal
**complex-symmetric** matrix

```
H = D − (c²/2)(I − X),   D_nn = l_n(l_n+2),  l_n = 2n + |m₁| + |m₂|
```

Truncation is the only approximation. Measured convergence: machine precision by
`N ≈ 10` for `|c²| ≲ 100`; geometric decay confirmed up to `|c²| ~ 4×10³`.

Cross-checked against a wholly independent uniform-grid finite-difference
discretization (`mp5d.angular.finite_difference`) sharing no basis, recurrence,
or expansion with the spectral construction.

# Large-n structure of the MP5D recurrence, and the Nollert-type tail

**Status: leading order derived and verified. Higher orders NOT derived.**

## Dominant balance

The raw recurrence is `sum_i c_{n,i} a_{n+1-i} = 0` with

```
c_{n,i} = A_{i+1}(n+1-i)(n-i) + B_i(n+1-i) + C_{i-1}
        = A_{i+1} n^2 + [A_{i+1}(1-2i) + B_i] n + [A_{i+1} i(i-1) + B_i(1-i) + C_{i-1}]
```

so the leading (`n^2`) balance is `sum_i A_{i+1} R^{1-i} = 0` for `R = a_{n+1}/a_n`,
i.e.

> **A(1/R) = 0** — the characteristic roots are the reciprocals of the roots of
> the leading polynomial `A`, which are exactly the ODE's singular points in `x`.

The relevant root is `x = 1` (the irregular point, which sets the radius of
convergence), giving `R -> 1`.

## Why the corrections are half-integer powers

`R = 1` is a **degenerate** root: measured multiplicities of the root `x = 1`
are

| polynomial | multiplicity at `x = 1` |
| --- | --- |
| `A` | >= 5 |
| `B` | exactly 4 |
| `C` | exactly 4 |

(measured at three parameter points including two rotating ones; the clearing
factor contributes `(1-x)^4` and `P2 ~ (1-x)^2`.)

Carrying the expansion to next order, the `n^{-3/2}` equation reads
`-u1 A'(1) = 0`. Because `A'(1) = 0` the coefficient `u1` is **not** forced to
vanish, and the balance moves to the next order — which is precisely the
condition for half-integer powers. Had `A'(1)` been nonzero, `u1 = 0` and the
expansion would run in integer powers of `1/n`.

So

```
R_n = 1 + u1 n^{-1/2} + u2 n^{-1} + u3 n^{-3/2} + ...
```

### Independent confirmation of the exponent

Two checks, neither of them a fit of the coefficients:

1. **Exponent test.** `(R_n - 1) n^{1/2}` converges while `(R_n - 1) n`
   diverges. For the ST5D fundamental: `n^{1/2}` scaling gives
   `-0.599+0.975i`, `-0.579+0.990i` at `n = 200, 400` (converging), whereas
   `n^1` scaling gives `-11.6+19.8i`, `-8.5+13.8i` (growing).
2. **Depth-convergence signature.** The truncation error decays like
   `exp(-k sqrt(N))` — measured earlier as ~5-7 digits per depth doubling,
   which is the hallmark of `exp(-2 sqrt(C n))` coefficient behaviour and is
   incompatible with integer powers.

## Leading coefficient

For coefficient asymptotics `a_n ~ exp(-2 sqrt(C n))`,
`R_n = 1 - sqrt(C/n) + ...`. Matching against the peeled exponential
`exp(c/(1-x))` with `c = i Omega (r_+ - r_-)` gives

```
u1 = -sqrt(-2c),      c = i Omega (r_+ - r_-),   Omega = sqrt(w^2 - mu^2)
```

with the sign fixed by minimality (`Re u1 < 0`, so `|a_n|` decays).

Numerical check, ST5D fundamental: `c = 0.38338 + 0.53384i`, so
`-sqrt(-2c) = -0.52332 + 1.02010i`. The measured `(R_n - 1) sqrt(n)` at
`n = 200, 400` is `-0.5995 + 0.9755i`, `-0.5786 + 0.9903i`, moving toward that
value as the `n^{-1/2}` correction dies away.

### Branch choices

* `Omega = sqrt(w^2 - mu^2)` is taken with `Re Omega >= 0`; the tail therefore
  inherits the same branch cut at `Re w = 0` that the solver already clamps
  away from.
* `sqrt(-2c)` is taken on the principal branch and then sign-flipped if
  `Re u1 > 0`, which selects the minimal (decaying) solution.
* The expansion assumes `r_+ > r_-` strictly. It is **not** valid at
  extremality, where `c -> 0` and the whole balance degenerates.

## What is NOT derived

`u2`, `u3`, ... have **not** been derived. `asymptotic_ratio(..., order > 1)`
raises `NotImplementedError` rather than returning a silently wrong tail. The
measured gain from the leading term alone is modest (see
`docs/TAIL_VALIDATION.md`); the large accelerations reported in the Nollert
literature use several terms.

## A limitation found while doing this

Reading `R_n` cleanly at large `n` is itself limited. The recursive Gaussian
reduction accumulates cancellation: past `n ~ 500-800` the reduced coefficients
degrade, and the backward ratio recursion collapses onto the spurious,
problem-independent value `R_n - 1 = 2/n` (identical for two completely
different parameter sets — the giveaway that it is an artifact). The clean
window used above is `n = 200-400`. This does not affect the CF *value* much,
because deep terms contribute little to the backward evaluation, but it does
cap how far the asymptotics can be probed and it should be fixed (higher working
precision inside the reduction, or a reduction that avoids the long recursive
chain) before `u2` is extracted numerically.

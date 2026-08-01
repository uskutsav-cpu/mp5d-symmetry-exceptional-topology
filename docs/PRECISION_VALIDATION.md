# Arbitrary-precision validation

**Status: the double-precision floor is removed.** Recurrence coefficients, the
angular eigenvalue, the Gaussian reduction, the continued fraction and the root
finder now all run at the requested working precision
(`src/mp5d/radial/highprec.py`). Running only the root solver in high precision
would have proved nothing; the coefficients themselves were the floor.

## What was actually changed

| Stage | Double path | High-precision path |
| --- | --- | --- |
| ODE coefficient functions | NumPy complex128 | `mpmath.mpc` |
| Polynomial extraction | NumPy FFT | explicit `O(N²)` DFT on the unit circle, same contour and half-sample phase offset |
| Angular eigenvalue | LAPACK eigensolve of the Jacobi matrix | Jacobi three-term recurrence solved by the same continued fraction |
| Horizons `r₊`, `r₋` | from `MPGeometry` (double) | recomputed at working precision |
| Reduction / CF / Muller | NumPy | `mpmath` |

The angular step matters: a LAPACK eigensolve would have silently reinjected
double precision into the middle of the chain.

The horizon step was a real leak found during this work. Seeding `rp`, `rm` from
`geo.z_plus` imports a double-precision root and caps everything near 16 digits.
It is now recomputed from `a, b, M` at working precision. (For the
Schwarzschild–Tangherlini case `z₊ = M` exactly, so that case was unaffected —
which is why the leak only showed up once rotating cases were run.)

## Results

Three modes, ladder `(dps, depth) = (30, 400), (50, 800), (80, 1600)`:

| Case | `a` | `b` | `μ` | `l` | stable digits | beyond double? |
| --- | --- | --- | --- | --- | --- | --- |
| `st5d_l0n0` | 0 | 0 | 0 | 0 | **25.8** | yes |
| `hh_table_iii_row1` | 0.2 | 0.3 | 0.1 | 2 | **49.2** | yes |
| `hh_table_ii_massive_singly` | 0.3 | 0 | 0.3 | 4 | **49.1** | yes |

Continued-fraction residuals reach `1e-27` at 30 digits, `1e-47` at 50, `1e-77`
at 80, and `1e-116` at 120 — against a double-precision floor of `1e-12`.

Full ladders, runtimes and residuals: `results/precision_validation.json`.

## The bottleneck is now truncation depth, not arithmetic

This is the load-bearing observation of the section. At fixed `dps = 80`,
varying only the continued-fraction truncation depth for the
Schwarzschild–Tangherlini fundamental:

| depth | CF residual | digits agreeing with previous depth |
| --- | --- | --- |
| 200 | `4.5e-78` | — |
| 400 | `5.1e-77` | 12.9 |
| 800 | `1.1e-76` | 18.3 |
| 1600 | `3.3e-77` | 25.8 |

The residual is flat at `~1e-77` throughout — the *equation* is being solved to
77 digits at every depth — yet the root itself only stabilizes to 26 digits at
depth 1600, gaining roughly 5–7 digits per doubling. So the remaining error is
the difference between the truncated continued fraction and the true infinite
one, not arithmetic.

That is precisely the error a Nollert-style asymptotic tail removes.
**The tail is not implemented**, and this table is the quantitative case for
implementing it: without it, high overtones and weakly damped massive modes will
need impractical depths.

## Reproduction

```bash
PYTHONPATH=src .venv/bin/python scripts/rotating_validation.py
```

Writes `results/precision_validation.json` and `results/rotating_validation.json`.

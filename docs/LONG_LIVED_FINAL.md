# Long-lived and quasi-bound branches

Gate 12. Produced by `scripts/track_long_lived.py`; raw data in
`results/long_lived_branches_final.json`.

## Method

Continuation in the scalar mass `μ` at fixed background, with an **adaptive
step**. A uniform `μ` grid fails near the threshold: `dω/dμ` grows sharply and
the corrector lands on a neighbouring branch — the first atlas build reported
exactly that as a `predictor_miss` at `μ = 1.8`. Here the step is halved
whenever the predictor misses or the frequency moves more than `0.15 |ω|`, down
to `--min-step`, and the walk stops with the reason recorded rather than
accepting a jump.

Regime labels: `ordinary` (`−Im ω > 0.05`), `long_lived` (`−Im ω ≤ 0.05`),
`quasi_bound` (`Re ω < μ`, i.e. past the `ω² = μ²` threshold).

## Result

Eleven branches tracked across six sectors. The damping collapses by up to
**8.6 orders of magnitude** along a single continuous branch.

| sector | `l` | `N` | `s` | `δ` | min `−Im ω` | at `μ` | regimes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (1,1) | 2 | 0 | 0.00 | 0.00 | **8.65e-10** | 2.968 | ordinary → long-lived → quasi-bound |
| (0,0) | 2 | 0 | 0.24 | 0.00 | 5.77e-09 | 2.964 | ordinary → long-lived → quasi-bound |
| (1,0) | 1 | 0 | 0.24 | 0.00 | 5.84e-09 | 1.950 | ordinary → quasi-bound |
| (1,1) | 2 | 0 | 0.24 | 0.10 | 1.61e-07 | 2.964 | ordinary → quasi-bound |
| (1,1) | 2 | 0 | 0.24 | 0.00 | 2.52e-07 | 2.963 | ordinary → quasi-bound |
| (2,0) | 2 | 0 | 0.24 | 0.00 | 1.79e-06 | 2.960 | ordinary → quasi-bound |
| (0,0) | 0 | 0 | 0.24 | 0.00 | 1.95e-04 | 1.217 | ordinary → quasi-bound |
| (0,0) | 0 | 0 | 0.00 | 0.00 | 1.34e-04 | 1.077 | ordinary → quasi-bound |
| (1,1) | 2 | 1 | 0.24 | 0.00 | 9.36e-02 | 3.000 | ordinary → quasi-bound |
| (2,1) | 3 | 0 | 0.24 | 0.00 | 1.32e-01 | 3.000 | ordinary only |
| (2,2) | 4 | 0 | 0.24 | 0.00 | 2.11e-01 | 3.000 | ordinary only |

Representative trajectory, sector (1,1), `l = 2`, `N = 0`, `s = δ = 0`:

| `μ` | `ω` | `−Im ω` | regime |
| --- | --- | --- | --- |
| 0.000 | `1.510567 − 3.575e-01 i` | 3.58e-01 | ordinary |
| 0.900 | `1.632365 − 3.137e-01 i` | 3.14e-01 | ordinary |
| 1.800 | `2.009984 − 1.887e-01 i` | 1.89e-01 | ordinary |
| 2.250 | `2.307313 − 1.045e-01 i` | 1.05e-01 | ordinary |
| 2.700 | `2.691450 − 2.208e-02 i` | 2.21e-02 | quasi-bound |
| 2.949 | `2.948643 − 2.021e-05 i` | 2.02e-05 | quasi-bound |

`Re ω → μ` while `−Im ω → 0`: the standard massive-scalar **quasiresonance**.
The mode approaches the `ω² = μ²` threshold, `Ω = √(ω²−μ²)` turns imaginary,
and the outgoing condition at infinity becomes a decaying one.

Fundamental `N = 0` modes reach the threshold within the tracked range;
overtones (`N = 1`) and higher-`l` branches (`l = 3, 4`) need larger `μ` and
remain ordinary out to `μ = 3`.

## Why every walk ends in `stalled`

`stalled` means the adaptive step reached `--min-step` without a converged,
smoothly-connected root. It is the correct terminus here, not a defect: as
`−Im ω → 0` the mode becomes arbitrarily narrow, the continued fraction needs
ever greater depth to resolve it, and double precision runs out. The reported
minima are therefore **lower bounds on how far the branch was followed**, not
the true infimum of the damping — which is presumably zero at the threshold.

## Relation to the exceptional-point search

Two of these branches (`(1,1) l=2` and `(0,0) l=2`) pass through the long-lived
regime with damping below `1e-8` while a neighbouring overtone remains at
`−Im ω ~ 0.1`. That is a **large** separation in the imaginary direction, so the
quasiresonant transition is not a site of branch coalescence. It is a branch
point in the `μ` dependence — the same conclusion reached earlier for the
`μ ≈ 2.0` transition — and is classified as such, not as an exceptional point.

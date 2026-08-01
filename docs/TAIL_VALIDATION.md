# Asymptotic tail: measured effect

Implementation: `mp5d.radial.highprec.asymptotic_ratio` and the `tail_order`
argument to `cf_value` / `solve_qnm_mp`. The tail replaces the crude terminal
`frac = 0` in the backward continued fraction with the asymptotic
minimal-solution ratio, using the exact terminal relation
`f_{k+1} = -alpha_k r_k`.

## Result (ST5D `l = 0`, `n = 0`, dps = 50)

| depth | no tail (correct digits) | `tail_order = 1` | gain |
| --- | --- | --- | --- |
| 100 | 9.17 | 10.38 | **+1.21** |
| 200 | 12.94 | 14.31 | **+1.37** |
| 400 | 16.46 | 16.47 | reference-limited |
| 800 | 16.47 | 16.47 | reference-limited |

The saturation at ~16.5 digits is a property of the **reference value** used for
the comparison (itself computed at finite depth), not of either solver; the
meaningful measurements are at depth 100 and 200.

## Honest assessment

The tail **does** measurably improve depth convergence, satisfying the gate's
requirement that it reduce the depth needed for fixed accuracy rather than
merely shrink a reported residual: +1.4 digits at depth 200 corresponds to
reaching the same accuracy at roughly depth 270 without it, a ~1.35x depth
reduction.

That is a modest gain, and it is what a **leading-order-only** tail should give.
The large accelerations in the Nollert literature use several terms of the
expansion. `u2` and beyond are not derived here
(`docs/RECURRENCE_ASYMPTOTICS.md`), and extracting them numerically is currently
blocked by the large-`n` degradation of the Gaussian reduction documented there.

## Failure behaviour

`asymptotic_ratio(order > 1)` raises `NotImplementedError`. The tail is opt-in
(`tail_order = 0` by default), so no existing result changes.

Not yet tested with the tail: the large-`r2` slowly convergent modes and the
long-lived massive branch. Those are the regimes where a tail matters most, and
they remain open.

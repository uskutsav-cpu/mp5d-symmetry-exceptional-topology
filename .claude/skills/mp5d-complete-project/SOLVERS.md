# Solver inventory and validated domains

`z_minus` = inner horizon radius squared = the near-extremality difficulty
parameter. Extremality = `M - (|a|+|b|)^2`.

| solver | module | formulation | recurrence-free | validated domain |
| --- | --- | --- | --- | --- |
| A | `radial/qnm.py` (`solve_qnm_cf`) | Leaver continued fraction | no | `z_minus <= 0.20` |
| B | `radial/qnm.py` (`solve_qnm_hill`) | raw recurrence + Hill + Wynn | no | as A |
| C | `radial/solver_c.py` | exterior complex scaling | **yes** | ordinary modes to `z_minus = 0.19` |
| D | `radial/multidomain_near_horizon.py` | multidomain near-horizon spectral | **yes** | ordinary modes only |
| HP | `radial/highprec.py` | full mpmath chain | no | arbitrary precision |

Independence rule: two wrappers around the same discretization are not
independent. A/B share the recurrence family; C and D are the recurrence-free
pair. Any major claim needs at least one recurrence-free confirmation.

Known failure, preserved: at `z_minus = 0.44` no solver converges; best
three-solver spread is `3.5e-3`. Neither the near-horizon scale (refuted by D)
nor double-precision conditioning (refuted by measurement) explains it.

## Seeding

`results/static_seed_table.json` holds static (`a=b=mu=0`) roots indexed by `l`
only -- valid because `Omega_a = Omega_b = 0` and `Lambda = l(l+2)` there.
Overtones are ordered by increasing `|Im omega|` at that point and then carried
by continuation. The continued-fraction **inversion index does not select an
overtone** -- it changes conditioning, not the zero set. Relying on it yields
duplicate "overtones".

`l` must satisfy `l = 2n + |m1| + |m2|`, `n >= 0`. Even `l` needs `(0,0)`-type
labels, odd `l` needs e.g. `(1,0)`; using `(0,0)` for odd `l` silently returns
no roots.

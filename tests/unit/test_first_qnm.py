"""Regression tests on the first trustworthy MP5D quasinormal modes.

These lock in the result of the radial-solver rescue.  If a convention, a
prefactor, a branch choice or the recurrence reduction is ever perturbed, these
fail immediately.
"""

import numpy as np
import pytest

from mp5d.radial.qnm import solve_qnm_cf, solve_qnm_hill

TWO_PI = 2.0 * np.pi

# Matyjasek, arXiv:2107.04815, Tables I-II, converted from omega/T_H with
# T_H = 1/(2 pi).  Metric f(r) = 1 - r^{3-D}, D = 5, so r_+ = 1, matching M = 1.
PUBLISHED = {
    (0, 0, (0, 0)): (3.35418783669, -2.40881848257),
    (0, 1, (0, 0)): (2.33646259225, -8.31019918700),
    (1, 0, (1, 0)): (6.38382253011, -2.27657411582),
    (1, 1, (1, 0)): (5.38079295983, -7.27345089157),
}


@pytest.mark.parametrize("key", list(PUBLISHED))
def test_st5d_scalar_matches_published(key):
    ell, n, (m1, m2) = key
    wref = complex(*PUBLISHED[key]) / TWO_PI
    seed = wref * 1.04 + 0.01  # deliberately not the answer
    # overtones converge more slowly in the truncation depth, so n = 1 needs a
    # deeper schedule to reach the same absolute accuracy as the fundamental
    sched = (100, 200, 400) if n == 0 else (100, 200, 400, 800)
    sol = solve_qnm_cf(
        0.0, 0.0, 0.0, m1, m2, ell, overtone=n,
        initial_frequency=seed, depth_schedule=sched,
    )
    assert sol.converged
    assert abs(sol.omega - wref) < 1e-9, f"got {sol.omega}, published {wref}"
    assert sol.ode_residual < 1e-9
    assert sol.omega.imag < 0, "damped mode must have Im(omega) < 0"


def test_depth_convergence_is_monotone_then_flat():
    wref = complex(*PUBLISHED[(0, 0, (0, 0))]) / TWO_PI
    sol = solve_qnm_cf(
        0.0, 0.0, 0.0, 0, 0, 0, overtone=0,
        initial_frequency=0.60 - 0.30j, depth_schedule=(100, 200, 400, 800),
    )
    errs = [abs(w - wref) for _, w in sol.depth_table]
    assert errs[0] > errs[-1]
    assert errs[-1] < 1e-11
    # stable between the two deepest truncations
    assert abs(sol.depth_table[-1][1] - sol.depth_table[-2][1]) < 1e-11


def test_solver_A_and_B_agree_on_the_fundamental():
    """Independent pair: CF+Gaussian reduction vs raw-recurrence Hill+Wynn."""
    seed = 0.60 - 0.30j
    A = solve_qnm_cf(0.0, 0.0, 0.0, 0, 0, 0, 0, initial_frequency=seed,
                     depth_schedule=(100, 200, 400))
    B = solve_qnm_hill(0.0, 0.0, 0.0, 0, 0, 0, 0, initial_frequency=seed,
                       sizes=(120, 180, 240))
    assert abs(A.omega - B.omega) < 1e-9


@pytest.mark.parametrize("seed", [0.60 - 0.30j, 0.70 - 0.25j, 0.50 - 0.34j])
def test_independent_seeds_reach_the_same_root(seed):
    wref = complex(*PUBLISHED[(0, 0, (0, 0))]) / TWO_PI
    sol = solve_qnm_cf(0.0, 0.0, 0.0, 0, 0, 0, 0, initial_frequency=seed,
                       depth_schedule=(100, 200, 400))
    assert abs(sol.omega - wref) < 1e-10


def test_exchange_symmetry_is_exact_for_two_spins():
    """(a, m1) <-> (b, m2) is an exact identity; this is the sharpest diagnostic."""
    sch = (100, 200, 400)
    P = solve_qnm_cf(0.30, 0.15, 0.20, 1, 0, 1, 0,
                     initial_frequency=1.02 - 0.36j, depth_schedule=sch)
    Q = solve_qnm_cf(0.15, 0.30, 0.20, 0, 1, 1, 0,
                     initial_frequency=1.02 - 0.36j, depth_schedule=sch)
    assert abs(P.omega - Q.omega) < 1e-12


def test_u2_multiplet_degeneracy_is_numerically_exact():
    """Numerical confirmation of claims C6/C7, derived analytically in session 1.

    At ``delta = 0`` the full separated problem depends on ``(m1, m2)`` only
    through ``m = m1 + m2``, so the three ``l = 2, m = 2`` states must share one
    frequency.
    """
    sch = (100, 200, 400)
    omegas = [
        solve_qnm_cf(0.25, 0.25, 0.15, m1, m2, 2, 0,
                     initial_frequency=1.30 - 0.35j, depth_schedule=sch).omega
        for m1, m2 in [(1, 1), (2, 0), (0, 2)]
    ]
    spread = max(abs(x - y) for x in omegas for y in omegas)
    assert spread < 1e-12, f"U(2) multiplet split by {spread:.2e} at delta = 0"


def test_spin_asymmetry_lifts_the_multiplet():
    """The same states must separate once delta != 0 -- otherwise the previous
    test would be passing for a trivial reason."""
    sch = (100, 200, 400)
    omegas = [
        solve_qnm_cf(0.35, 0.15, 0.15, m1, m2, 2, 0,
                     initial_frequency=1.30 - 0.35j, depth_schedule=sch).omega
        for m1, m2 in [(1, 1), (2, 0)]
    ]
    assert abs(omegas[0] - omegas[1]) > 1e-3

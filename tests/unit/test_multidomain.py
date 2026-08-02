"""Solver D: multidomain near-horizon solver."""

import numpy as np
import pytest

from mp5d.radial.multidomain_near_horizon import solve_qnm_multidomain

TWO_PI = 2.0 * np.pi
W_FUND = (3.35418783669 - 2.40881848257j) / TWO_PI
W_L1 = (6.38382253011 - 2.27657411582j) / TWO_PI


def test_solver_d_is_recurrence_free():
    import inspect

    import mp5d.radial.multidomain_near_horizon as md

    code = inspect.getsource(md)
    for forbidden in ("reduce_to_three_term", "continued_fraction",
                      "hill_determinant", "wynn_epsilon", "recurrence_row"):
        assert forbidden not in code, f"Solver D must not use {forbidden}"


@pytest.mark.parametrize("m1,m2,ell,ref", [(0, 0, 0, W_FUND), (1, 0, 1, W_L1)])
def test_static_modes(m1, m2, ell, ref):
    r = solve_qnm_multidomain(0.0, 0.0, 0.0, m1, m2, ell,
                              initial_frequency=ref * 1.01, theta=np.radians(60),
                              Y_m=8.0, n_inner=90, n_outer=200)
    assert r.converged
    assert abs(r.omega - ref) < 1e-4


def test_rotating_two_spin():
    ref = 1.68112 - 0.3472j
    r = solve_qnm_multidomain(0.2, 0.3, 0.1, 1, 1, 2, initial_frequency=ref * 1.01,
                              theta=np.radians(60), Y_m=8.0, n_inner=90, n_outer=200)
    assert abs(r.omega - ref) < 2e-4


@pytest.mark.slow
def test_near_extremal_does_not_converge_and_that_is_recorded():
    """Regression on a documented NEGATIVE result (FAILED_APPROACHES #10).

    At r2 = 0.44 the root wanders with resolution instead of settling.  If a
    future change makes it converge, this test fails and the documentation must
    be revisited.
    """
    vals = []
    for nin in (90, 200):
        r = solve_qnm_multidomain(0.38418, 0.61136, 1.9, 1, 1, 2,
                                  initial_frequency=2.141 - 0.491j,
                                  theta=np.radians(60), Y_m=8.0,
                                  n_inner=nin, n_outer=260)
        vals.append(r.omega)
    assert abs(vals[0] - vals[1]) > 1e-4, (
        "near-extremal point now appears convergent; revisit FAILED_APPROACHES #10"
    )

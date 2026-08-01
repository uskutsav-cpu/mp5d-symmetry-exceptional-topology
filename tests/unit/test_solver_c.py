"""Solver C: complex-scaled spectral solver, recurrence-free.

The key physical test is *invariance*: a genuine quasinormal mode must not move
when the complex-scaling angle, the contour length, or the resolution change,
whereas rotated-continuum artifacts move with the angle.
"""

import numpy as np
import pytest

from mp5d.geometry import MPGeometry
from mp5d.radial.solver_c import (
    SolverCProblem,
    min_scaling_angle,
    solve_qnm_c,
)

TWO_PI = 2.0 * np.pi
W_FUND = (3.35418783669 - 2.40881848257j) / TWO_PI
W_OVER = (2.33646259225 - 8.31019918700j) / TWO_PI


def test_solver_c_uses_no_recurrence_machinery():
    """Independence is structural: assert the module imports none of it."""
    import mp5d.radial.solver_c as sc

    src = sc.__doc__ or ""
    assert "recurrence-free" in src.lower() or "Recurrence-free" in src
    import inspect

    code = inspect.getsource(sc)
    for forbidden in ("reduce_to_three_term", "continued_fraction",
                      "hill_determinant", "wynn_epsilon", "recurrence_row"):
        assert forbidden not in code, f"Solver C must not use {forbidden}"


def test_minimum_scaling_angle_formula():
    """Outgoing decays iff tan(theta) > -Om_I/Om_R."""
    th = min_scaling_angle(W_FUND, 0.0)
    assert 0.60 < th < 0.65          # ~35.7 degrees
    assert abs(np.degrees(th) - 35.68) < 0.1


def test_sigma_min_dips_only_at_the_qnm():
    geo = MPGeometry(a=0.0, b=0.0, M=1.0)
    p = SolverCProblem(geo, 0.0, 0, 0, 0, theta=np.radians(60), L=60.0, resolution=200)
    at_root, _ = p.sigma_min(W_FUND)
    for generic in (0.70 - 0.30j, 0.40 - 0.50j, 0.90 - 0.20j):
        away, _ = p.sigma_min(generic)
        assert away > 1e4 * at_root, f"no contrast at {generic}"
    assert at_root < 1e-12


def test_static_fundamental():
    r = solve_qnm_c(0, 0, 0.0, 0, 0, 0, initial_frequency=W_FUND * 1.02,
                    theta=np.radians(60), L=60.0, resolution=220)
    assert abs(r.omega - W_FUND) < 1e-5


def test_static_overtone_requires_a_large_enough_angle():
    """Below the criterion the outgoing solution grows and the root is wrong.

    This is the theory predicting its own failure mode, so it is pinned.
    """
    need = np.degrees(min_scaling_angle(W_OVER, 0.0))
    assert 74 < need < 75

    bad = solve_qnm_c(0, 0, 0.0, 0, 0, 0, initial_frequency=W_OVER * 1.02,
                      theta=np.radians(60), L=90.0, resolution=300)
    assert abs(bad.omega - W_OVER) > 1e-2, "expected failure below the criterion"

    good = solve_qnm_c(0, 0, 0.0, 0, 0, 0, initial_frequency=W_OVER * 1.02,
                       theta=np.radians(82), L=90.0, resolution=300)
    assert abs(good.omega - W_OVER) < 1e-6


@pytest.mark.parametrize("theta_deg", [45, 55, 65, 75])
def test_physical_mode_is_invariant_under_scaling_angle(theta_deg):
    """A resonance must not move with theta; a continuum artifact would."""
    r = solve_qnm_c(0, 0, 0.0, 0, 0, 0, initial_frequency=W_FUND * 1.02,
                    theta=np.radians(theta_deg), L=60.0, resolution=220)
    assert abs(r.omega - W_FUND) < 1e-5


@pytest.mark.parametrize("L", [40.0, 90.0])
def test_invariant_under_contour_length(L):
    r = solve_qnm_c(0, 0, 0.0, 0, 0, 0, initial_frequency=W_FUND * 1.02,
                    theta=np.radians(60), L=L, resolution=260)
    assert abs(r.omega - W_FUND) < 1e-5


@pytest.mark.parametrize("a,b,mu,m1,m2,ell,ref", [
    (0.2, 0.3, 0.1, 1, 1, 2, 1.68112 - 0.3472j),
    (0.4, 0.2, 0.9, 1, 1, 2, 1.82222 - 0.311152j),
    (0.3, 0.1, 0.1, 1, 1, 4, 2.6344 - 0.349097j),
])
def test_rotating_two_spin_modes(a, b, mu, m1, m2, ell, ref):
    """Huang-Huang Table III, reproduced by a solver with no recurrence in it."""
    r = solve_qnm_c(a, b, mu, m1, m2, ell, initial_frequency=ref * 1.01,
                    theta=np.radians(55), L=70.0, resolution=280)
    assert abs(r.omega - ref) < 2e-4


@pytest.mark.slow
def test_large_r2_beyond_the_published_cfm_limit():
    """r2 = 0.140, past the r2 ~ 0.1 limit Huang-Huang state for their CFM."""
    ref = 2.77349 - 0.332904j
    r = solve_qnm_c(0.3, 0.4, 0.1, 1, 1, 4, initial_frequency=ref * 1.01,
                    theta=np.radians(55), L=70.0, resolution=280)
    assert abs(r.omega - ref) < 3e-4

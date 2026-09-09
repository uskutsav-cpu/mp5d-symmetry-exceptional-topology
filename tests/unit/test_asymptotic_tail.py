"""Tests for the Nollert-type asymptotic tail."""

import mpmath as mp
import pytest

from mp5d.geometry import MPGeometry
from mp5d.radial.highprec import (
    HighPrecisionProblem,
    asymptotic_ratio,
    solve_qnm_mp,
)

REF_ST5D = (
    "0.533835574268276465288894914625582866112313851"
    " - 0.383375368512460197238361437790785932867003806j"
)


def test_tail_order_zero_is_the_identity():
    assert asymptotic_ratio(100, mp.mpc(0.4, 0.5), order=0) == 1


def test_undserived_orders_fail_loudly():
    """Better to raise than to return a silently wrong tail."""
    with pytest.raises(NotImplementedError):
        asymptotic_ratio(100, mp.mpc(0.4, 0.5), order=2)


def test_leading_coefficient_matches_minus_sqrt_minus_2c():
    """u1 = -sqrt(-2c), sign fixed by minimality (Re u1 < 0)."""
    with mp.workdps(40):
        c = mp.mpc("0.3833753685", "0.5338355743")
        n = 10_000
        u1 = (asymptotic_ratio(n, c, 1) - 1) * mp.sqrt(n)
        expect = -mp.sqrt(-2 * c)
        if expect.real > 0:
            expect = -expect
        assert abs(u1 - expect) < mp.mpf("1e-25")
        assert u1.real < 0, "minimal solution requires Re(u1) < 0"


def test_tail_ratio_tends_to_one():
    with mp.workdps(30):
        c = mp.mpc("0.38", "0.53")
        assert abs(asymptotic_ratio(10**8, c, 1) - 1) < mp.mpf("1e-3")


def test_exp_coefficient_is_i_Omega_times_horizon_gap():
    with mp.workdps(30):
        geo = MPGeometry(a=0.2, b=0.3, M=1.0)
        p = HighPrecisionProblem(geo, 0.1, 1, 1, 2, dps=30)
        om = mp.mpc("1.681109", "-0.347227")
        c = p.exp_coefficient(om)
        Om = mp.sqrt(om**2 - p.mu**2)
        if Om.real < 0:
            Om = -Om
        assert abs(c - mp.mpc(0, 1) * Om * (p.rp - p.rm)) < mp.mpf("1e-25")


@pytest.mark.slow
def test_tail_improves_depth_convergence():
    """The gate requirement: fewer terms for the same accuracy, not just a
    smaller reported residual."""
    with mp.workdps(50):
        ref = mp.mpmathify(REF_ST5D)

    def digits(depth, tail_order):
        s = solve_qnm_mp(
            0.0,
            0.0,
            0.0,
            0,
            0,
            0,
            0,
            initial_frequency=0.53 - 0.38j,
            dps=50,
            depth=depth,
            tail_order=tail_order,
        )
        with mp.workdps(50):
            d = abs(mp.mpmathify(s.omega_str) - ref)
            return float(-mp.log10(d)) if d > 0 else 50.0

    for depth in (100, 200):
        assert digits(depth, 1) > digits(depth, 0) + 0.5, (
            f"tail gave no measurable gain at depth {depth}"
        )


@pytest.mark.slow
def test_tail_does_not_move_the_converged_root():
    """At large depth the tail must not shift a converged answer."""
    base = solve_qnm_mp(
        0.0, 0.0, 0.0, 0, 0, 0, 0, initial_frequency=0.53 - 0.38j, dps=40, depth=600, tail_order=0
    )
    tailed = solve_qnm_mp(
        0.0, 0.0, 0.0, 0, 0, 0, 0, initial_frequency=0.53 - 0.38j, dps=40, depth=600, tail_order=1
    )
    assert abs(base.omega - tailed.omega) < 1e-12

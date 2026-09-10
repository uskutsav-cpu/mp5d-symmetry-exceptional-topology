"""Validate the certification pipeline on systems with known answers.

A certificate is worthless if the certifier cannot be shown to reject as well as
accept.  Every routine is therefore exercised on both sides: a box that does
contain a unique root, and a box that provably contains none.
"""

from __future__ import annotations

import numpy as np
import pytest

from mp5d.certification.intervals import (
    arb_available,
    enclose_angular_eigenvalue,
    interval_newton_exclude,
    krawczyk_test,
)

# --------------------------------------------------------------------------
# Krawczyk
# --------------------------------------------------------------------------


def _circle_system():
    """f(x,y) = (x^2 + y^2 - 1, x - y); roots at +-(1,1)/sqrt(2)."""

    def f(v):
        x, y = v
        return np.array([x * x + y * y - 1.0, x - y])

    def jac(c, r):
        x, y = c
        rx, ry = r
        lo = np.array([[2 * (x - rx), 2 * (y - ry)], [1.0, -1.0]])
        hi = np.array([[2 * (x + rx), 2 * (y + ry)], [1.0, -1.0]])
        return np.minimum(lo, hi), np.maximum(lo, hi)

    return f, jac


def test_krawczyk_certifies_a_unique_root():
    f, jac = _circle_system()
    root = 1.0 / np.sqrt(2.0)
    res = krawczyk_test(f, jac, centre=[root + 1e-4, root - 1e-4],
                        radius=[1e-2, 1e-2])
    assert res.contains_unique_root, res.detail
    assert not res.excludes_root
    assert res.contraction < 1.0


def test_krawczyk_excludes_a_root_free_box():
    f, jac = _circle_system()
    res = krawczyk_test(f, jac, centre=[5.0, 5.0], radius=[1e-2, 1e-2])
    assert res.excludes_root, res.detail
    assert not res.contains_unique_root


def test_krawczyk_rejects_a_midpoint_only_jacobian():
    """A midpoint Jacobian is a heuristic, not a certificate."""
    f, _ = _circle_system()

    def bad_jac(c, r):
        x, y = c
        return np.array([[2 * x, 2 * y], [1.0, -1.0]]), np.array([1.0])

    with pytest.raises(ValueError):
        krawczyk_test(f, bad_jac, centre=[0.7, 0.7], radius=[1e-3, 1e-3])


def test_krawczyk_reports_inconclusive_rather_than_guessing():
    """A box far too large for contraction must not be certified either way."""
    f, jac = _circle_system()
    res = krawczyk_test(f, jac, centre=[0.5, 0.5], radius=[2.0, 2.0])
    assert res.inconclusive
    assert not res.contains_unique_root


# --------------------------------------------------------------------------
# Interval exclusion
# --------------------------------------------------------------------------


def test_interval_exclusion_detects_root_free_box():
    def f_int(c, r):
        # f(x) = x^2 + 1 over [c-r, c+r]; always >= 1
        lo = np.array([max(0.0, abs(c[0]) - r[0]) ** 2 + 1.0])
        hi = np.array([(abs(c[0]) + r[0]) ** 2 + 1.0])
        return lo, hi

    assert interval_newton_exclude(f_int, [0.0], [1.0])


def test_interval_exclusion_is_silent_when_a_root_may_be_present():
    def f_int(c, r):
        return np.array([c[0] - r[0]]), np.array([c[0] + r[0]])

    assert not interval_newton_exclude(f_int, [0.0], [1.0])


# --------------------------------------------------------------------------
# Arb enclosure of the angular eigenvalue
# --------------------------------------------------------------------------


@pytest.mark.skipif(not arb_available(), reason="python-flint not installed")
def test_arb_encloses_the_exact_equal_spin_eigenvalue():
    """At ``c2 = 0`` the answer is exactly ``l(l+2)``; the ball must contain it.

    This is the one place where an independent exact value exists, so it is the
    only honest way to validate the enclosure machinery.
    """
    for (m1, m2, k) in [(1, 1, 0), (0, 0, 1), (2, 0, 0)]:
        ell = 2 * k + abs(m1) + abs(m2)
        out = enclose_angular_eigenvalue(m1, m2, k, 0.0, N=40, prec=200)
        exact = ell * (ell + 2)
        assert abs(out["eigenvalue_mid"].real - exact) < 1e-8, (m1, m2, k, out)
        assert abs(out["eigenvalue_mid"].imag) < 1e-8
        assert "NOT a continuum statement" in out["scope"]


@pytest.mark.skipif(not arb_available(), reason="python-flint not installed")
def test_enclosure_scope_is_labelled_not_promoted():
    out = enclose_angular_eigenvalue(1, 1, 0, 0.3 + 0.1j, N=30, prec=160)
    assert "E_truncation" in out["scope"]
    assert out["truncation_N"] == 30

import math

import pytest

from mp5d.geometry import MPGeometry, ab_to_sdelta, sdelta_to_ab


def test_schwarzschild_tangherlini_limit():
    g = MPGeometry(a=0.0, b=0.0, M=1.0)
    assert g.r_plus == pytest.approx(1.0)
    assert g.r_minus == pytest.approx(0.0)
    # T_H = (D-3)/(4 pi r_+) = 1/(2 pi r_+) in D=5
    assert g.T_H == pytest.approx(1.0 / (2 * math.pi * g.r_plus))
    assert g.Omega_a == 0.0 and g.Omega_b == 0.0


def test_horizon_polynomial_roots():
    g = MPGeometry(a=0.31, b=0.17, M=1.0)
    for z in (g.z_plus, g.z_minus):
        assert g.P(z) == pytest.approx(0.0, abs=1e-13)
    # Vieta: z_+ z_- = a^2 b^2 and z_+ + z_- = M - a^2 - b^2
    assert g.z_plus * g.z_minus == pytest.approx((g.a * g.b) ** 2)
    assert g.z_plus + g.z_minus == pytest.approx(g.M - g.a**2 - g.b**2)


def test_extremality_boundary():
    # M = (|a| + |b|)^2 is extremal: z_+ = z_- , kappa = 0
    a, b = 0.3, 0.2
    g = MPGeometry(a=a, b=b, M=(a + b) ** 2)
    assert g.extremality == pytest.approx(0.0)
    assert g.z_plus == pytest.approx(g.z_minus)
    assert g.kappa == pytest.approx(0.0, abs=1e-12)


def test_exchange_symmetry_of_horizon_quantities():
    g = MPGeometry(a=0.41, b=0.13)
    h = g.exchanged()
    assert h.z_plus == pytest.approx(g.z_plus)
    assert h.kappa == pytest.approx(g.kappa)
    assert h.Omega_a == pytest.approx(g.Omega_b)
    assert h.Omega_b == pytest.approx(g.Omega_a)
    # delta -> -delta at fixed s
    assert h.s == pytest.approx(g.s)
    assert h.delta == pytest.approx(-g.delta)


def test_sdelta_roundtrip():
    for a, b in [(0.3, 0.1), (0.25, 0.25), (0.0, 0.4)]:
        s, d = ab_to_sdelta(a, b)
        aa, bb = sdelta_to_ab(s, d)
        assert (aa, bb) == pytest.approx((a, b))
        g = MPGeometry.from_sdelta(s, d)
        assert (g.a, g.b) == pytest.approx((a, b))


def test_equal_spin_surface():
    g = MPGeometry.from_sdelta(s=0.3, delta=0.0)
    assert g.a == pytest.approx(g.b)
    assert g.Omega_a == pytest.approx(g.Omega_b)
    # z_+ z_- = a^4 at equal spin
    assert g.z_plus * g.z_minus == pytest.approx(g.a**4)

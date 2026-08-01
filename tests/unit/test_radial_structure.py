"""Exact checks on the radial equation's singular structure."""

import cmath

import pytest
import sympy as sp

from mp5d.geometry import MPGeometry
from mp5d.radial.structure import (
    ASYMPTOTIC_POWER,
    G_const,
    W_of_z,
    horizon_exponent,
)

z, r = sp.symbols("z r", positive=True)
a, b, M, w, mu, m1, m2, Lam = sp.symbols("a b M omega mu m1 m2 Lambda")

P = (z + a**2) * (z + b**2) - M * z
W = (z + a**2) * (z + b**2) * w - m1 * a * (z + b**2) - m2 * b * (z + a**2)
G = a * b * w - a * m2 - b * m1
c0 = -Lam - (a**2 + b**2) * w**2 + 2 * w * (a * m1 + b * m2)


def test_origin_is_an_ordinary_point():
    """The 1/z terms cancel identically: W(0) = a b G and P(0) = a^2 b^2."""
    assert sp.expand(W.subs(z, 0) - a * b * G) == 0
    assert sp.expand(P.subs(z, 0) - a**2 * b**2) == 0
    # residue of [W^2/(z P) - G^2/z] at z = 0
    residue = sp.cancel(W.subs(z, 0) ** 2 / P.subs(z, 0) - G**2)
    assert residue == 0


def test_potential_is_even_in_r_so_no_coulomb_tail():
    """No 1/r term at large r, hence the exponent is exactly -3/2."""
    V = W**2 / (z * P) - G**2 / z + c0 - mu**2 * z
    V_r = V.subs(z, r**2)
    # V(r) = V(-r) identically  =>  the 1/r expansion has only even powers
    assert sp.cancel(sp.together(V_r - V_r.subs(r, -r))) == 0

    # confirm the exponent: R ~ e^{i Om r} r^beta solves
    # R'' + (3/r) R' + Om^2 R = 0 to O(1/r) only for beta = -3/2
    Om, beta = sp.symbols("Om beta")
    R = sp.exp(sp.I * Om * r) * r**beta
    expr = sp.diff(R, r, 2) + 3 / r * sp.diff(R, r) + Om**2 * R
    lead = sp.simplify(sp.cancel(sp.expand(expr / R) * r))
    sol = sp.solve(sp.limit(lead, r, sp.oo), beta)
    assert sol == [sp.Rational(-3, 2)]
    assert ASYMPTOTIC_POWER == -1.5


def test_horizon_indicial_exponent_derivation():
    """Frobenius at z_+, derived from the ODE itself rather than asserted.

    Substitute ``z = z_+ + x`` and ``R = x^s`` into
    ``4[P R'' + P' R'] + [W^2/(zP) - G^2/z + c0 - mu^2 z] R = 0``
    and read off the coefficient of the leading power ``x^(s-1)``.
    """
    zp, zm, x, s = sp.symbols("z_p z_m x s", positive=True)
    Wp = sp.Symbol("W_p")

    Pf = x * (x + zp - zm)                       # P(z_+ + x)
    Wf = Wp + sp.Symbol("W1") * x                # W is analytic at z_+
    zf = zp + x
    R = x**s

    ode = 4 * (Pf * sp.diff(R, x, 2) + sp.diff(Pf, x) * sp.diff(R, x)) + (
        Wf**2 / (zf * Pf)
        - sp.Symbol("Gsq") / zf
        + sp.Symbol("c0")
        - sp.Symbol("musq") * zf
    ) * R

    # leading behaviour is x^(s-1); its coefficient is the indicial polynomial
    indicial = sp.simplify(sp.limit(sp.expand(ode / x ** (s - 1)), x, 0))
    roots = sp.solve(sp.Eq(indicial, 0), s)
    expected = {
        sp.I * Wp / (2 * sp.sqrt(zp) * (zp - zm)),
        -sp.I * Wp / (2 * sp.sqrt(zp) * (zp - zm)),
    }
    assert {sp.simplify(sp.radsimp(rt)) for rt in roots} == {
        sp.simplify(sp.radsimp(e)) for e in expected
    }


def test_W_at_horizon_factorizes_through_angular_velocities():
    """W(z_+) = (z_+ + a^2)(z_+ + b^2) (w - m1 Om_a - m2 Om_b)."""
    zp = sp.Symbol("z_p", positive=True)
    Wp = W.subs(z, zp)
    Om_a = a / (zp + a**2)
    Om_b = b / (zp + b**2)
    target = (zp + a**2) * (zp + b**2) * (w - m1 * Om_a - m2 * Om_b)
    assert sp.expand(sp.cancel(sp.together(Wp - target))) == 0


@pytest.mark.parametrize(
    "aa,bb,m1v,m2v", [(0.3, 0.2, 1, 0), (0.25, 0.25, 2, 2), (0.4, 0.1, -1, 3)]
)
def test_horizon_exponent_equals_surface_gravity_form(aa, bb, m1v, m2v):
    """Numeric identity sigma_+ = W(z_+)/(2 r_+ (z_+-z_-)) = (w - m.Om)/(2 kappa)."""
    geo = MPGeometry(a=aa, b=bb)
    omega = 0.83 - 0.21j
    Wp = W_of_z(geo, geo.z_plus, omega, m1v, m2v)
    direct = Wp / (2 * geo.r_plus * (geo.z_plus - geo.z_minus))
    assert direct == pytest.approx(horizon_exponent(geo, omega, m1v, m2v))


def test_horizon_exponent_respects_equal_spin_m_sum_rule():
    """At delta = 0 the horizon exponent depends only on m1 + m2."""
    geo = MPGeometry.from_sdelta(s=0.3, delta=0.0)
    omega = 0.7 - 0.15j
    s1 = horizon_exponent(geo, omega, 3, -1)
    s2 = horizon_exponent(geo, omega, 1, 1)
    s3 = horizon_exponent(geo, omega, 0, 2)
    assert s1 == pytest.approx(s2) == pytest.approx(s3)


def test_G_matches_symbolic():
    geo = MPGeometry(a=0.31, b=0.17)
    omega = 0.5 + 0.2j
    got = G_const(geo, omega, 2, -1)
    want = complex(
        G.subs({a: geo.a, b: geo.b, w: omega, m1: 2, m2: -1}).evalf()
    )
    assert got == pytest.approx(want)


def test_superradiant_threshold_is_where_horizon_exponent_vanishes():
    """sigma_+ = 0 exactly at w = m1 Om_a + m2 Om_b, the superradiance bound."""
    geo = MPGeometry(a=0.35, b=0.15)
    m1v, m2v = 2, 1
    w_c = m1v * geo.Omega_a + m2v * geo.Omega_b
    assert abs(horizon_exponent(geo, complex(w_c), m1v, m2v)) < 1e-14
    assert cmath.isclose(
        horizon_exponent(geo, complex(w_c + 1e-6), m1v, m2v).real,
        1e-6 / (2 * geo.kappa),
        rel_tol=1e-6,
    )

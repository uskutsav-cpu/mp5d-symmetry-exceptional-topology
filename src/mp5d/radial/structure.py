"""Singular structure and characteristic exponents of the radial equation.

Every radial method (continued fraction, collocation, shooting) needs these, so
they are derived once here and pinned symbolically in
``tests/unit/test_radial_structure.py``.

In ``z = r^2`` the radial equation of ``docs/CONVENTIONS.md`` §4 reads

    4 [ P(z) R'' + P'(z) R' ] + [ W^2/(z P) - G^2/z + c0 - mu^2 z ] R = 0

with

    P(z) = (z + a^2)(z + b^2) - M z = (z - z_+)(z - z_-)
    W(z) = (z + a^2)(z + b^2) w - m1 a (z + b^2) - m2 b (z + a^2)
    G    = a b w - a m2 - b m1
    c0   = -Lambda - (a^2 + b^2) w^2 + 2 w (a m1 + b m2)

Singular structure
------------------
* ``z = 0`` (i.e. ``r = 0``) is an **ordinary point** whenever ``a b != 0``.
  The two apparently singular terms cancel identically, because
  ``W(0) = a b G`` and ``P(0) = a^2 b^2`` give
  ``W(0)^2 / P(0) - G^2 = 0`` exactly.  If either spin vanishes the
  cancellation still holds but ``P(0) = 0`` and ``z = 0`` becomes singular, so
  the singly-rotating and non-rotating cases need separate treatment.
* ``z = z_+`` and ``z = z_-`` are regular singular points with exponents
  ``-/+ i sigma_+`` and ``-/+ i sigma_-``.
* ``z = infinity`` is irregular; in ``r`` the solutions behave as
  ``e^{+/- i Omega r} r^{-3/2}`` with ``Omega = sqrt(w^2 - mu^2)``.

So the problem has **two** regular singular points plus an irregular point at
infinity: it is of confluent Heun type, which is what makes a Leaver-style
three-term recurrence in ``x = (z - z_+)/(z - z_-)`` viable.

Note on the massive asymptotics
-------------------------------
The radial potential is an even function of ``r``, so its large-``r`` expansion
contains **no** ``1/r`` term.  Consequently there is no Coulomb-like
logarithmic phase: the exponent is exactly ``-3/2``, with corrections starting
at ``O(1/r)`` from the even structure only.  This is a genuine simplification
relative to the four-dimensional massive-scalar problem, where the ``1/r`` tail
of the Newtonian potential forces an ``r^{i nu}`` factor.  It follows from
``D - 3 = 2``.
"""

from __future__ import annotations

import cmath

from ..geometry import MPGeometry

__all__ = [
    "horizon_exponent",
    "inner_horizon_exponent",
    "Omega_infinity",
    "W_of_z",
    "G_const",
    "c0_const",
    "ASYMPTOTIC_POWER",
]

#: Power of ``r`` multiplying ``exp(+/- i Omega r)`` at spatial infinity.
ASYMPTOTIC_POWER = -1.5


def G_const(geo: MPGeometry, omega: complex, m1: int, m2: int) -> complex:
    """``G = a b w - a m2 - b m1``."""
    return geo.a * geo.b * omega - geo.a * m2 - geo.b * m1


def W_of_z(geo: MPGeometry, z: complex, omega: complex, m1: int, m2: int) -> complex:
    a, b = geo.a, geo.b
    return (z + a**2) * (z + b**2) * omega - m1 * a * (z + b**2) - m2 * b * (z + a**2)


def c0_const(
    geo: MPGeometry, omega: complex, m1: int, m2: int, Lambda: complex
) -> complex:
    a, b = geo.a, geo.b
    return -Lambda - (a**2 + b**2) * omega**2 + 2 * omega * (a * m1 + b * m2)


def horizon_exponent(
    geo: MPGeometry, omega: complex, m1: int, m2: int
) -> complex:
    """``sigma_+ = (w - m1 Omega_a - m2 Omega_b) / (2 kappa)``.

    The physical (ingoing) QNM solution behaves as ``(z - z_+)^{-i sigma_+}``.

    This identity is *derived*, not assumed: the indicial equation at ``z_+``
    gives ``sigma_+ = W(z_+) / (2 r_+ (z_+ - z_-))``, and
    ``W(z_+) = (z_+ + a^2)(z_+ + b^2)(w - m1 Omega_a - m2 Omega_b)``, which
    combines with ``kappa = r_+ (z_+ - z_-)/[(z_+ + a^2)(z_+ + b^2)]`` to give
    the expression above.  See ``test_horizon_exponent_equals_surface_gravity_form``.
    """
    return (omega - m1 * geo.Omega_a - m2 * geo.Omega_b) / (2.0 * geo.kappa)


def inner_horizon_exponent(
    geo: MPGeometry, omega: complex, m1: int, m2: int
) -> complex:
    """Indicial exponent scale at the inner horizon ``z_-``.

    ``sigma_- = W(z_-) / (2 r_- (z_- - z_+))``; solutions go as
    ``(z - z_-)^{+/- i sigma_-}``.  Needed by the Leaver ansatz, not by the
    boundary conditions.
    """
    zm, zp = geo.z_minus, geo.z_plus
    W = W_of_z(geo, zm, omega, m1, m2)
    return W / (2.0 * cmath.sqrt(zm) * (zm - zp))


def Omega_infinity(omega: complex, mu: float) -> complex:
    """``Omega = sqrt(w^2 - mu^2)``, the wavenumber at spatial infinity.

    Branch: chosen so that ``Re Omega >= 0``, which selects ``e^{+i Omega r}``
    as the outgoing solution for the QNM boundary condition.  For ``|w| < mu``
    this becomes the bound-state / quasi-bound branch and the caller must decide
    which sheet is wanted; that choice is recorded per-run, not here.
    """
    Om = cmath.sqrt(omega**2 - mu**2)
    return Om if Om.real >= 0 else -Om

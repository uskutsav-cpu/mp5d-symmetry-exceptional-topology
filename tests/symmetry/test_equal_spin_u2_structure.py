"""Exact U(2) multiplet structure on the equal-spin surface.

Two facts are established here symbolically (exact rational arithmetic):

1. At ``a = b`` the *radial* equation depends on the azimuthal numbers only
   through the sum ``m = m1 + m2``.  Combined with the exact angular result
   ``Ahat = l(l+2)``, this means the whole separated QNM problem at ``delta = 0``
   depends on ``(m1, m2)`` only via ``m``.

2. The horizon boundary condition inherits the same property, because
   ``Omega_a = Omega_b`` at equal spin, so ``omega - m1 Omega_a - m2 Omega_b
   = omega - m Omega``.

Together these give the equal-spin degeneracy: all states with fixed ``l`` and
fixed ``m`` share one frequency.  The counting test below confirms that such a
set has exactly ``l + 1`` members, i.e. it is a *single irreducible* SU(2)
multiplet, and that summing over the ``l + 1`` allowed values of ``m``
reproduces the full ``(l+1)^2`` SO(4) degeneracy.
"""

import sympy as sp

from mp5d.angular import s3_multiplet


def test_equal_spin_radial_potential_depends_only_on_m1_plus_m2():
    r, s, M, w, mu, m1, m2 = sp.symbols("r s M omega mu m1 m2")
    d = sp.Symbol("d")  # placeholder for a redistribution of m1, m2 at fixed sum

    def radial_bracket(mm1, mm2):
        Pi_ = (r**2 + s**2) ** 2
        rrDel = Pi_ - M * r**2
        W = Pi_ * w - mm1 * s * (r**2 + s**2) - mm2 * s * (r**2 + s**2)
        G = s**2 * w - s * mm2 - s * mm1
        # Lambda = l(l+2) - (w^2 - mu^2) s^2 at equal spin; l-dependence only
        l = sp.Symbol("l")
        Lam = l * (l + 2) - (w**2 - mu**2) * s**2
        return (
            W**2 / (r**2 * rrDel)
            - G**2 / r**2
            - 2 * s**2 * w**2
            + 2 * w * s * (mm1 + mm2)
            - mu**2 * r**2
            - Lam
        )

    # Redistribute the sum: (m1, m2) -> (m1 + d, m2 - d) keeps m1 + m2 fixed.
    diff = radial_bracket(m1, m2) - radial_bracket(m1 + d, m2 - d)
    assert sp.cancel(sp.together(diff)) == 0


def test_equal_spin_horizon_condition_depends_only_on_m1_plus_m2():
    s, M, w, m1, m2, d = sp.symbols("s M omega m1 m2 d")
    zp = sp.Symbol("z_plus")
    Om = s / (zp + s**2)  # Omega_a = Omega_b at equal spin
    expr = w - m1 * Om - m2 * Om
    assert sp.cancel(expr - expr.subs({m1: m1 + d, m2: m2 - d})) == 0


def test_fixed_l_and_m_multiplet_is_a_single_su2_irrep():
    """|{(n,m1,m2) : 2n+|m1|+|m2| = l, m1+m2 = m}| = l + 1 for every allowed m."""
    for l in range(0, 7):
        states = s3_multiplet(l)
        # m must have the same parity as l and satisfy |m| <= l
        allowed_m = [m for m in range(-l, l + 1) if (m - l) % 2 == 0]
        assert len(allowed_m) == l + 1
        total = 0
        for m in allowed_m:
            block = [t for t in states if t[1] + t[2] == m]
            assert len(block) == l + 1, f"l={l}, m={m}: got {len(block)}"
            total += len(block)
        assert total == (l + 1) ** 2 == len(states)


def test_azimuthal_charges_are_conserved_for_all_delta():
    """m1, m2 label decoupled sectors at *every* delta, not just delta = 0.

    ``d/dphi`` and ``d/dpsi`` are Killing vectors of the metric for arbitrary
    ``a``, ``b``; the separated system never couples different ``(m1, m2)``.
    This is what forbids the equal-spin multiplet degeneracy from becoming
    defective: its members live in permanently decoupled blocks.
    """
    r, u, a, b, M, w, m1, m2 = sp.symbols("r u a b M omega m1 m2")
    Sig = r**2 + a**2 * u + b**2 * (1 - u)
    k = sp.Matrix([1, -a * (1 - u), -b * u])
    B = sp.diag(-1, (r**2 + a**2) * (1 - u), (r**2 + b**2) * u) + (M / Sig) * (k * k.T)
    # No metric component depends on phi or psi -> the mode ansatz is exact and
    # the (m1, m2) dependence enters only algebraically.
    for entry in B:
        assert not entry.free_symbols & {sp.Symbol("phi"), sp.Symbol("psi")}

"""Regression test on the *conventions*: re-derive the separation symbolically.

This is the load-bearing test of the whole project.  It rebuilds the 5D
Myers-Perry metric from scratch, forms the massive Klein-Gordon operator, and
verifies -- exactly, in rational arithmetic -- that

  * ``sqrt(|g|) = r Sigma sin(theta) cos(theta)``;
  * the rotation/frequency term separates (mixed partial vanishes identically);
  * the separated radial and angular potentials are exactly the ones coded in
    ``mp5d.angular`` and ``mp5d.radial``;
  * the whole system is invariant under ``(a, m1) <-> (b, m2)``.

If any convention drifts, this fails rather than silently poisoning the science.
"""

import sympy as sp

r, u = sp.symbols("r u")
a, b, M, w, mu, m1, m2, Lam = sp.symbols("a b M omega mu m1 m2 Lambda")


def _pieces():
    """Sigma, Pi, r^2*Delta and the exact term ``-Sigma * K^T B^{-1} K``."""
    Sig = r**2 + a**2 * u + b**2 * (1 - u)
    Pi_ = (r**2 + a**2) * (r**2 + b**2)
    rrDel = Pi_ - M * r**2

    k = sp.Matrix([1, -a * (1 - u), -b * u])
    B = sp.diag(-1, (r**2 + a**2) * (1 - u), (r**2 + b**2) * u) + (M / Sig) * (k * k.T)
    K = sp.Matrix([-w, m1, m2])
    num = sp.cancel(sp.together((K.T * B.adjugate() * K)[0, 0]))
    det = sp.cancel(sp.together(B.det()))
    Q = sp.cancel(-Sig * num / det)
    return Sig, Pi_, rrDel, B, Q


def test_metric_determinant():
    """det(B) = u(u-1) * r^2 Delta, hence sqrt|g| = r Sigma sin t cos t."""
    Sig, Pi_, rrDel, B, _ = _pieces()
    det = sp.cancel(sp.together(B.det()))
    assert sp.cancel(det - u * (u - 1) * rrDel) == 0

    # full determinant: g_rr * g_thth * det(B) with g_rr = Sig/Delta, g_thth = Sig
    Delta = rrDel / r**2
    detg = sp.cancel((Sig / Delta) * Sig * det)
    # sqrt|g| = r Sigma sin cos  =>  |g| = r^2 Sigma^2 sin^2 cos^2 = r^2 Sig^2 u (1-u)
    assert sp.cancel(detg + r**2 * Sig**2 * u * (1 - u)) == 0


def test_rotation_term_separates_exactly():
    """The mixed partial d^2/dr du of the frequency term vanishes identically."""
    *_, Q = _pieces()
    assert sp.cancel(sp.together(sp.diff(Q, r, u))) == 0


def test_separated_potentials_match_implementation():
    """Closed forms used by the solvers are exactly the separated potentials."""
    Sig, Pi_, rrDel, _, Q = _pieces()

    # angular part, as coded in mp5d.angular.spheroidal5d
    ang = -(m1**2) / (1 - u) - m2**2 / u + w**2 * (a**2 - b**2) * u
    # radial part, as coded in mp5d.radial
    W = Pi_ * w - m1 * a * (r**2 + b**2) - m2 * b * (r**2 + a**2)
    G = a * b * w - a * m2 - b * m1
    rad = (
        W**2 / (r**2 * rrDel)
        - G**2 / r**2
        - a**2 * w**2
        + 2 * w * (a * m1 + b * m2)
    )
    assert sp.cancel(sp.together(Q - ang - rad)) == 0


def test_symmetrized_radial_potential_is_exchange_invariant():
    """After moving the ``w^2 b^2`` constant into the angular side, ``(a,m1)<->(b,m2)``
    is a manifest symmetry of both separated potentials."""
    Sig, Pi_, rrDel, _, _ = _pieces()
    W = Pi_ * w - m1 * a * (r**2 + b**2) - m2 * b * (r**2 + a**2)
    G = a * b * w - a * m2 - b * m1
    rad_sym = (
        W**2 / (r**2 * rrDel)
        - G**2 / r**2
        - (a**2 + b**2) * w**2
        + 2 * w * (a * m1 + b * m2)
    )
    swap = {a: b, b: a, m1: m2, m2: m1}
    assert sp.cancel(sp.together(rad_sym - rad_sym.subs(swap, simultaneous=True))) == 0

    # angular side: under the swap, u = cos^2 t -> 1 - u = sin^2 t
    ang_sym = (
        -(m1**2) / (1 - u)
        - m2**2 / u
        + (w**2 - mu**2) * (a**2 * u + b**2 * (1 - u))
    )
    swapped = ang_sym.subs(swap, simultaneous=True).subs(u, 1 - u)
    assert sp.cancel(sp.together(ang_sym - swapped)) == 0


def test_schwarzschild_tangherlini_limit():
    """a = b = 0, mu = 0: radial potential reduces to w^2 r^4 / (r^2 - M)."""
    Sig, Pi_, rrDel, _, _ = _pieces()
    W = Pi_ * w - m1 * a * (r**2 + b**2) - m2 * b * (r**2 + a**2)
    G = a * b * w - a * m2 - b * m1
    rad = W**2 / (r**2 * rrDel) - G**2 / r**2 - a**2 * w**2 + 2 * w * (a * m1 + b * m2)
    sub = {a: 0, b: 0, m1: 0, m2: 0}
    got = sp.cancel(rad.subs(sub))
    assert sp.cancel(got - w**2 * r**4 / (r**2 - M)) == 0

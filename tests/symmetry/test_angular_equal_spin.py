"""Exactness of the equal-spin angular degeneracy, and its unfolding in delta.

These tests establish, numerically, the statement that drives the project:
on ``delta = 0`` the shifted angular eigenvalue ``Ahat`` equals ``l(l+2)``
*exactly* for every ``(n, m1, m2)`` in the ``S^3`` multiplet, for arbitrary
``s`` and ``mu``.
"""

import pytest

from mp5d.angular import (
    angular_eigenvalue,
    angular_spectrum,
    angular_spectrum_fd,
    l_of_n,
    s3_multiplet,
    spheroidicity,
)
from mp5d.geometry import MPGeometry


@pytest.mark.parametrize("s", [0.0, 0.15, 0.3, 0.45])
@pytest.mark.parametrize("mu", [0.0, 0.2, 0.7])
def test_spheroidicity_vanishes_identically_on_equal_spin(s, mu):
    g = MPGeometry.from_sdelta(s=s, delta=0.0)
    omega = 0.83 - 0.19j
    assert spheroidicity(omega, mu, g.a, g.b) == 0


@pytest.mark.parametrize("l", [0, 1, 2, 3, 4])
def test_s3_multiplet_dimension(l):
    """The equal-spin degenerate multiplet has dimension (l+1)^2."""
    assert len(s3_multiplet(l)) == (l + 1) ** 2


@pytest.mark.parametrize("l", [0, 1, 2, 3])
def test_equal_spin_eigenvalue_is_exactly_l_l_plus_2(l):
    """Every member of the multiplet shares Ahat = l(l+2) at c2 = 0, exactly."""
    for n, m1, m2 in s3_multiplet(l):
        Ahat = angular_eigenvalue(m1, m2, n, c2=0.0, N=30)
        assert Ahat == pytest.approx(l * (l + 2), abs=1e-10)
        assert l_of_n(n, m1, m2) == l


def test_degeneracy_is_lifted_only_by_the_product_s_delta():
    """c2 = 4 s delta (w^2 - mu^2): both s = 0 and delta = 0 restore degeneracy."""
    omega, mu = 0.7 - 0.1j, 0.3
    a_l = []
    for s, d in [(0.0, 0.2), (0.3, 0.0), (0.3, 0.2)]:
        g = MPGeometry.from_sdelta(s, d)
        c2 = spheroidicity(omega, mu, g.a, g.b)
        # two distinct members of the l = 2 multiplet
        A1 = angular_eigenvalue(2, 0, 0, c2, N=40)
        A2 = angular_eigenvalue(1, 1, 0, c2, N=40)
        a_l.append(abs(A1 - A2))
    assert a_l[0] == pytest.approx(0.0, abs=1e-10)  # s = 0
    assert a_l[1] == pytest.approx(0.0, abs=1e-10)  # delta = 0
    assert a_l[2] > 1e-3  # generic point: split


@pytest.mark.parametrize(
    "m1,m2,c2",
    [
        (0, 0, 1.3),
        (1, 0, -2.1),
        (2, 1, 0.9 + 0.4j),
        (1, 1, -1.7 + 0.8j),
        (3, 2, 2.5 - 1.1j),
    ],
)
def test_spectral_matches_independent_finite_difference(m1, m2, c2):
    """Tridiagonal Jacobi construction vs. a plain FD discretization."""
    spec, _ = angular_spectrum(m1, m2, c2, N=60)
    fd = angular_spectrum_fd(m1, m2, c2, N=1500)
    for k in range(3):  # lowest three angular branches
        best = min(abs(fd - spec[k]))
        assert best < 5e-3, f"branch {k}: spectral {spec[k]}, closest FD off by {best}"


@pytest.mark.parametrize("c2", [0.8, -1.5, 1.1 - 0.6j])
def test_angular_exchange_symmetry(c2):
    """(a,m1)<->(b,m2) flips the sign of c2 and swaps m1,m2.

    Sending ``u -> 1 - u`` in the angular equation maps
    ``(m1, m2, c2, Ahat) -> (m2, m1, -c2, Ahat + c2)``, hence the exact relation

        Ahat(m1, m2, c2) = Ahat(m2, m1, -c2) - c2.
    """
    for n in range(3):
        A = angular_eigenvalue(m1=2, m2=1, n=n, c2=c2, N=50)
        B = angular_eigenvalue(m1=1, m2=2, n=n, c2=-c2, N=50) - c2
        assert A == pytest.approx(B, rel=1e-9, abs=1e-9)


def test_convergence_in_truncation():
    """Geometric convergence of Ahat with basis size.

    ``c2`` is taken large enough that the truncation error is still resolvable
    at small ``N``; at ``|c2| ~ 1`` the expansion is already at machine
    precision by ``N ~ 12``, which would make the monotonicity check vacuous.
    """
    # Physical regime: |c2| = O(1).  Converged to machine precision by N ~ 10,
    # so the meaningful assertion is stability of the value, not a decay rate.
    c2 = 2.0 - 1.0j
    ref = angular_eigenvalue(1, 1, 1, c2, N=200, n_steps=8)
    for N in (10, 16, 24, 60):
        assert abs(angular_eigenvalue(1, 1, 1, c2, N=N, n_steps=8) - ref) < 1e-12

    # Stress regime: |c2| large enough that truncation error is resolvable, so
    # the geometric decay itself can be checked.
    c2 = 4000.0 - 1500.0j
    ref = angular_eigenvalue(1, 1, 1, c2, N=400, n_steps=40)
    errs = [abs(angular_eigenvalue(1, 1, 1, c2, N=N, n_steps=40) - ref) for N in (12, 24, 36, 60)]
    assert errs[0] > 1e-6, f"stress c2 too small to resolve truncation: {errs}"
    assert errs[-1] < 1e-8, f"not converged: {errs}"
    assert all(x > y for x, y in zip(errs, errs[1:], strict=False)), f"non-monotone: {errs}"

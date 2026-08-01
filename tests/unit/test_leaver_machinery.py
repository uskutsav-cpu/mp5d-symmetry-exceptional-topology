"""Tests for the pieces of the radial solvers that ARE working.

The radial eigenvalue solvers are not yet functional (see
``docs/FAILED_APPROACHES.md`` entries 3 and 5).  These tests pin the components
that are correct, so the working parts do not regress while the root-finding
condition is fixed, and they pin the *diagnoses* so the known failure modes are
detected if anyone reintroduces them.
"""

import numpy as np
import pytest

from mp5d.geometry import MPGeometry
from mp5d.radial.leaver import LeaverProblem, polynomial_coefficients


def test_fft_polynomial_recovery_is_exact():
    rng = np.random.default_rng(20260801)
    true = rng.normal(size=9) + 1j * rng.normal(size=9)

    def f(u):
        return sum(c * u**j for j, c in enumerate(true))

    got = polynomial_coefficients(f, 16)
    assert np.abs(got[:9] - true).max() < 1e-12
    assert np.abs(got[9:]).max() < 1e-12


def test_fft_recovery_rejects_an_uncleared_pole():
    """The tail check must catch a function that is not a polynomial."""
    with pytest.raises((ValueError, FloatingPointError)):
        polynomial_coefficients(lambda u: 1.0 / (u - 1.05), 16)


def test_fft_recovery_is_not_index_reversed():
    """Regression: numpy's ifft convention returns coefficients reversed."""
    got = polynomial_coefficients(lambda u: 3.0 + 0.0 * u, 8)
    assert abs(got[0] - 3.0) < 1e-12
    assert np.abs(got[1:]).max() < 1e-12


@pytest.mark.parametrize(
    "a,b", [(0.0, 0.0), (0.3, 0.0), (0.25, 0.25), (0.35, 0.12)]
)
def test_leaver_ode_coefficients_are_polynomial_after_clearing(a, b):
    """The clearing factor removes every pole, for single, equal and unequal spins.

    This is the self-check built into the FFT recovery: if a pole survived, the
    recovered tail would not vanish and ``polynomial_coefficients`` would raise.
    """
    geo = MPGeometry(a=a, b=b, M=1.0)
    p = LeaverProblem(geo, 0.0, 0, 0, 0, depth=30)
    omega = 0.9 - 0.9j
    A, B, C = p.poly_ABC(omega, p.Lambda_of(omega), degree_bound=32)
    assert len(A) > 0 and len(B) > 0 and len(C) > 0
    assert np.isfinite(np.abs(A)).all()


def test_leaver_degrees_are_stable_under_degree_bound():
    geo = MPGeometry(a=0.0, b=0.0, M=1.0)
    p = LeaverProblem(geo, 0.0, 0, 0, 0, depth=30)
    omega = 0.9 - 0.9j
    Lam = p.Lambda_of(omega)

    def degrees(db):
        out = []
        for arr in p.poly_ABC(omega, Lam, degree_bound=db):
            big = np.abs(arr).max()
            nz = np.nonzero(np.abs(arr) > 1e-11 * big)[0]
            out.append(int(nz[-1]) if len(nz) else -1)
        return tuple(out)

    assert degrees(20) == degrees(32)


def test_common_u_power_is_stripped():
    """Regression: leaving the clearing factor's u^k makes the Hill matrix
    singular for every omega (observed as a uniform rel sigma_min ~ 1e-32)."""
    geo = MPGeometry(a=0.0, b=0.0, M=1.0)
    p = LeaverProblem(geo, 0.0, 0, 0, 0, depth=30)
    omega = 0.9 - 0.9j
    A, B, C = p.poly_ABC(omega, p.Lambda_of(omega))
    # after stripping, at least one of A, B, C must have a nonzero constant term
    consts = [abs(arr[0]) / max(np.abs(arr).max(), 1e-300) for arr in (A, B, C)]
    assert max(consts) > 1e-8, "common power of u was not stripped"

    M = p.hill_matrix(omega)
    rows = np.abs(M).max(axis=1)
    assert rows.min() > 0, "Hill matrix still has an identically zero row"


@pytest.mark.parametrize("mu", [0.0, 0.3])
def test_lambda_is_recomputed_not_frozen(mu):
    """Lambda must track omega, since c2 = (w^2 - mu^2)(a^2 - b^2)."""
    geo = MPGeometry(a=0.35, b=0.12, M=1.0)
    p = LeaverProblem(geo, mu, 1, 1, 2, depth=20)
    L1 = p.Lambda_of(0.6 - 0.2j)
    L2 = p.Lambda_of(0.9 - 0.4j)
    assert abs(L1 - L2) > 1e-6

    # ...and on the equal-spin surface it must NOT depend on omega beyond the
    # exact b^2 shift, because c2 vanishes identically there.
    eq = MPGeometry.from_sdelta(s=0.3, delta=0.0)
    q = LeaverProblem(eq, mu, 1, 1, 2, depth=20)
    for om in (0.6 - 0.2j, 0.9 - 0.4j):
        Ahat = q.Lambda_of(om) + (om**2 - mu**2) * eq.b**2
        assert abs(Ahat - 2 * 4) < 1e-9  # l(l+2) with l = 2

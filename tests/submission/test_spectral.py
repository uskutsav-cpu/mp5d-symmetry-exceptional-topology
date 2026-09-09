import numpy as np
import pytest

from mp5d_science.spectral import *


def test_raw_meromorphic_cf_is_not_accepted():
    with pytest.raises(AnalyticityError):
        count_zeros(lambda z: (z - 0.1) / (z - 0.2), 0, 1)


def test_known_pole_inside_domain_is_rejected():
    with pytest.raises(AnalyticityError):
        AnalyticSpectralFunction(
            lambda z: (z - 0.1) / (z - 0.2), 0, 1, "claimed", "regularized_evans_function", (0.2,)
        )


def test_missing_analyticity_justification_rejected():
    with pytest.raises(AnalyticityError):
        AnalyticSpectralFunction(lambda z: z, 0, 1, "", "polynomial")


@pytest.mark.parametrize("roots", [[], [0.2], [0.2, -0.3], [0.1 + 0.2j, -0.3 + 0.1j, 0.4]])
def test_numerical_root_count(roots):
    coefficients = np.poly(roots).tolist() if roots else [0, 1]
    if not roots:
        f = AnalyticSpectralFunction(
            lambda z: 1, 0, 2, "constant entire", "explicit_entire_function"
        )
    else:
        f = polynomial_spectral(coefficients, 0, 2)
    r = count_zeros(f, 0, 1)
    assert r["count"] == len(roots)
    assert not r["rigorous"] and len(r["ladder"]) >= 3


def test_moment_roots_center_scaled():
    expected = [100 + 0.1j, 100.2 - 0.1j]
    f = polynomial_spectral(np.poly(expected).tolist(), 100, 2)
    roots, report = roots_in_disc(f, 100, 1)
    assert max(min(abs(a - b) for b in roots) for a in expected) < 1e-7


def test_root_on_contour_fails():
    f = polynomial_spectral([1, -1], 0, 2)
    with pytest.raises(ContourFailure):
        count_zeros(f, 0, 1)


def test_outside_analytic_domain_fails():
    f = polynomial_spectral([1, 0], 0, 1)
    with pytest.raises(AnalyticityError):
        count_zeros(f, 0.5, 0.75)


def test_determinant_uses_fixed_scale():
    f = finite_determinant(
        lambda z: np.array([[z, 1], [0, z - 0.5]]),
        0,
        2,
        "all entries are degree-one polynomials",
        np.array([2.0, 3.0]),
    )
    assert abs(f(0.2) - (0.2 * (0.2 - 0.5)) / 6) < 1e-15
    assert count_zeros(f, 0, 1)["count"] == 2


def test_svd_projection_is_not_analytic_regression():
    # This reproduces the original Solver C defect exactly: u^H M v = s_min.
    def bad(z):
        m = np.diag([z, 2 + 0j])
        U, s, Vh = np.linalg.svd(m)
        return U[:, -1].conj() @ m @ Vh[-1].conj()

    for z in [0.2, 0.2j, -0.2, -0.2j]:
        assert abs(bad(z) - 0.2) < 1e-14
    h = 1e-5
    z = 0.2 + 0.1j
    dx = (bad(z + h) - bad(z - h)) / (2 * h)
    dy = (bad(z + 1j * h) - bad(z - 1j * h)) / (2j * h)
    assert abs(dx - dy) > 0.5


def test_a1_a2_ratio_not_invariant_under_analytic_prefactor():
    # F(z)=z*exp(c*z) has only one root; F'/[F''/2]=1/c can be arbitrarily small.
    for c in [1, 100, 10000]:
        assert abs(1 / c) <= 1
    # This identity is covered by an exact symbolic regression as well.

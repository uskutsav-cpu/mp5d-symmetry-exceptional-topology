import pytest
import sympy as sp

from mp5d_science.certification import *


def test_exact_interval_arithmetic():
    a = Interval("1/3", "2/3")
    b = Interval("-2", "1")
    assert (a + b) == Interval("-5/3", "5/3")
    assert (a * b) == Interval("-4/3", "2/3")
    assert a.reciprocal() == Interval("3/2", "3")


def test_interval_zero_division():
    with pytest.raises(ZeroDivisionError):
        Interval("-1", "1").reciprocal()


def test_no_float_coefficients_in_rigorous_path():
    with pytest.raises(TypeError):
        krawczyk_polynomial([[1.0, 0], [0, 0]], ["0", "0"], "0.1", problem_id="bad")


def test_real_sqrt_two_certified():
    result = krawczyk_polynomial(
        [["1", "0"], ["0", "0"], ["-2", "0"]],
        ["1.4142135623730950488", "0"],
        "0.00000001",
        problem_id="sqrt2",
    )
    assert result["status"] == "CERTIFIED_UNIQUE_FINITE_ROOT"
    assert not result["continuum_certificate"]


def test_complex_root_certified():
    result = krawczyk_polynomial(
        [["1", "0"], ["0", "0"], ["1", "0"]], ["0", "1"], "0.000001", problem_id="imaginary-unit"
    )
    assert result["status"] == "CERTIFIED_UNIQUE_FINITE_ROOT"


def test_multiple_root_not_certified_unique_by_singular_newton():
    result = krawczyk_polynomial(
        [["1", "0"], ["0", "0"], ["0", "0"]], ["0", "0"], "0.1", problem_id="ep"
    )
    assert result["status"] == "INCONCLUSIVE"


def test_wrong_center_does_not_certify():
    result = krawczyk_polynomial(
        [["1", "0"], ["0", "0"], ["-2", "0"]], ["2", "0"], "0.01", problem_id="wrong"
    )
    assert result["status"] == "INCONCLUSIVE"


def test_radial_finite_hessenberg_determinant_against_dense_sympy():
    from mp5d_science.radial_polynomial import coefficient_formula

    size = 4
    w = sp.Symbol("w")
    A, B, C, _ = coefficient_formula(
        sp.Integer(0),
        sp.Integer(0),
        sp.Integer(0),
        0,
        0,
        w,
        sp.Integer(0),
        sp.Integer(1),
        sp.Integer(0),
        w,
        sp.I,
    )

    def get(p, i):
        return p[i] if 0 <= i < len(p) else 0

    matrix = sp.Matrix(
        size,
        size,
        lambda n, j: get(A, n - j + 2) * j * (j - 1) + get(B, n - j + 1) * j + get(C, n - j),
    )
    dense = sp.Poly(matrix.det(), w).monic()
    data = radial_finite_polynomial(0, size)
    result = sp.Poly.from_list(
        [sp.Rational(a) + sp.I * sp.Rational(b) for a, b in data["coefficients"]], w
    )
    assert sp.expand(dense.as_expr() - result.as_expr()) == 0

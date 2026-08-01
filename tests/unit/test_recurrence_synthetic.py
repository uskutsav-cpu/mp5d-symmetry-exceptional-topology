"""Validate the recurrence machinery on synthetic ODEs, with no physics involved.

If these pass, the extraction / reduction / continued-fraction code is correct
as *numerics*.  Any later disagreement with a black-hole benchmark is then a
physics or convention problem, not a bug in this layer.
"""

import numpy as np
import pytest

from mp5d.radial.recurrence import (
    ReductionDiagnostics,
    check_frobenius_form,
    continued_fraction,
    continued_fraction_lentz,
    recurrence_row,
    recurrence_width,
    reduce_to_three_term,
)


def poly_mul(p, q):
    out = np.zeros(len(p) + len(q) - 1, dtype=complex)
    for i, a in enumerate(p):
        out[i : i + len(q)] += a * np.asarray(q)
    return out


# A synthetic ODE with a genuine 3-term recurrence.
# A(x) = x - x^2, B(x) = b0 + b1 x, C(x) = c0 + c1 x  ->  A_0 = 0 as required.
A3 = np.array([0.0, 1.0, -1.0], dtype=complex)
B3 = np.array([0.7, -1.3], dtype=complex)
C3 = np.array([0.25, 0.4], dtype=complex)


def test_width_and_frobenius_guard():
    assert recurrence_width(A3, B3, C3) == 3
    check_frobenius_form(A3)  # must not raise
    with pytest.raises(ValueError):
        check_frobenius_form(np.array([1.0, 1.0, 1.0], dtype=complex))


def test_recurrence_row_matches_direct_series_substitution():
    """Independent check: expand A y'' + B y' + C y for an explicit series."""
    rng = np.random.default_rng(11)
    a = rng.normal(size=14) + 1j * rng.normal(size=14)

    # direct: coefficient of x^k in A y'' + B y' + C y
    def direct(k):
        tot = 0.0 + 0.0j
        for j, Aj in enumerate(A3):
            m = k - j + 2
            if 0 <= m < len(a):
                tot += Aj * m * (m - 1) * a[m]
        for j, Bj in enumerate(B3):
            m = k - j + 1
            if 0 <= m < len(a):
                tot += Bj * m * a[m]
        for j, Cj in enumerate(C3):
            m = k - j
            if 0 <= m < len(a):
                tot += Cj * a[m]
        return tot

    w = recurrence_width(A3, B3, C3)
    for k in range(2, 9):
        row = recurrence_row(k, A3, B3, C3, w)
        viarow = sum(
            row[i] * a[k + 1 - i] for i in range(w) if 0 <= k + 1 - i < len(a)
        )
        assert abs(viarow - direct(k)) < 1e-11


def test_reduction_is_identity_on_a_three_term_recurrence():
    alpha, beta, gamma = reduce_to_three_term(A3, B3, C3, depth=20)
    w = recurrence_width(A3, B3, C3)
    for n in range(1, 20):
        row = recurrence_row(n, A3, B3, C3, w)
        assert abs(alpha[n] - row[0]) < 1e-12
        assert abs(beta[n] - row[1]) < 1e-12
        assert abs(gamma[n] - row[2]) < 1e-12


@pytest.mark.parametrize(
    "g",
    [
        np.array([1.0, 2.0], dtype=complex),
        np.array([1.0, -0.5, 0.3], dtype=complex),
        np.array([2.0, 0.1, -0.7, 0.4], dtype=complex),
    ],
)
def test_widening_preserves_the_forward_solution(g):
    """Multiplying the ODE by g(x) changes nothing physically.

    The raw recurrence gets wider, but the sequence generated forward from
    ``a_0 = 1`` must be identical.  (Note the reduced three-term rows are NOT
    proportional between the two: with ``A_0 = 0`` the raw recurrence already
    determines ``a_{n+1}`` uniquely, so its solution space is one-dimensional,
    while a three-term recurrence has a two-dimensional one.  The reduction
    embeds the former in the latter, and the extra direction depends on the
    elimination history.  Only the forward solution is invariant.)
    """
    Aw, Bw, Cw = poly_mul(g, A3), poly_mul(g, B3), poly_mul(g, C3)
    assert recurrence_width(Aw, Bw, Cw) > 3

    def forward(A, B, C, depth=30):
        w = recurrence_width(A, B, C)
        a = np.zeros(depth + 2, dtype=complex)
        a[0] = 1.0
        for n in range(depth):
            row = recurrence_row(n, A, B, C, w)
            acc = sum(
                row[i] * a[n + 1 - i] for i in range(1, w) if 0 <= n + 1 - i < len(a)
            )
            a[n + 1] = -acc / row[0]
        return a

    a_raw = forward(A3, B3, C3)
    a_wide = forward(Aw, Bw, Cw)
    assert np.abs(a_raw - a_wide).max() < 1e-9


def test_generated_sequence_satisfies_both_recurrences():
    """Solve the raw N-term recurrence forward; the reduced rows must hold too."""
    g = np.array([1.0, -0.5, 0.3], dtype=complex)
    Aw, Bw, Cw = poly_mul(g, A3), poly_mul(g, B3), poly_mul(g, C3)
    w = recurrence_width(Aw, Bw, Cw)
    depth = 40

    a = np.zeros(depth + 2, dtype=complex)
    a[0] = 1.0
    for n in range(depth):
        row = recurrence_row(n, Aw, Bw, Cw, w)
        acc = sum(
            row[i] * a[n + 1 - i] for i in range(1, w) if 0 <= n + 1 - i < len(a)
        )
        a[n + 1] = -acc / row[0]

    alpha, beta, gamma = reduce_to_three_term(Aw, Bw, Cw, depth=depth)
    for n in range(1, depth - 1):
        res = alpha[n] * a[n + 1] + beta[n] * a[n] + gamma[n] * a[n - 1]
        scale = max(
            abs(alpha[n] * a[n + 1]), abs(beta[n] * a[n]), abs(gamma[n] * a[n - 1])
        )
        assert abs(res) < 1e-8 * max(scale, 1e-300), f"3-term fails at n={n}"


def test_lentz_agrees_with_backward_recurrence():
    alpha, beta, gamma = reduce_to_three_term(A3, B3, C3, depth=300)
    back = continued_fraction(alpha, beta, gamma, 300)
    lentz = continued_fraction_lentz(alpha, beta, gamma, 300)
    assert abs(back - lentz) < 1e-8 * max(abs(back), 1.0)


def test_reduction_diagnostics_are_recorded():
    g = np.array([2.0, 0.1, -0.7, 0.4], dtype=complex)
    Aw, Bw, Cw = poly_mul(g, A3), poly_mul(g, B3), poly_mul(g, C3)
    diag = ReductionDiagnostics()
    reduce_to_three_term(Aw, Bw, Cw, depth=30, diagnostics=diag)
    assert diag.n_steps > 0
    assert np.isfinite(diag.min_pivot)

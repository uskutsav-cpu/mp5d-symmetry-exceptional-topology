"""Frobenius recurrence machinery: extraction, N-term -> 3-term reduction, CF.

This module is deliberately **physics-free**.  Everything here operates on the
polynomial coefficients of a linear ODE

    A(x) y'' + B(x) y' + C(x) y = 0,      y = sum_{n>=0} a_n x^n,

so it can be validated against synthetic recurrences with known solutions,
independently of any black-hole problem.  ``mp5d.radial.qnm`` supplies the
Myers-Perry ``A, B, C``.

Recurrence
----------
Matching the coefficient of ``x^k`` gives

    sum_j [ A_j (k-j+2)(k-j+1) a_{k-j+2} + B_j (k-j+1) a_{k-j+1} + C_j a_{k-j} ] = 0

Re-indexing by ``i`` such that the unknown is ``a_{k+1-i}``:

    c[k, i] = A_{i+1} (k+1-i)(k-i) + B_i (k+1-i) + C_{i-1}

so equation ``k`` involves ``a_{k+1}`` down to ``a_{k+1-p}`` where
``p = max(deg A - 1, deg B, deg C + 1)``.  Note ``A_0`` multiplies ``a_{k+2}``;
for a Frobenius expansion at a regular singular point ``A_0 = 0``, so the
leading unknown is ``a_{k+1}`` and the recurrence is genuinely forward-solving.

Reduction
---------
:func:`reduce_to_three_term` implements the general Gaussian elimination that
turns a ``(p+1)``-term recurrence into a three-term one.  The mathematics
follows Leaver's classical elimination as generalized by Karikos, Saes, Wagle
and Yunes, *Beyond Three Terms: Continued Fractions for Rotating Black Holes in
Modified Gravity*, arXiv:2604.18680 (CC BY 4.0).  Reimplemented from the
mathematics with attribution; no source was copied.

The idea: process ``n`` upward.  Once the row at index ``m`` has already been
reduced to ``alpha_m a_{m+1} + beta_m a_m + gamma_m a_{m-1} = 0``, it expresses
``a_{m-1}`` in terms of ``a_{m+1}`` and ``a_m``.  Substituting that into a wider
row at index ``n`` removes its lowest entry and raises the row's floor by one.
Repeating until only three entries survive is exact and preserves the solution
space of the infinite system.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "check_frobenius_form",
    "recurrence_row",
    "recurrence_width",
    "reduce_to_three_term",
    "continued_fraction",
    "continued_fraction_inverted",
    "cf_residual",
    "continued_fraction_lentz",
    "ReductionDiagnostics",
]


class ReductionDiagnostics:
    """Records pivot behaviour during the Gaussian reduction."""

    def __init__(self) -> None:
        self.min_pivot: float = float("inf")
        self.min_pivot_index: int = -1
        self.n_steps: int = 0

    def note(self, n: int, pivot: complex) -> None:
        self.n_steps += 1
        p = abs(pivot)
        if p < self.min_pivot:
            self.min_pivot = p
            self.min_pivot_index = n

    def __repr__(self) -> str:
        return (
            f"ReductionDiagnostics(steps={self.n_steps}, "
            f"min|pivot|={self.min_pivot:.3e} at n={self.min_pivot_index})"
        )


def _trim(arr, rel: float = 1e-13):
    """Drop numerically-zero trailing coefficients so degrees are honest."""
    a = np.asarray(arr, dtype=complex)
    big = np.abs(a).max() if a.size else 0.0
    if big == 0:
        return a[:1]
    nz = np.nonzero(np.abs(a) > rel * big)[0]
    return a[: nz[-1] + 1]


def check_frobenius_form(A, rel: float = 1e-10) -> None:
    """Require ``A_0 = 0``.

    If ``A(0) != 0`` the coefficient of ``x^k`` also contains ``a_{k+2}``, and a
    three-term reduction of the form used here is not applicable.  ``A_0 = 0`` is
    exactly the statement that ``x = 0`` is a regular singular point whose
    indicial exponent has already been peeled off, which is the situation this
    module is built for.  Failing loudly beats silently dropping a term.
    """
    A = np.asarray(A, dtype=complex)
    big = np.abs(A).max()
    if big > 0 and abs(A[0]) > rel * big:
        raise ValueError(
            f"A(0) = {A[0]:.3e} is not negligible relative to max|A| = {big:.3e}; "
            "the expansion point is not a peeled regular singular point and the "
            "recurrence would also involve a_{k+2}"
        )


def recurrence_width(A, B, C) -> int:
    """``p+1``: number of terms in the raw recurrence."""
    A, B, C = _trim(A), _trim(B), _trim(C)
    p = max(len(A) - 2, len(B) - 1, len(C))
    return max(p, 1) + 1


def recurrence_row(k: int, A, B, C, width: int):
    """Coefficients ``c[i]`` of ``a_{k+1-i}`` in equation ``k``, ``i = 0..width-1``."""

    def g(arr, j):
        return arr[j] if 0 <= j < len(arr) else 0.0

    out = np.zeros(width, dtype=complex)
    for i in range(width):
        m = k + 1 - i  # index of the unknown a_m
        out[i] = g(A, i + 1) * m * (m - 1) + g(B, i) * m + g(C, i - 1)
    return out


def reduce_to_three_term(A, B, C, depth: int, diagnostics: ReductionDiagnostics | None = None):
    """Reduce the raw recurrence to ``alpha_n a_{n+1} + beta_n a_n + gamma_n a_{n-1} = 0``.

    Returns ``(alpha, beta, gamma)`` arrays of length ``depth``.  ``gamma[0]`` is
    unused because ``a_{-1} = 0``.

    Exactness: at every step the row being modified is replaced by an exact
    linear combination of itself and an already-reduced row, so the solution
    space of the infinite system is preserved.
    """
    check_frobenius_form(A)
    width = recurrence_width(A, B, C)
    alpha = np.zeros(depth, dtype=complex)
    beta = np.zeros(depth, dtype=complex)
    gamma = np.zeros(depth, dtype=complex)

    for n in range(depth):
        row = recurrence_row(n, A, B, C, width)
        # row[i] multiplies a_{n+1-i}; entries with n+1-i < 0 are unreachable
        for i in range(width):
            if n + 1 - i < 0:
                row[i] = 0.0

        # eliminate everything below a_{n-1}, i.e. indices i >= 3
        i = width - 1
        while i >= 3:
            if row[i] != 0.0:
                m = n + 1 - i  # the index we are eliminating: a_m
                # use the reduced row at index m+1:
                #   alpha_{m+1} a_{m+2} + beta_{m+1} a_{m+1} + gamma_{m+1} a_m = 0
                #   => a_m = -(alpha_{m+1} a_{m+2} + beta_{m+1} a_{m+1}) / gamma_{m+1}
                j = m + 1
                if j < 0 or m < 0:
                    row[i] = 0.0
                    i -= 1
                    continue
                piv = gamma[j]
                if diagnostics is not None:
                    diagnostics.note(j, piv)
                if piv == 0:
                    raise ZeroDivisionError(f"zero pivot gamma[{j}] during reduction at n={n}")
                f = row[i] / piv
                # a_m -> -(alpha_j a_{m+2} + beta_j a_{m+1})/gamma_j
                # a_{m+2} sits at position i-2, a_{m+1} at position i-1
                row[i - 2] -= f * alpha[j]
                row[i - 1] -= f * beta[j]
                row[i] = 0.0
            i -= 1

        alpha[n] = row[0]
        beta[n] = row[1]
        gamma[n] = row[2] if width > 2 else 0.0

    return alpha, beta, gamma


def continued_fraction(alpha, beta, gamma, depth: int | None = None, tiny: float = 1e-300):
    """Leaver's condition, evaluated bottom-up (backward recurrence).

    Returns ``beta_0 - alpha_0 gamma_1 / (beta_1 - alpha_1 gamma_2 / (...))``.
    A quasinormal mode is a zero of this function.
    """
    n = len(alpha) if depth is None else min(depth, len(alpha))
    frac = 0.0 + 0.0j
    for k in range(n - 1, 0, -1):
        den = beta[k] - frac
        if abs(den) < tiny:
            den = tiny
        frac = alpha[k - 1] * gamma[k] / den
    return beta[0] - frac


def continued_fraction_lentz(alpha, beta, gamma, depth: int, tiny: float = 1e-300):
    """Modified Lentz evaluation of the same continued fraction (top-down).

    Provided as an independent evaluation path: Lentz and the backward
    recurrence must agree, which catches indexing and convergence errors.
    """
    f = beta[0]
    if abs(f) < tiny:
        f = tiny
    Cv = f
    Dv = 0.0 + 0.0j
    for k in range(1, depth):
        an = -alpha[k - 1] * gamma[k]
        bn = beta[k]
        Dv = bn + an * Dv
        if abs(Dv) < tiny:
            Dv = tiny
        Cv = bn + an / Cv
        if abs(Cv) < tiny:
            Cv = tiny
        Dv = 1.0 / Dv
        delta = Cv * Dv
        f = f * delta
        if abs(delta - 1.0) < 1e-16:
            break
    return f


def continued_fraction_inverted(
    alpha, beta, gamma, inversion: int, depth: int | None = None, tiny: float = 1e-300
):
    """The ``inversion``-th inversion of Leaver's condition.

    Algebraically equivalent to :func:`continued_fraction` (same zero set), but
    numerically it isolates the ``inversion``-th root, which is what makes
    overtone-specific root finding stable.  ``inversion = 0`` reduces to the
    plain condition.

    Two practical reasons this is needed rather than optional:

    * a plain Newton solve on the ``0``-th form wanders between branches when
      the seed is not already very close to the intended root;
    * when the Frobenius series *terminates* (a polynomial solution) the plain
      form degenerates to ``0/0``, because ``gamma`` and ``beta`` both vanish at
      the truncation index.  The inverted form evaluates the finite upward part
      separately and stays well defined.
    """
    n = len(alpha) if depth is None else min(depth, len(alpha))
    k = int(inversion)
    if k <= 0:
        return continued_fraction(alpha, beta, gamma, n, tiny)

    # downward (infinite) tail:  alpha_k gamma_{k+1} / (beta_{k+1} - ...)
    frac = 0.0 + 0.0j
    for j in range(n - 1, k, -1):
        den = beta[j] - frac
        if abs(den) < tiny:
            den = tiny
        frac = alpha[j - 1] * gamma[j] / den
    down = frac

    # upward (finite) part: continued fraction running from index 0 up to k
    up = 0.0 + 0.0j
    for j in range(0, k):
        den = beta[j] - up
        if abs(den) < tiny:
            den = tiny
        up = alpha[j] * gamma[j + 1] / den

    return beta[k] - down - up


def cf_residual(alpha, beta, gamma, depth: int) -> float:
    """|CF(depth) - CF(depth//2)|, a truncation-stability estimate."""
    a = continued_fraction(alpha, beta, gamma, depth)
    b = continued_fraction(alpha, beta, gamma, max(depth // 2, 2))
    return float(abs(a - b))

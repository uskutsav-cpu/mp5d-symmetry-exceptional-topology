"""Angular (5D spheroidal) operator for the massive scalar on 5D Myers-Perry.

Separated angular equation (derivation and conventions: ``docs/CONVENTIONS.md``,
regression-tested by ``tests/unit/test_separation.py``)::

    1/(sin t cos t) d/dt ( sin t cos t dS/dt )
      + [ Lambda + (w^2 - mu^2)(a^2 cos^2 t + b^2 sin^2 t)
          - m1^2 / sin^2 t - m2^2 / cos^2 t ] S = 0

In ``u = cos^2 theta`` this is

    4 d/du ( u(1-u) dS/du )
      + [ Ahat + c2 * u - m1^2/(1-u) - m2^2/u ] S = 0

with the two quantities that control everything downstream::

    Ahat = Lambda + (w^2 - mu^2) b^2          (shifted angular eigenvalue)
    c2   = (w^2 - mu^2) (a^2 - b^2)           (spheroidicity)

Because ``a^2 - b^2 = 4 s delta``, the spheroidicity is *exactly*

    c2 = 4 s delta (w^2 - mu^2),

so ``c2`` vanishes identically on the equal-spin surface ``delta = 0`` for every
``s`` and every ``mu``.  There the operator is the round ``S^3`` Laplacian and

    Ahat = l (l + 2),   l = 2 n + |m1| + |m2|,   n = 0, 1, 2, ...

exactly, with the full ``SO(4)`` degeneracy.  This is an exact statement, not a
perturbative one, and it is the origin of the enhanced-symmetry structure this
project studies.

Discretization
--------------
Substituting ``S = u^(alpha/2) (1-u)^(beta/2) f(u)`` with ``alpha = |m2|``,
``beta = |m1|`` turns the ``c2 = 0`` problem into Jacobi's equation, so the
orthonormal Jacobi polynomials ``p_n^(alpha, beta)(x)``, ``x = 1 - 2u``, are the
exact eigenbasis at ``delta = 0``.  Multiplication by ``u = (1 - x)/2`` is
tridiagonal in that basis, hence the angular operator is an *exactly
tridiagonal, complex-symmetric* matrix

    H = D - (c2/2) (I - X),    D_nn = l_n (l_n + 2),  l_n = 2n + alpha + beta,

where ``X`` is the Jacobi matrix of ``p_n^(alpha, beta)``.  Truncation is the
only approximation, and it converges geometrically in ``N``.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "jacobi_recurrence",
    "angular_matrix",
    "angular_spectrum",
    "angular_eigenvalue",
    "l_of_n",
    "s3_multiplet",
    "spheroidicity",
]


def spheroidicity(omega: complex, mu: float, a: float, b: float) -> complex:
    """``c2 = (w^2 - mu^2)(a^2 - b^2) = 4 s delta (w^2 - mu^2)``."""
    return (omega**2 - mu**2) * (a**2 - b**2)


def l_of_n(n: int, m1: int, m2: int) -> int:
    """``S^3`` angular momentum of the ``n``-th radial-node angular mode."""
    return 2 * n + abs(m1) + abs(m2)


def s3_multiplet(l: int) -> list[tuple[int, int, int]]:
    """All ``(n, m1, m2)`` with ``2n + |m1| + |m2| = l``.

    These are exactly the states that are degenerate in ``Ahat`` on the
    equal-spin surface.  The list has ``(l+1)**2`` entries, the dimension of the
    ``SO(4)`` scalar harmonic representation on ``S^3``.
    """
    out = []
    for m1 in range(-l, l + 1):
        for m2 in range(-l, l + 1):
            rest = l - abs(m1) - abs(m2)
            if rest >= 0 and rest % 2 == 0:
                out.append((rest // 2, m1, m2))
    return out


def jacobi_recurrence(alpha: int, beta: int, N: int):
    """Diagonal/off-diagonal of the Jacobi matrix of orthonormal ``p_n^(a,b)``.

    Returns ``(diag, offdiag)`` with ``len(diag) == N`` and
    ``len(offdiag) == N - 1``, such that ``x p_n = e_{n-1} p_{n-1} + d_n p_n +
    e_n p_{n+1}``.
    """
    n = np.arange(N, dtype=float)
    Nn = 2 * n + alpha + beta
    with np.errstate(divide="ignore", invalid="ignore"):
        diag = (beta**2 - alpha**2) / (Nn * (Nn + 2))
    # 0/0 occurs only for n = 0 with alpha = beta = 0, where the limit is 0.
    diag = np.where(np.isfinite(diag), diag, 0.0)

    k = np.arange(N - 1, dtype=float)
    Nk = 2 * k + alpha + beta
    A = 2 * (k + 1) * (k + alpha + beta + 1) / ((Nk + 1) * (Nk + 2))
    C = 2 * (k + 1 + alpha) * (k + 1 + beta) / ((Nk + 2) * (Nk + 3))
    off = np.sqrt(A * C)
    return diag, off


def angular_matrix(m1: int, m2: int, c2: complex, N: int = 40) -> np.ndarray:
    """Complex-symmetric tridiagonal matrix whose eigenvalues are ``Ahat``."""
    alpha, beta = abs(m2), abs(m1)
    xdiag, xoff = jacobi_recurrence(alpha, beta, N)
    n = np.arange(N)
    l = 2 * n + alpha + beta

    H = np.zeros((N, N), dtype=complex)
    H[n, n] = l * (l + 2) - 0.5 * c2 * (1.0 - xdiag)
    idx = np.arange(N - 1)
    H[idx, idx + 1] = 0.5 * c2 * xoff
    H[idx + 1, idx] = 0.5 * c2 * xoff
    return H


def angular_spectrum(m1: int, m2: int, c2: complex, N: int = 40):
    """Full truncated spectrum, ordered by continuity from ``c2 = 0``.

    Returns ``(Ahat, vectors)``.  Ordering is obtained by tracking the
    eigenvectors' dominant Jacobi index rather than by sorting eigenvalues,
    so labels remain meaningful for complex ``c2``.
    """
    H = angular_matrix(m1, m2, c2, N)
    vals, vecs = np.linalg.eig(H)
    # label each eigenvector by its dominant Jacobi index n
    order = np.argsort(np.argmax(np.abs(vecs), axis=0), kind="stable")
    return vals[order], vecs[:, order]


def angular_eigenvalue(
    m1: int, m2: int, n: int, c2: complex, N: int = 40, n_steps: int = 0
) -> complex:
    """``Ahat`` for the branch that reduces to ``l(l+2)`` at ``c2 = 0``.

    ``n_steps > 0`` performs homotopy continuation in ``c2`` from 0, which is
    required once ``|c2|`` is large enough that dominant-index labelling is
    ambiguous (i.e. near angular-sector level repulsion).
    """
    if n_steps <= 0:
        vals, _ = angular_spectrum(m1, m2, c2, N)
        return vals[n]

    l = l_of_n(n, m1, m2)
    val = complex(l * (l + 2))
    for k in range(1, n_steps + 1):
        vals, _ = angular_spectrum(m1, m2, c2 * k / n_steps, N)
        val = vals[int(np.argmin(np.abs(vals - val)))]
    return val

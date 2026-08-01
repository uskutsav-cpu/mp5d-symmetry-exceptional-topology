"""Independent finite-difference angular solver (cross-check for the spectral one).

This deliberately shares *no* machinery with ``spheroidal5d``: no Jacobi
three-term recurrence, no orthonormal basis, no spectral expansion.  It
discretizes the regularized angular ODE directly on a uniform grid in
``u = cos^2 theta``.

Writing ``S = u^(alpha/2) (1-u)^(beta/2) f(u)`` with ``alpha = |m2|``,
``beta = |m1|``, the angular equation becomes the Jacobi-type ODE

    u(1-u) f'' + [ (alpha+1) - (alpha+beta+2) u ] f' + (c2 u / 4) f
        = -[ (Ahat - L(L+2)) / 4 ] f,     L = alpha + beta,

so ``Ahat = L(L+2) - 4 * eig(Op)``.  The ``u(1-u)`` factor vanishes at both
endpoints, so the first-order terms there act as the natural boundary
conditions and no explicit BC has to be imposed.
"""

from __future__ import annotations

import numpy as np

__all__ = ["angular_spectrum_fd"]


def angular_spectrum_fd(m1: int, m2: int, c2: complex, N: int = 400):
    """Return ``Ahat`` eigenvalues, ordered by ``Re`` at ``c2 = 0`` continuity.

    ``N`` is the number of grid points on ``u in [0, 1]`` inclusive.
    """
    alpha, beta = abs(m2), abs(m1)
    L = alpha + beta

    u = np.linspace(0.0, 1.0, N)
    h = u[1] - u[0]

    p = u * (1.0 - u)                      # coefficient of f''
    q = (alpha + 1) - (alpha + beta + 2) * u   # coefficient of f'
    rr = c2 * u / 4.0                       # zeroth order

    A = np.zeros((N, N), dtype=complex)
    # interior: standard second-order central differences
    i = np.arange(1, N - 1)
    A[i, i - 1] = p[i] / h**2 - q[i] / (2 * h)
    A[i, i] = -2 * p[i] / h**2 + rr[i]
    A[i, i + 1] = p[i] / h**2 + q[i] / (2 * h)
    # endpoints: p = 0 exactly, so only the first-order term survives.
    # Use second-order one-sided differences.
    A[0, 0] = q[0] * (-3.0) / (2 * h) + rr[0]
    A[0, 1] = q[0] * 4.0 / (2 * h)
    A[0, 2] = q[0] * (-1.0) / (2 * h)
    A[N - 1, N - 1] = q[-1] * 3.0 / (2 * h) + rr[-1]
    A[N - 1, N - 2] = q[-1] * (-4.0) / (2 * h)
    A[N - 1, N - 3] = q[-1] * 1.0 / (2 * h)

    vals = np.linalg.eigvals(A)
    Ahat = L * (L + 2) - 4.0 * vals
    # Discard the spurious high-frequency end of the FD spectrum and order the
    # resolved low-lying modes.
    Ahat = Ahat[np.argsort(np.abs(Ahat))]
    return Ahat

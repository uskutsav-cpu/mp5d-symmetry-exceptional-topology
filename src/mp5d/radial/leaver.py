"""Radial Solver B: Leaver / Hill-determinant method in the Leaver variable.

Why this and not collocation
----------------------------
Boundary-factored Chebyshev collocation fails for this problem, and the failure
is structural rather than a tuning issue (recorded in
``docs/FAILED_APPROACHES.md``).  After peeling ``exp(i Om r)``, the *ingoing*
solution behaves as ``exp(-2 i Om r)``, which for ``Im w < 0`` **decays** at
infinity instead of diverging.  Requiring the peeled function to be a polynomial
therefore fails to exclude it: both solutions are representable, and the
collocation matrix is numerically singular throughout the lower half plane.

Leaver's method selects the physical solution differently, by demanding that the
Frobenius coefficient sequence about the horizon be the **minimal** solution of
the recurrence.  That is the correct characterization of a QNM.

Construction
------------
With ``u = (r - r_+)/(r - r_-)`` the singular points map as

    r = r_+  -> u = 0        r = r_-   -> u = infinity
    r = inf  -> u = 1        r = -r_+  -> u = 2 r_+/(r_+ + r_-)
                             r = -r_-  -> u = (r_+ + r_-)/(2 r_-)

Both spurious points satisfy ``|u| > 1`` because ``r_+ > r_-``, so the series
about ``u = 0`` has radius of convergence exactly 1, with the irregular point
sitting on the boundary at ``u = 1``.  That is precisely Leaver's setting.

The singular behaviour factors as

    R = u^{-i sig} (1 - u)^{3/2} exp( i Om (r_+ - r_-) / (1 - u) ) * y(u)

using ``r - r_+ = u (r_+ - r_-)/(1 - u)`` and ``r - r_- = (r_+ - r_-)/(1 - u)``,
which are exact.  The remaining ``y = sum a_n u^n`` obeys a linear ODE with
*rational* coefficients in ``u``.

Recovering the recurrence
-------------------------
The polynomial coefficients ``A(u), B(u), C(u)`` are recovered **numerically and
exactly** rather than symbolically: the coefficient functions are multiplied by
a clearing factor ``Q(u)`` that removes every pole, evaluated on a circle in the
complex ``u`` plane, and inverse-FFT'd.  For a polynomial this recovers the
coefficients to machine precision, and the decay of the high-order tail is a
*self-check* that ``Q`` really cleared all poles (:func:`polynomial_coefficients`
raises if it did not).  Fully symbolic derivation with all parameters symbolic
was attempted first and did not terminate; see ``docs/FAILED_APPROACHES.md``.

The QNM condition is then that the banded (Hill) matrix built from the
recurrence be singular.
"""

from __future__ import annotations

import cmath

import numpy as np

from ..angular import angular_eigenvalue
from ..geometry import MPGeometry

__all__ = [
    "LeaverProblem",
    "polynomial_coefficients",
    "SOLVER_VERSION",
]

SOLVER_VERSION = "leaver-hill-B/1.0"


def polynomial_coefficients(
    fn, degree_bound: int, tol: float = 1e-8
) -> np.ndarray:
    """Exact polynomial coefficients of ``fn`` via FFT on the unit circle.

    ``fn`` must be a polynomial of degree <= ``degree_bound``.  Raises if the
    recovered tail does not vanish, which is the built-in check that the
    clearing factor removed every pole.

    Sampling uses the **unit** circle with a half-sample phase offset,
    ``u_k = exp(2 pi i (k + 1/2)/N)``.  Two reasons: a radius ``rho != 1``
    forces a division by ``rho^n`` during recovery, which for ``rho < 1`` and
    ``n ~ 100`` amplifies roundoff by ``10^50`` and destroys the result; and the
    half-sample offset guarantees no sample point ever lands on ``u = 1``, where
    the uncleared coefficient functions have their irregular singular point.
    """
    N = 1
    while N < 2 * (degree_bound + 1):
        N *= 2
    k = np.arange(N)
    pts = np.exp(2j * np.pi * (k + 0.5) / N)
    vals = fn(pts)
    if not np.all(np.isfinite(vals)):
        raise FloatingPointError("non-finite values on the evaluation circle")
    # f(u_k) = sum_j g_j e^{+2 pi i j k/N} with g_j = c_j e^{i pi j/N}.
    # Applying the forward transform: fft(f)_m = N g_m.  (ifft has the opposite
    # sign convention and would return the coefficients index-reversed.)
    coef = (np.fft.fft(vals) / N) * np.exp(-1j * np.pi * np.arange(N) / N)
    scale = np.max(np.abs(coef))
    if scale == 0:
        return np.zeros(degree_bound + 1, dtype=complex)
    tail = np.max(np.abs(coef[degree_bound + 1 :])) / scale
    if tail > tol:
        raise ValueError(
            f"coefficient tail {tail:.2e} exceeds {tol:.1e}: the clearing factor "
            "did not remove all poles, or degree_bound is too small"
        )
    return coef[: degree_bound + 1]


class LeaverProblem:
    """Hill-determinant formulation for one ``(m1, m2, l)`` sector."""

    def __init__(
        self,
        geo: MPGeometry,
        mu: float,
        m1: int,
        m2: int,
        ell: int,
        depth: int = 60,
        angular_N: int = 60,
    ):
        if geo.a == 0.0 and geo.b == 0.0:
            # r_- = 0 exactly; the construction still holds (u = 1 - r_+/r) but
            # r = 0 becomes a singular point of the radial equation.
            pass
        self.geo = geo
        self.mu = mu
        self.m1, self.m2, self.ell = m1, m2, ell
        self.depth = depth
        self.angular_N = angular_N

        rest = ell - abs(m1) - abs(m2)
        if rest < 0 or rest % 2 != 0:
            raise ValueError(f"l={ell} incompatible with (m1,m2)=({m1},{m2})")
        self.n_ang = rest // 2

    # -- coupled angular eigenvalue --------------------------------------
    def Lambda_of(self, omega: complex) -> complex:
        g = self.geo
        c2 = (omega**2 - self.mu**2) * (g.a**2 - g.b**2)
        Ahat = angular_eigenvalue(
            self.m1, self.m2, self.n_ang, c2, N=self.angular_N, n_steps=6
        )
        return Ahat - (omega**2 - self.mu**2) * g.b**2

    # -- ODE for y in the Leaver variable --------------------------------
    def _P(self, uu: np.ndarray, omega: complex, Lam: complex):
        g = self.geo
        a, b, M, mu = g.a, g.b, g.M, self.mu
        m1, m2 = self.m1, self.m2
        rp, rm = g.r_plus, g.r_minus
        d = rp - rm

        r = (rp - rm * uu) / (1.0 - uu)
        r2 = r * r
        q = a**2 + b**2 - M
        w = a**2 * b**2

        Delta = r2 + q + w / r2
        dDelta = 2.0 * r - 2.0 * w / (r2 * r)

        W = (r2 + a**2) * (r2 + b**2) * omega - m1 * a * (r2 + b**2) - m2 * b * (r2 + a**2)
        G = a * b * omega - a * m2 - b * m1
        V = (
            W**2 / (r2 * r2 * Delta)
            - G**2 / r2
            - (a**2 + b**2) * omega**2
            + 2.0 * omega * (a * m1 + b * m2)
            - mu**2 * r2
            - Lam
        )

        dudr = (1.0 - uu) ** 2 / d
        d2udr2 = -2.0 * (1.0 - uu) ** 3 / d**2

        A_u = Delta * dudr**2
        B_u = Delta * d2udr2 + (Delta / r + dDelta) * dudr
        C_u = V

        Om = cmath.sqrt(omega**2 - mu**2)
        if Om.real < 0:
            Om = -Om
        sig = (omega - m1 * g.Omega_a - m2 * g.Omega_b) / (2.0 * g.kappa)
        c = 1j * Om * d

        LF = -1j * sig / uu - 1.5 / (1.0 - uu) + c / (1.0 - uu) ** 2
        LFp = 1j * sig / uu**2 - 1.5 / (1.0 - uu) ** 2 + 2.0 * c / (1.0 - uu) ** 3
        FppF = LFp + LF * LF

        P2 = A_u
        P1 = 2.0 * A_u * LF + B_u
        P0 = A_u * FppF + B_u * LF + C_u
        return P2, P1, P0

    def _clearing(self, uu: np.ndarray) -> np.ndarray:
        """Q(u) removing every pole of P2, P1, P0."""
        g = self.geo
        rp, rm = g.r_plus, g.r_minus
        # poles: u = 0 (from Delta zero and from u^{-i sig}); u = 1 (exponential
        # and the compactification); r = 0 i.e. (rp - rm u) = 0; and the zeros of
        # Delta at r = -rp, -rm.
        Nr = rp - rm * uu
        r = Nr / (1.0 - uu)
        r2 = r * r
        q = g.a**2 + g.b**2 - g.M
        w = g.a**2 * g.b**2
        Delta_num = (r2 * r2 + q * r2 + w) * (1.0 - uu) ** 4  # polynomial in u
        return uu**2 * (1.0 - uu) ** 4 * Delta_num * Nr**4

    def poly_ABC(self, omega: complex, Lam: complex, degree_bound: int = 40):
        def mk(idx):
            def f(uu):
                P = self._P(uu, omega, Lam)
                return P[idx] * self._clearing(uu)
            return f

        A = polynomial_coefficients(mk(0), degree_bound)
        B = polynomial_coefficients(mk(1), degree_bound)
        C = polynomial_coefficients(mk(2), degree_bound)

        # The clearing factor Q carries a power of u (it must, to cancel the
        # u^{-i sig} prefactor and the zero of Delta at the horizon).  That power
        # is a common factor of A, B and C, and leaving it in would zero the
        # leading rows of the recurrence and make the Hill matrix singular for
        # every omega.  Divide it out; the ODE is unchanged.
        def valuation(arr: np.ndarray) -> int:
            big = np.abs(arr).max()
            if big == 0:
                return len(arr)
            nz = np.nonzero(np.abs(arr) > 1e-11 * big)[0]
            return int(nz[0]) if len(nz) else len(arr)

        k = min(valuation(A), valuation(B), valuation(C))
        if k > 0:
            A, B, C = A[k:], B[k:], C[k:]
        return A, B, C

    def hill_matrix(self, omega: complex, Lam: complex | None = None) -> np.ndarray:
        if Lam is None:
            Lam = self.Lambda_of(omega)
        A, B, C = self.poly_ABC(omega, Lam)
        K = self.depth

        n = np.arange(K)[:, None]
        m = np.arange(K)[None, :]
        Mat = np.zeros((K, K), dtype=complex)
        idxA = n - m + 2
        idxB = n - m + 1
        idxC = n - m
        Apad = np.concatenate([A, np.zeros(2 * K, dtype=complex)])
        Bpad = np.concatenate([B, np.zeros(2 * K, dtype=complex)])
        Cpad = np.concatenate([C, np.zeros(2 * K, dtype=complex)])
        okA = (idxA >= 0) & (idxA < len(Apad))
        okB = (idxB >= 0) & (idxB < len(Bpad))
        okC = (idxC >= 0) & (idxC < len(Cpad))
        Mat += np.where(okA, Apad[np.clip(idxA, 0, len(Apad) - 1)], 0) * (m * (m - 1))
        Mat += np.where(okB, Bpad[np.clip(idxB, 0, len(Bpad) - 1)], 0) * m
        Mat += np.where(okC, Cpad[np.clip(idxC, 0, len(Cpad) - 1)], 0)
        scale = np.max(np.abs(Mat), axis=1)
        scale[scale == 0] = 1.0
        return Mat / scale[:, None]

    def smallest_singular_value(self, omega: complex) -> float:
        S = np.linalg.svd(self.hill_matrix(omega), compute_uv=False)
        return float(S[-1] / S[0])

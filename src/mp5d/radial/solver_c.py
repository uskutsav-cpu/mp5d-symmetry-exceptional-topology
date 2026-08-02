"""Solver C: complex-scaled spectral solver.  Recurrence-free by construction.

Uses **none** of: the Frobenius coefficient recurrence, the Gaussian reduction,
the Leaver continued fraction, the Hill determinant, Wynn acceleration, or
Solver A/B's spectral residual.  Its root condition is the smallest singular
value of an independently assembled collocation matrix.

Why complex scaling, and why the earlier collocation failed
-----------------------------------------------------------
``docs/FAILED_APPROACHES.md`` #3 records that boundary-factored collocation on
the *real* axis fails: after peeling ``exp(i Om r)`` the ingoing solution goes as
``exp(-2 i Om r)``, which for ``Im w < 0`` **decays**, so demanding boundedness
does not exclude it and the matrix is singular everywhere.

Rotating the contour repairs exactly this.  On

    r(rho) = r_+ + rho e^{i theta},     rho in [0, L]

the ``rho``-dependent part of the outgoing exponent has real part

    Re[i Om e^{i theta}] rho = rho ( -Om_I cos(theta) - Om_R sin(theta) )

so the outgoing solution **decays** once

    tan(theta) > -Om_I / Om_R          (with Om_I < 0 for a damped mode)

while the ingoing solution, carrying the opposite sign, **grows**.  Requiring
decay at the far end therefore selects the physical solution.  It is the complex
scaling, not boundedness, that enforces the radiation condition.

Contour safety
--------------
``r = 0`` and ``r = -r_+`` lie at ``rho e^{i theta}`` real and negative, i.e.
``theta = pi``; ``+-r_-`` likewise sit on the real axis.  So for
``0 < theta < pi`` the contour never meets a singular point of the radial
equation, and ``Delta != 0`` all along it.  The branch of
``Om = sqrt(w^2 - mu^2)`` is fixed to ``Re Om >= 0``, matching Solvers A/B.

Discretization
--------------
``R = (r - r_+)^{-i sigma} f(rho)`` peels the ingoing horizon behaviour, leaving
``f`` analytic at ``rho = 0``.  ``f`` is expanded in Chebyshev polynomials on
``[0, L]``, collocated at interior Chebyshev-Gauss points, with a single explicit
row imposing ``f(L) = 0`` -- legitimate here because the physical solution is
exponentially small at the far end of the rotated contour, whereas the spurious
one is exponentially large.

The angular eigenvalue is recomputed at every ``omega`` and never frozen.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import chebyshev as Cheb

from ..angular import angular_eigenvalue
from ..geometry import MPGeometry

__all__ = ["SolverCProblem", "solve_qnm_c", "SolverCResult", "min_scaling_angle"]

SOLVER_C = "complex-scaled-collocation-C/1.0"


def min_scaling_angle(omega: complex, mu: float) -> float:
    """Smallest ``theta`` for which the outgoing solution decays on the contour."""
    Om = cmath.sqrt(omega**2 - mu**2)
    if Om.real < 0:
        Om = -Om
    return float(np.arctan2(-Om.imag, Om.real))


@dataclass
class SolverCResult:
    omega: complex
    Lambda: complex
    smallest_singular_value: float
    singular_value_gap: float
    residual: float
    theta: float
    L: float
    resolution: int
    converged: bool
    iterations: int
    solver: str = SOLVER_C


class SolverCProblem:
    """Assembles the complex-scaled collocation matrix for one sector."""

    def __init__(self, geo: MPGeometry, mu: float, m1: int, m2: int, ell: int,
                 theta: float = 1.0, L: float = 60.0, resolution: int = 200,
                 angular_N: int = 60):
        self.geo, self.mu = geo, mu
        self.m1, self.m2, self.ell = m1, m2, ell
        self.theta, self.L, self.N = theta, L, resolution
        self.angular_N = angular_N
        rest = ell - abs(m1) - abs(m2)
        if rest < 0 or rest % 2 != 0:
            raise ValueError(f"l={ell} incompatible with (m1,m2)=({m1},{m2})")
        self.k_ang = rest // 2

        # Chebyshev basis on [0, L]; interior Gauss points for the ODE rows
        j = np.arange(self.N - 1)
        t = np.cos(np.pi * (2 * j + 1) / (2 * (self.N - 1)))   # in (-1,1)
        self.rho = 0.5 * self.L * (t + 1.0)
        ident = np.eye(self.N)
        D1 = np.zeros((self.N, self.N))
        for k in range(self.N):
            dk = Cheb.chebder(ident[k])
            D1[: len(dk), k] = dk
        D2 = D1 @ D1
        # map d/dt -> d/drho
        self.sc = 2.0 / self.L
        V = Cheb.chebvander(t, self.N - 1)
        self.V, self.D1, self.D2 = V, V @ D1 * self.sc, V @ D2 * self.sc**2
        self.Vend = Cheb.chebvander(np.array([1.0]), self.N - 1)  # rho = L

    def Lambda_of(self, omega: complex) -> complex:
        g = self.geo
        c2 = (omega**2 - self.mu**2) * (g.a**2 - g.b**2)
        Ahat = angular_eigenvalue(self.m1, self.m2, self.k_ang, c2,
                                  N=self.angular_N, n_steps=6)
        return Ahat - (omega**2 - self.mu**2) * g.b**2

    def matrix(self, omega: complex, Lam: complex | None = None) -> np.ndarray:
        if Lam is None:
            Lam = self.Lambda_of(omega)
        g = self.geo
        a, b, M, mu = g.a, g.b, g.M, self.mu
        m1, m2 = self.m1, self.m2
        rp = g.r_plus
        e = np.exp(1j * self.theta)

        rho = self.rho
        r = rp + rho * e
        r2 = r * r
        Delta = r2 + (a**2 + b**2 - M) + (a**2 * b**2) / r2
        dDelta = 2.0 * r - 2.0 * (a**2 * b**2) / (r2 * r)
        W = ((r2 + a**2) * (r2 + b**2) * omega
             - m1 * a * (r2 + b**2) - m2 * b * (r2 + a**2))
        G = a * b * omega - a * m2 - b * m1
        Vpot = (W**2 / (r2 * r2 * Delta) - G**2 / r2
                - (a**2 + b**2) * omega**2 + 2.0 * omega * (a * m1 + b * m2)
                - mu**2 * r2 - Lam)

        sig = (omega - m1 * g.Omega_a - m2 * g.Omega_b) / (2.0 * g.kappa)
        p = -1j * sig  # R = (r-r_+)^{-i sigma} f = (rho e^{i th})^{p} f

        alpha = Delta
        beta = Delta / r + dDelta
        # ODE in rho:  alpha f_rr + [2 alpha p/rho + beta e] f_r
        #              + [alpha p(p-1)/rho^2 + beta e p/rho + V e^2] f = 0
        C2 = alpha
        C1 = 2.0 * alpha * p / rho + beta * e
        C0 = alpha * p * (p - 1.0) / rho**2 + beta * e * p / rho + Vpot * e * e

        Mrows = C2[:, None] * self.D2 + C1[:, None] * self.D1 + C0[:, None] * self.V
        scale = np.max(np.abs(Mrows), axis=1)
        scale[scale == 0] = 1.0
        Mrows = Mrows / scale[:, None]
        # far-end decay row: f(L) = 0.  Physical solution is exponentially small
        # there on the rotated contour; the spurious one is exponentially large.
        end = self.Vend / np.max(np.abs(self.Vend))
        return np.vstack([Mrows, end])

    def sigma_min(self, omega: complex) -> tuple[float, float]:
        S = np.linalg.svd(self.matrix(omega), compute_uv=False)
        return float(S[-1] / S[0]), float(S[-2] / S[-1])


def solve_qnm_c(a: float, b: float, mu: float, m1: int, m2: int, ell: int,
                initial_frequency: complex, theta: float = 1.0, L: float = 60.0,
                resolution: int = 200, M: float = 1.0, tol: float = 1e-11,
                maxiter: int = 60, angular_N: int = 60) -> SolverCResult:
    """Root of the complex-scaled collocation condition, by Muller on sigma_min."""
    geo = MPGeometry(a=a, b=b, M=M)
    prob = SolverCProblem(geo, mu, m1, m2, ell, theta, L, resolution, angular_N)

    def f(w: complex) -> complex:
        Mat = prob.matrix(w)
        U, S, Vh = np.linalg.svd(Mat)
        u = U[:, -1]
        v = Vh[-1, :].conj()
        return u.conj() @ (Mat @ v)   # analytic near a simple root

    omega = initial_frequency
    h = 1e-5 * max(abs(omega), 1.0)
    xs = [omega - h, omega + h, omega]
    fs = [f(x) for x in xs]
    ok = False
    it = 0
    while it < maxiter:
        it += 1
        x0, x1, x2 = xs[-3:]
        f0, f1, f2 = fs[-3:]
        q = (x2 - x1) / (x1 - x0)
        A = q * f2 - q * (1 + q) * f1 + q * q * f0
        B = (2 * q + 1) * f2 - (1 + q) ** 2 * f1 + q * q * f0
        Cq = (1 + q) * f2
        disc = cmath.sqrt(B * B - 4 * A * Cq)
        den = B + disc if abs(B + disc) > abs(B - disc) else B - disc
        if den == 0:
            break
        x3 = x2 - (x2 - x1) * 2 * Cq / den
        if x3.real < 0.02:
            x3 = complex(0.02, x3.imag)
        f3 = f(x3)
        xs.append(x3)
        fs.append(f3)
        if abs(x3 - x2) < tol * max(abs(x3), 1.0) or abs(f3) < tol:
            ok = True
            break
    omega = xs[-1]
    smin, gap = prob.sigma_min(omega)
    return SolverCResult(
        omega=omega, Lambda=prob.Lambda_of(omega), smallest_singular_value=smin,
        singular_value_gap=gap, residual=float(abs(fs[-1])), theta=theta, L=L,
        resolution=resolution, converged=ok, iterations=it,
    )

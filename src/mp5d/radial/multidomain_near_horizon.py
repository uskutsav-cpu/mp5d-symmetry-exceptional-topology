"""Multidomain near-horizon solver (Solver D): exterior complex scaling with a
resolved near-horizon domain.

Why this exists
---------------
Single-domain Solver C (``solver_c.py``) loses accuracy as ``r_+ -> r_-``.  The
diagnosis (``docs/NEAR_EXTREMAL_CANDIDATE.md``): the inner-horizon singularity
sits at distance ``Delta_r = r_+ - r_-`` from the domain endpoint ``rho = 0``,
while the domain has length ``L ~ 60``.  At ``r_2 = 0.44``, ``Delta_r = 0.094``,
so the singularity is ``1.6e-3`` of a domain length from the endpoint: the
Bernstein ellipse collapses and Chebyshev convergence with it.  Measured
consequence: Solvers A, B, C agreed only to ``3.5e-3`` there.

The fix is domain decomposition with a **scaled near-horizon coordinate**

    y = (r - r_+) / (r_+ - r_-)

on which the inner-horizon singularity sits at ``y = -1``, a distance of order
unity from the inner domain regardless of how small ``Delta_r`` becomes.

Construction
------------
* **Inner domain** ``y in [0, Y_m]``, kept on the real axis.  The horizon
  behaviour ``(r - r_+)^{-i sigma}`` is peeled analytically, so the remaining
  function is analytic at ``y = 0``.
* **Outer domain** ``t in [0, L]`` with ``r = r_m + t e^{i theta}``, complex
  scaled so the outgoing solution decays (same criterion as Solver C:
  ``tan theta > -Om_I / Om_R``).  Scaling begins only *after* the near-horizon
  real region, which is exterior complex scaling proper.
* **Matching** at ``r_m = r_+ + Y_m Delta_r``: continuity of the function and of
  its radial derivative, both expressed in ``r`` so the two parametrizations are
  compared on equal footing.

This is recurrence-free: no Frobenius recurrence, no Gaussian reduction, no
continued fraction, no Hill determinant, no Wynn acceleration.  ``Lambda`` is
recomputed from the angular solver at every ``omega``.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import chebyshev as Cheb

from ..angular import angular_eigenvalue
from ..geometry import MPGeometry

__all__ = ["MultidomainProblem", "solve_qnm_multidomain", "MultidomainResult"]

SOLVER_D = "multidomain-near-horizon-D/1.0"


def _cheb_blocks(N: int, length: float):
    """Value / d/ds / d2/ds2 at interior Gauss points on ``s in [0, length]``."""
    j = np.arange(N - 1)
    t = np.cos(np.pi * (2 * j + 1) / (2 * (N - 1)))
    s = 0.5 * length * (t + 1.0)
    ident = np.eye(N)
    D1c = np.zeros((N, N))
    for k in range(N):
        dk = Cheb.chebder(ident[k])
        D1c[: len(dk), k] = dk
    D2c = D1c @ D1c
    sc = 2.0 / length
    V = Cheb.chebvander(t, N - 1)
    return s, V, (V @ D1c) * sc, (V @ D2c) * sc**2, D1c, sc


def _edge(N: int, at: float, D1c, sc):
    """Row vectors for value and d/ds at a domain end (``at = -1`` or ``+1``)."""
    v = Cheb.chebvander(np.array([at]), N - 1)
    return v, (v @ D1c) * sc


@dataclass
class MultidomainResult:
    omega: complex
    Lambda: complex
    smallest_singular_value: float
    residual: float
    theta: float
    Y_m: float
    L: float
    n_inner: int
    n_outer: int
    converged: bool
    iterations: int
    solver: str = SOLVER_D


class MultidomainProblem:
    """Two-domain complex-scaled collocation with a resolved near-horizon region."""

    def __init__(self, geo: MPGeometry, mu: float, m1: int, m2: int, ell: int,
                 theta: float = 1.0, Y_m: float = 8.0, L: float = 60.0,
                 n_inner: int = 90, n_outer: int = 200, angular_N: int = 60):
        self.geo, self.mu = geo, mu
        self.m1, self.m2, self.ell = m1, m2, ell
        self.theta, self.Y_m, self.L = theta, Y_m, L
        self.n_in, self.n_out = n_inner, n_outer
        self.angular_N = angular_N
        rest = ell - abs(m1) - abs(m2)
        if rest < 0 or rest % 2 != 0:
            raise ValueError(f"l={ell} incompatible with (m1,m2)=({m1},{m2})")
        self.k_ang = rest // 2

        self.dr = geo.r_plus - geo.r_minus          # the resolved scale
        self.r_m = geo.r_plus + Y_m * self.dr       # matching radius

        # inner: y in [0, Y_m];  r = r_+ + y * dr
        (self.y, self.Vi, self.D1i, self.D2i, self.D1ci, self.sci) = _cheb_blocks(
            n_inner, Y_m)
        # outer: t in [0, L];    r = r_m + t e^{i theta}
        (self.t, self.Vo, self.D1o, self.D2o, self.D1co, self.sco) = _cheb_blocks(
            n_outer, L)

        self.Vi_end, self.D1i_end = _edge(n_inner, 1.0, self.D1ci, self.sci)   # y = Y_m
        self.Vo_0, self.D1o_0 = _edge(n_outer, -1.0, self.D1co, self.sco)      # t = 0
        self.Vo_end, _ = _edge(n_outer, 1.0, self.D1co, self.sco)              # t = L

    def Lambda_of(self, omega: complex) -> complex:
        g = self.geo
        c2 = (omega**2 - self.mu**2) * (g.a**2 - g.b**2)
        Ahat = angular_eigenvalue(self.m1, self.m2, self.k_ang, c2,
                                  N=self.angular_N, n_steps=6)
        return Ahat - (omega**2 - self.mu**2) * g.b**2

    def _potential(self, r, omega, Lam):
        g = self.geo
        a, b, M, mu = g.a, g.b, g.M, self.mu
        m1, m2 = self.m1, self.m2
        r2 = r * r
        Delta = r2 + (a**2 + b**2 - M) + (a**2 * b**2) / r2
        dDelta = 2.0 * r - 2.0 * (a**2 * b**2) / (r2 * r)
        W = ((r2 + a**2) * (r2 + b**2) * omega
             - m1 * a * (r2 + b**2) - m2 * b * (r2 + a**2))
        G = a * b * omega - a * m2 - b * m1
        V = (W**2 / (r2 * r2 * Delta) - G**2 / r2
             - (a**2 + b**2) * omega**2 + 2.0 * omega * (a * m1 + b * m2)
             - mu**2 * r2 - Lam)
        return Delta, Delta / r + dDelta, V

    def matrix(self, omega: complex, Lam: complex | None = None) -> np.ndarray:
        if Lam is None:
            Lam = self.Lambda_of(omega)
        g = self.geo
        rp = g.r_plus
        sig = (omega - self.m1 * g.Omega_a - self.m2 * g.Omega_b) / (2.0 * g.kappa)
        p = -1j * sig
        e = np.exp(1j * self.theta)
        Ni, No = self.n_in, self.n_out

        # ---- inner domain: r = r_+ + y*dr, R = (y*dr)^p f(y)
        y = self.y
        r_in = rp + y * self.dr
        al, be, V = self._potential(r_in, omega, Lam)
        # d/dr = (1/dr) d/dy
        C2 = al / self.dr**2
        C1 = (2.0 * al * p / y) / self.dr**2 + be / self.dr
        C0 = al * p * (p - 1.0) / (y * self.dr) ** 2 + be * p / (y * self.dr) + V
        Min = C2[:, None] * self.D2i + C1[:, None] * self.D1i + C0[:, None] * self.Vi
        Min /= np.maximum(np.max(np.abs(Min), axis=1), 1e-300)[:, None]

        # ---- outer domain: r = r_m + t e^{i th}, R = g(t) with no extra peel
        t = self.t
        r_out = self.r_m + t * e
        al_o, be_o, V_o = self._potential(r_out, omega, Lam)
        C2o = al_o / e**2
        C1o = be_o / e
        Mout = C2o[:, None] * self.D2o + C1o[:, None] * self.D1o + V_o[:, None] * self.Vo
        Mout /= np.maximum(np.max(np.abs(Mout), axis=1), 1e-300)[:, None]

        # ---- assemble
        rows = []
        rows.append(np.hstack([Min[: Ni - 1], np.zeros((Ni - 1, No), dtype=complex)]))
        rows.append(np.hstack([np.zeros((No - 2, Ni), dtype=complex), Mout[: No - 2]]))

        # matching at r_m.  Inner R = (Y_m*dr)^p f, dR/dr = (Y_m dr)^p (f' /dr + p f/(Y_m dr))
        pref = (self.Y_m * self.dr) ** p
        val_in = pref * self.Vi_end
        der_in = pref * (self.D1i_end / self.dr + (p / (self.Y_m * self.dr)) * self.Vi_end)
        val_out = self.Vo_0
        der_out = self.D1o_0 / e
        mrow1 = np.hstack([val_in, -val_out])
        mrow2 = np.hstack([der_in, -der_out])
        for mr in (mrow1, mrow2):
            rows.append(mr / max(np.max(np.abs(mr)), 1e-300))

        decay = np.hstack([np.zeros((1, Ni), dtype=complex), self.Vo_end])
        rows.append(decay / max(np.max(np.abs(decay)), 1e-300))
        return np.vstack(rows)

    def sigma_min(self, omega: complex) -> float:
        S = np.linalg.svd(self.matrix(omega), compute_uv=False)
        return float(S[-1] / S[0])


def solve_qnm_multidomain(a: float, b: float, mu: float, m1: int, m2: int, ell: int,
                          initial_frequency: complex, theta: float = 1.0,
                          Y_m: float = 8.0, L: float = 60.0, n_inner: int = 90,
                          n_outer: int = 200, M: float = 1.0, tol: float = 1e-11,
                          maxiter: int = 50, angular_N: int = 60) -> MultidomainResult:
    geo = MPGeometry(a=a, b=b, M=M)
    prob = MultidomainProblem(geo, mu, m1, m2, ell, theta, Y_m, L,
                              n_inner, n_outer, angular_N)

    def f(w: complex) -> complex:
        Mat = prob.matrix(w)
        U, S, Vh = np.linalg.svd(Mat)
        return U[:, -1].conj() @ (Mat @ Vh[-1, :].conj())

    omega = initial_frequency
    h = 1e-5 * max(abs(omega), 1.0)
    xs = [omega - h, omega + h, omega]
    fs = [f(x) for x in xs]
    ok, it = False, 0
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
    return MultidomainResult(
        omega=omega, Lambda=prob.Lambda_of(omega),
        smallest_singular_value=prob.sigma_min(omega), residual=float(abs(fs[-1])),
        theta=theta, Y_m=Y_m, L=L, n_inner=n_inner, n_outer=n_outer,
        converged=ok, iterations=it,
    )

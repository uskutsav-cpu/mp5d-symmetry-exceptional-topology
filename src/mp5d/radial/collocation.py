"""Radial Solver A: boundary-factored Chebyshev collocation.

.. warning::
   **NOT FUNCTIONAL as a QNM eigenvalue solver.**  The selection principle is
   wrong for ``Im w < 0``; see the "Why this fails" section below and
   ``docs/FAILED_APPROACHES.md`` entry 3.  ``solve_qnm`` will not converge to a
   quasinormal frequency.  The coefficient assembly is correct and is reused by
   ``mp5d.radial.leaver``.

Formulation
-----------
The radial equation of ``docs/CONVENTIONS.md`` §4 is

    Delta R'' + (Delta/r + Delta') R' + V(r) R = 0

with (this closed form is used because it is manifestly finite at large ``r``)

    Delta = r^2 + (a^2 + b^2 - M) + a^2 b^2 / r^2
    V     = W^2/(r^4 Delta) - G^2/r^2 - (a^2+b^2) w^2 + 2w(a m1 + b m2)
            - mu^2 r^2 - Lambda

The two physical behaviours are factored out analytically:

    R(r) = exp(i Om r) (r - r_+)^{-i sig} r^{-3/2 + i sig} f(x)

* ``(r - r_+)^{-i sig}`` with ``sig = (w - m1 Om_a - m2 Om_b)/(2 kappa)`` is the
  **ingoing** horizon behaviour (derived in ``structure.py``, not assumed);
* ``exp(i Om r) r^{-3/2}`` with ``Om = sqrt(w^2 - mu^2)`` is the **outgoing**
  behaviour at infinity.  The power is exactly ``-3/2`` with no Coulomb phase,
  because the radial potential is exactly even in ``r`` (claim C18).  The extra
  ``r^{+i sig}`` compensates the large-``r`` tail of the horizon factor so that
  ``f -> const`` at infinity.

The compactified coordinate is

    x = 1 - 2 r_+ / r,      x = -1 at the horizon,  x = +1 at infinity.

Discretization
--------------
``f`` is expanded in Chebyshev polynomials and collocated at the **interior**
Chebyshev-Gauss points.  Interior collocation is deliberate: both endpoints are
singular points of the equation for ``f``, and evaluating there would require
the exact indicial cancellation to happen in floating point.  Requiring ``f`` to
be a polynomial already enforces analyticity at both ends, which *is* the
boundary condition once the singular factors are peeled off.  Nothing is imposed
by hand.

The result is a nonlinear eigenvalue problem ``M(w) c = 0``.  ``Lambda`` is
**not** frozen: it is recomputed from the angular solver at every iterate, with
spheroidicity ``c2 = (w^2 - mu^2)(a^2 - b^2)``, so the coupled dependence
``Lambda = Lambda(w, a, b, mu, m1, m2, l)`` is respected throughout.

Why this fails
--------------
After peeling ``exp(i Om r)``, the *ingoing* solution behaves as
``exp(-2 i Om r)``, which for ``Im w < 0`` **decays** at infinity rather than
diverging.  Demanding that the peeled function be a bounded polynomial therefore
fails to exclude it: both solutions are representable and the matrix is
numerically singular throughout the lower half plane, with no isolated root.
This is structural and cannot be tuned away.

Root finding is Newton on the smallest singular triplet: with ``M v = s u``,
the analytic function ``g(w) = u^H M(w) v`` vanishes at a QNM and
``w <- w - g / (u^H M'(w) v)`` converges quadratically at a simple root.
"""

from __future__ import annotations

import cmath
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
from numpy.polynomial import chebyshev as C

from ..angular import angular_eigenvalue
from ..geometry import MPGeometry

__all__ = ["QNMResult", "RadialCollocation", "solve_qnm", "SOLVER_VERSION"]

SOLVER_VERSION = "collocation-A/1.0"


@dataclass
class QNMResult:
    """Structured QNM result.  Serializable via :meth:`to_dict`."""

    omega: complex
    Lambda: complex
    residual: float
    boundary_residual_horizon: float
    boundary_residual_infinity: float
    convergence: list[float]
    condition_estimate: float
    smallest_singular_value: float
    singular_value_gap: float
    precision: str
    resolution: int
    mode_labels: dict[str, int]
    parameters: dict[str, float]
    solver: str
    commit: str
    converged: bool
    eigenvector: np.ndarray | None = field(default=None, repr=False)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.pop("eigenvector", None)
        d["omega"] = [self.omega.real, self.omega.imag]
        d["Lambda"] = [self.Lambda.real, self.Lambda.imag]
        return d


def _cheb_operators(N: int, x: np.ndarray):
    """Value / first / second derivative matrices in the Chebyshev coefficient basis."""
    ident = np.eye(N)
    D1 = np.zeros((N, N))
    for k in range(N):
        dk = C.chebder(ident[k])
        D1[: len(dk), k] = dk
    D2 = D1 @ D1
    V = C.chebvander(x, N - 1)
    return V, V @ D1, V @ D2


class RadialCollocation:
    """Assembles ``M(w)`` for one ``(m1, m2, l)`` sector."""

    def __init__(
        self,
        geo: MPGeometry,
        mu: float,
        m1: int,
        m2: int,
        ell: int,
        resolution: int = 80,
        angular_N: int = 60,
    ):
        self.geo = geo
        self.mu = mu
        self.m1, self.m2, self.ell = m1, m2, ell
        self.N = resolution
        self.angular_N = angular_N

        # angular radial-node index n from l = 2n + |m1| + |m2|
        rest = ell - abs(m1) - abs(m2)
        if rest < 0 or rest % 2 != 0:
            raise ValueError(
                f"l={ell} incompatible with (m1,m2)=({m1},{m2}): "
                "need l - |m1| - |m2| >= 0 and even"
            )
        self.n_ang = rest // 2

        j = np.arange(self.N)
        self.x = np.cos(np.pi * (2 * j + 1) / (2 * self.N))  # interior Gauss points
        self.V, self.D1, self.D2 = _cheb_operators(self.N, self.x)

    # -- physics ----------------------------------------------------------
    def Lambda_of(self, omega: complex) -> complex:
        """Coupled angular separation constant at this ``omega``."""
        g = self.geo
        c2 = (omega**2 - self.mu**2) * (g.a**2 - g.b**2)
        Ahat = angular_eigenvalue(
            self.m1, self.m2, self.n_ang, c2, N=self.angular_N, n_steps=6
        )
        return Ahat - (omega**2 - self.mu**2) * g.b**2

    def _coeffs(self, omega: complex, Lam: complex, x: np.ndarray):
        g = self.geo
        a, b, M, mu = g.a, g.b, g.M, self.mu
        m1, m2 = self.m1, self.m2
        rp = g.r_plus

        r = 2.0 * rp / (1.0 - x)
        r2 = r * r

        Delta = r2 + (a**2 + b**2 - M) + (a**2 * b**2) / r2
        dDelta = 2.0 * r - 2.0 * (a**2 * b**2) / (r2 * r)

        W = (r2 + a**2) * (r2 + b**2) * omega - m1 * a * (r2 + b**2) - m2 * b * (r2 + a**2)
        G = a * b * omega - a * m2 - b * m1

        Vpot = (
            W**2 / (r2 * r2 * Delta)
            - G**2 / r2
            - (a**2 + b**2) * omega**2
            + 2.0 * omega * (a * m1 + b * m2)
            - mu**2 * r2
            - Lam
        )

        Om = cmath.sqrt(omega**2 - mu**2)
        if Om.real < 0:
            Om = -Om
        sig = (omega - m1 * g.Omega_a - m2 * g.Omega_b) / (2.0 * g.kappa)

        # A'/A and A''/A for A = exp(i Om r) (r-rp)^{-i sig} r^{-3/2 + i sig}
        p = -1.5 + 1j * sig
        LA = 1j * Om - 1j * sig / (r - rp) + p / r
        LAp = 1j * sig / (r - rp) ** 2 - p / r2

        alpha = Delta
        beta = Delta / r + dDelta

        P2 = alpha
        P1 = 2.0 * alpha * LA + beta
        P0 = alpha * (LAp + LA * LA) + beta * LA + Vpot

        xr = (1.0 - x) ** 2 / (2.0 * rp)
        xrr = -((1.0 - x) ** 3) / (2.0 * rp**2)

        C2c = P2 * xr**2
        C1c = P2 * xrr + P1 * xr
        C0c = P0
        return C2c, C1c, C0c

    def matrix(self, omega: complex, Lam: complex | None = None) -> np.ndarray:
        if Lam is None:
            Lam = self.Lambda_of(omega)
        c2, c1, c0 = self._coeffs(omega, Lam, self.x)
        M = c2[:, None] * self.D2 + c1[:, None] * self.D1 + c0[:, None] * self.V
        # row equilibration: does not change the null space, improves conditioning
        scale = np.max(np.abs(M), axis=1)
        scale[scale == 0] = 1.0
        return M / scale[:, None]

    # -- residuals --------------------------------------------------------
    def residual_on_finer_grid(self, omega: complex, coeff: np.ndarray) -> float:
        """Evaluate the ODE residual at points that were *not* collocated."""
        Lam = self.Lambda_of(omega)
        k = np.arange(2 * self.N)
        xf = np.cos(np.pi * (2 * k + 1) / (4 * self.N))
        Vf, D1f, D2f = _cheb_operators(self.N, xf)
        c2, c1, c0 = self._coeffs(omega, Lam, xf)
        res = c2 * (D2f @ coeff) + c1 * (D1f @ coeff) + c0 * (Vf @ coeff)
        nrm = np.max(np.abs(Vf @ coeff))
        scale = np.maximum(np.abs(c2) + np.abs(c1) + np.abs(c0), 1e-300)
        return float(np.max(np.abs(res) / scale) / max(nrm, 1e-300))

    def boundary_residuals(self, coeff: np.ndarray) -> tuple[float, float]:
        """How well the peeled solution matches pure ingoing/outgoing behaviour.

        With the factors removed, ``f`` must be analytic at both ends; the
        diagnostic is the decay of the Chebyshev tail, which measures exactly
        that.  Reported separately for the two endpoints by evaluating the
        expansion's convergence toward ``x = -1`` and ``x = +1``.
        """
        tail = np.abs(coeff[-max(3, self.N // 10):]).max() / max(
            np.abs(coeff).max(), 1e-300
        )
        fm = abs(C.chebval(-1.0, coeff))
        fp = abs(C.chebval(1.0, coeff))
        nrm = max(np.abs(C.chebval(self.x, coeff)).max(), 1e-300)
        # a genuinely singular solution would blow up at an endpoint
        return float(tail * fm / nrm), float(tail * fp / nrm)


def solve_qnm(
    a: float,
    b: float,
    field_mass: float,
    m1: int,
    m2: int,
    ell: int,
    overtone: int = 0,
    initial_frequency: complex | None = None,
    precision: int = 80,
    resolution: int = 80,
    M: float = 1.0,
    angular_N: int = 60,
    tol: float = 1e-12,
    maxiter: int = 60,
    commit: str = "unknown",
) -> QNMResult:
    """Solve for one quasinormal frequency.

    ``precision`` is accepted for interface stability; Solver A runs in double
    precision and reports ``"double"``.  Precision escalation is provided by
    Solver B (``mp5d.radial.leaver``), which is the arbitrary-precision path.
    """
    geo = MPGeometry(a=a, b=b, M=M)
    if not geo.has_horizon:
        raise ValueError(f"no horizon for a={a}, b={b}, M={M}")

    prob = RadialCollocation(
        geo, field_mass, m1, m2, ell, resolution=resolution, angular_N=angular_N
    )

    omega = initial_frequency if initial_frequency is not None else _seed(geo, ell, overtone)

    history: list[float] = []
    converged = False
    u = v = None
    svals = np.array([1.0, 1.0])

    for _ in range(maxiter):
        Mat = prob.matrix(omega)
        U, S, Vh = np.linalg.svd(Mat)
        svals = S
        u = U[:, -1]
        v = Vh[-1, :].conj()
        g = u.conj() @ (Mat @ v)
        history.append(float(abs(g)))
        if abs(g) < tol:
            converged = True
            break
        h = 1e-7 * max(abs(omega), 1.0)
        dM = (prob.matrix(omega + h) - prob.matrix(omega - h)) / (2.0 * h)
        den = u.conj() @ (dM @ v)
        if den == 0 or not np.isfinite(abs(den)):
            break
        step = -g / den
        # damp wild steps
        if abs(step) > 0.5 * max(abs(omega), 1.0):
            step *= 0.5 * max(abs(omega), 1.0) / abs(step)
        omega = omega + step
        if abs(step) < 1e-15 * max(abs(omega), 1.0):
            converged = True
            break

    coeff = v if v is not None else np.zeros(prob.N)
    Lam = prob.Lambda_of(omega)
    res = prob.residual_on_finer_grid(omega, coeff)
    bh, bi = prob.boundary_residuals(coeff)
    gap = float(svals[-2] / svals[-1]) if svals[-1] > 0 else float("inf")

    return QNMResult(
        omega=complex(omega),
        Lambda=complex(Lam),
        residual=res,
        boundary_residual_horizon=bh,
        boundary_residual_infinity=bi,
        convergence=history,
        condition_estimate=float(svals[0] / max(svals[-1], 1e-300)),
        smallest_singular_value=float(svals[-1]),
        singular_value_gap=gap,
        precision="double",
        resolution=resolution,
        mode_labels={"m1": m1, "m2": m2, "l": ell, "n": prob.n_ang, "overtone": overtone},
        parameters={
            "a": a, "b": b, "s": geo.s, "delta": geo.delta,
            "mu": field_mass, "M": M,
        },
        solver=SOLVER_VERSION,
        commit=commit,
        converged=converged,
        eigenvector=coeff,
    )


def _seed(geo: MPGeometry, ell: int, overtone: int) -> complex:
    """Crude eikonal-flavoured seed; refined by continuation in practice."""
    rp = geo.r_plus
    return complex((ell + 1.0) / (1.6 * rp), -(overtone + 0.5) * 0.9 / rp)

"""Augmented exceptional-point conditions for the MP5D radial spectrum.

An EP2 is a *repeated root* of the spectral condition, not a pair of nearby
frequencies.  With ``F(omega, p)`` the spectral condition (Solver A's continued
fraction), the conditions are

    F = 0,        dF/domega = 0

Both are complex, so the augmented system is **4 real equations**.  The unknowns
are ``omega`` (2 real) plus **2 real parameters** drawn from ``(s, delta, mu)``.
Hence an EP2 is codimension 2: in the three-parameter physical space the
exceptional set is generically a **one-dimensional locus**, and one parameter
must be held fixed to isolate points.

Total derivative
----------------
The angular eigenvalue depends on ``omega``:

    dF/domega = partial_omega F + partial_Lambda F * dLambda/domega

This is handled automatically and exactly, because ``cf_value`` recomputes
``Lambda(omega)`` from the angular solver on every evaluation.  A finite
difference of ``cf_value`` in ``omega`` therefore already carries the
``dLambda/domega`` term -- ``Lambda`` is never frozen.  Freezing it is the
classic way to produce spurious "EPs", so this is asserted in the tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..geometry import MPGeometry
from ..radial.leaver import LeaverProblem
from ..radial.recurrence import continued_fraction_inverted, reduce_to_three_term

__all__ = ["EPCandidate", "spectral_condition", "d_spectral_condition",
           "solve_ep2", "PARAM_NAMES"]

PARAM_NAMES = ("s", "delta", "mu")


def _problem(s: float, delta: float, mu: float, m1: int, m2: int, ell: int,
             M: float = 1.0, angular_N: int = 50) -> LeaverProblem:
    geo = MPGeometry(a=s + delta, b=s - delta, M=M)
    return LeaverProblem(geo, mu, m1, m2, ell, depth=400, angular_N=angular_N)


def spectral_condition(omega: complex, s: float, delta: float, mu: float,
                       m1: int, m2: int, ell: int, overtone: int = 0,
                       depth: int = 300, degree_bound: int = 32) -> complex:
    """``F(omega, p)``: Solver A's continued-fraction condition.

    ``Lambda`` is recomputed inside, so the omega-dependence is total.
    """
    prob = _problem(s, delta, mu, m1, m2, ell)
    A, B, C = prob.poly_ABC(omega, prob.Lambda_of(omega), degree_bound=degree_bound)
    al, be, ga = reduce_to_three_term(A, B, C, depth=depth)
    return continued_fraction_inverted(al, be, ga, overtone, depth)


def d_spectral_condition(omega: complex, s: float, delta: float, mu: float,
                         m1: int, m2: int, ell: int, overtone: int = 0,
                         depth: int = 300, h: float = 1e-6) -> complex:
    """``dF/domega`` by a central difference.  ``F`` is analytic in ``omega``."""
    hh = h * max(abs(omega), 1.0)
    fp = spectral_condition(omega + hh, s, delta, mu, m1, m2, ell, overtone, depth)
    fm = spectral_condition(omega - hh, s, delta, mu, m1, m2, ell, overtone, depth)
    return (fp - fm) / (2.0 * hh)


@dataclass
class EPCandidate:
    omega: complex
    params: dict[str, float]
    fixed_param: str
    residual_F: float
    residual_dF: float
    converged: bool
    iterations: int
    m1: int
    m2: int
    ell: int
    overtone: int
    history: list[float] = field(default_factory=list)

    def to_dict(self):
        return {
            "omega": [self.omega.real, self.omega.imag],
            "params": self.params, "fixed_param": self.fixed_param,
            "residual_F": self.residual_F, "residual_dF": self.residual_dF,
            "converged": self.converged, "iterations": self.iterations,
            "mode_labels": {"m1": self.m1, "m2": self.m2, "l": self.ell,
                            "overtone": self.overtone},
            "history": self.history,
        }


def solve_ep2(omega0: complex, s0: float, delta0: float, mu0: float,
              m1: int, m2: int, ell: int, overtone: int = 0,
              free: tuple[str, str] = ("delta", "mu"),
              depth: int = 300, tol: float = 1e-10, maxiter: int = 40,
              step_limit: float = 0.25) -> EPCandidate:
    """Newton on ``[F, dF/domega] = 0`` in ``(omega, p_free1, p_free2)``.

    ``free`` names the two real parameters allowed to vary; the third is held
    fixed, which is what makes the EP2 an isolated point rather than a curve.
    """
    p = {"s": s0, "delta": delta0, "mu": mu0}
    fixed = [n for n in PARAM_NAMES if n not in free]
    if len(fixed) != 1:
        raise ValueError(f"exactly one parameter must be fixed, got free={free}")
    fixed_name = fixed[0]

    x = np.array([omega0.real, omega0.imag, p[free[0]], p[free[1]]], dtype=float)
    hist: list[float] = []
    converged = False
    it = 0

    def residual(vec):
        om = complex(vec[0], vec[1])
        q = dict(p)
        q[free[0]], q[free[1]] = vec[2], vec[3]
        if q["mu"] < 0:
            q["mu"] = 0.0
        geo = MPGeometry(a=q["s"] + q["delta"], b=q["s"] - q["delta"], M=1.0)
        if not geo.has_horizon or geo.extremality <= 1e-6:
            raise ValueError("left the sub-extremal region")
        F = spectral_condition(om, q["s"], q["delta"], q["mu"], m1, m2, ell,
                               overtone, depth)
        dF = d_spectral_condition(om, q["s"], q["delta"], q["mu"], m1, m2, ell,
                                  overtone, depth)
        # normalize dF by the local scale of F so the two blocks are comparable
        return np.array([F.real, F.imag, dF.real, dF.imag]), F, dF

    while it < maxiter:
        it += 1
        try:
            R, F, dF = residual(x)
        except Exception:
            break
        nrm = float(np.linalg.norm(R))
        hist.append(nrm)
        if nrm < tol:
            converged = True
            break
        # numerical Jacobian
        J = np.zeros((4, 4))
        for j in range(4):
            dx = 1e-6 * max(abs(x[j]), 1e-3)
            xp = x.copy()
            xp[j] += dx
            try:
                Rp, _, _ = residual(xp)
            except Exception:
                J[:, j] = 0.0
                continue
            J[:, j] = (Rp - R) / dx
        try:
            step = np.linalg.solve(J, -R)
        except np.linalg.LinAlgError:
            break
        n = float(np.linalg.norm(step))
        if n > step_limit:
            step *= step_limit / n
        x = x + step

    om = complex(x[0], x[1])
    q = dict(p)
    q[free[0]], q[free[1]] = x[2], x[3]
    try:
        F = spectral_condition(om, q["s"], q["delta"], q["mu"], m1, m2, ell,
                               overtone, depth)
        dF = d_spectral_condition(om, q["s"], q["delta"], q["mu"], m1, m2, ell,
                                  overtone, depth)
    except Exception:
        F = dF = complex("nan")
    return EPCandidate(
        omega=om, params=q, fixed_param=fixed_name,
        residual_F=float(abs(F)), residual_dF=float(abs(dF)),
        converged=converged, iterations=it, m1=m1, m2=m2, ell=ell,
        overtone=overtone, history=hist,
    )

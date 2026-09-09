"""Fail-closed EP2 evidence on finite numerical problems, not continuum proofs."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

import numpy as np

WITHDRAWN = frozenset({"min_root_separation", "min_abs_discriminant", "root_separation",
                       "a1_over_a2", "a1/a2", "discriminant_independent"})


class Verdict(StrEnum):
    INCONCLUSIVE = "INCONCLUSIVE"
    DISTINCT_AT_TESTED_POINT = "DISTINCT_AT_TESTED_POINT"
    EP2_NUMERICALLY_SUPPORTED = "EP2_NUMERICALLY_SUPPORTED"


def _matrix(a: Any, name: str) -> np.ndarray:
    a = np.asarray(a, dtype=complex)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or a.shape[0] < 2:
        raise ValueError(f"{name} must be square of dimension >= 2")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"non-finite {name}")
    return a


def jordan_test(T: Any, dT: Any, ddT: Any | None = None, *,
                relative_tol: float = 1e-8, absolute_tol: float = 1e-12) -> dict:
    """Test nullity one, a first generalized vector, and obstruction at length 3.

    For an analytic matrix pencil T(w), chains satisfy
    T u0=0, T u1=-T' u0, T u2=-T' u1-(T''/2)u0.
    A nonzero left-nullvector projection of the last RHS obstructs length 3.
    Setting ddT=None does NOT assert a linear pencil; exact order stays unknown.
    All comparisons are numerical, and no discretization error is enclosed.
    """
    T, dT = _matrix(T, "T"), _matrix(dT, "dT")
    if T.shape != dT.shape:
        raise ValueError("matrix shapes differ")
    if relative_tol <= 0 or absolute_tol <= 0:
        raise ValueError("tolerances must be positive")
    if ddT is not None:
        ddT = _matrix(ddT, "ddT")
        if ddT.shape != T.shape:
            raise ValueError("second derivative shape differs")
    U, s, Vh = np.linalg.svd(T)
    scale = max(float(s[0]), np.finfo(float).tiny)
    threshold = max(absolute_tol, relative_tol * scale)
    nullity = int(np.count_nonzero(s <= threshold))
    u0, v0 = Vh[-1].conj(), U[:, -1]
    rhs = -dT @ u0
    rhsnorm = float(np.linalg.norm(rhs))
    derivative_scale = max(float(np.linalg.norm(dT, 2)), np.finfo(float).tiny)
    solvability = float(abs(np.vdot(v0, rhs)) / derivative_scale)
    # Truncated SVD solve, with a gauge u0^H u1=0 at a one-dimensional kernel.
    inv = np.divide(1.0, s, out=np.zeros_like(s), where=s > threshold)
    u1 = Vh.conj().T @ (inv * (U.conj().T @ rhs))
    chain_residual = float(np.linalg.norm(T @ u1 - rhs) /
                           max(rhsnorm, derivative_scale, np.finfo(float).tiny))
    kernel_residual = float(np.linalg.norm(T @ u0) / max(scale, 1.0))
    chain_supported = bool(nullity == 1 and solvability <= relative_tol and
                           chain_residual <= relative_tol and kernel_residual <= relative_tol)
    obstruction = None
    order_two_supported = False
    if ddT is not None:
        rhs2 = dT @ u1 + ddT @ u0 / 2
        denom = max(float(np.linalg.norm(dT, 2) * np.linalg.norm(u1) +
                          np.linalg.norm(ddT, 2) / 2), np.finfo(float).tiny)
        obstruction = float(abs(np.vdot(v0, rhs2)) / denom)
        order_two_supported = chain_supported and obstruction > 100 * relative_tol
    return {"scope": "FINITE_ANALYTIC_MATRIX_PENCIL", "nullity": nullity,
            "geometric_multiplicity_one": nullity == 1,
            "kernel_residual": kernel_residual, "chain_residual": chain_residual,
            "solvability": solvability, "chain_supported": chain_supported,
            "order_two_obstruction": obstruction, "order_two_supported": bool(order_two_supported),
            "singular_values": s.tolist(), "continuum_certificate": False}


def puiseux_test(parameters: Any, gaps: Any, *, uncertainties: Any | None = None,
                 exponent_tol: float = 0.08, residual_tol: float = 0.05) -> dict:
    """Two shrinking-window fits; noisy/constant/linear gaps cannot pass."""
    t, y = np.asarray(parameters, dtype=float), np.asarray(gaps, dtype=float)
    if t.ndim != 1 or y.shape != t.shape or len(t) < 8:
        raise ValueError("need at least eight paired scalar samples")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)) or np.any(t <= 0) or np.any(y <= 0):
        raise ValueError("Puiseux samples must be finite and strictly positive")
    if len(np.unique(t)) != len(t) or max(t) / min(t) < 100:
        raise ValueError("need distinct perturbations spanning at least two decades")
    order = np.argsort(t)
    t, y = t[order], y[order]
    if uncertainties is not None:
        e = np.asarray(uncertainties, dtype=float)
        if e.shape != t.shape or not np.all(np.isfinite(e)) or np.any(e < 0):
            raise ValueError("invalid uncertainty estimates")
        if np.any(y <= 10 * e[order]):
            return {"supported": False, "reason": "splitting is not resolved above uncertainty"}
    windows = []
    for n in (len(t), max(4, len(t) // 2)):
        x, z = np.log(t[:n]), np.log(y[:n])
        slope, intercept = np.polyfit(x, z, 1)
        residual = float(np.sqrt(np.mean((z - slope * x - intercept) ** 2)))
        scores = {str(p): float(np.sqrt(np.mean((z - p*x - np.mean(z-p*x))**2)))
                  for p in (0.0, 0.5, 1.0)}
        windows.append({"n": n, "exponent": float(slope), "fit_residual": residual,
                        "fixed_exponent_scores": scores})
    supported = all(abs(w["exponent"] - 0.5) <= exponent_tol and
                    w["fit_residual"] <= residual_tol and
                    w["fixed_exponent_scores"]["0.5"] < min(w["fixed_exponent_scores"]["0.0"],
                                                             w["fixed_exponent_scores"]["1.0"])
                    for w in windows)
    return {"supported": bool(supported), "windows": windows,
            "scope": "DIRECTIONAL_NUMERICAL_FIT_NOT_A_THEOREM"}


@dataclass(frozen=True)
class SolverEvidence:
    solver: str
    discretization_lineage: str
    omega: complex
    residual: float
    uncertainty: float
    independently_located: bool
    converged: bool

    def usable(self, tol: float) -> bool:
        return bool(self.solver and self.discretization_lineage and self.independently_located and
                    self.converged and np.isfinite(self.omega) and np.isfinite(self.residual) and
                    0 <= self.residual <= tol and np.isfinite(self.uncertainty) and self.uncertainty >= 0)


@dataclass
class EP2Evidence:
    frequencies: list[SolverEvidence]
    partner_frequencies: list[SolverEvidence] = field(default_factory=list)
    matrix: dict | None = None
    puiseux: dict | None = None
    first_loop: list[int] | None = None
    second_loop: list[int] | None = None
    loops_converged: bool = False
    eigenvector_overlap: float | None = None
    eigenvector_trend_verified: bool = False
    coalescence_located: bool = False
    active_diagnostics: list[str] = field(default_factory=list)


def evaluate_ep2(evidence: EP2Evidence, *, residual_tol: float = 1e-8,
                 frequency_tol: float = 1e-7, gap_tol: float = 1e-6) -> dict:
    """All required witnesses must pass. Missing/error never means no EP exists."""
    if WITHDRAWN.intersection(evidence.active_diagnostics):
        raise ValueError("withdrawn diagnostics are forbidden in active EP evidence")
    reasons: list[str] = []
    good = [r for r in evidence.frequencies if r.usable(residual_tol)]
    if len(good) != len(evidence.frequencies) or len(good) < 2:
        reasons.append("insufficient converged independently located frequencies")
    if len({r.discretization_lineage for r in good}) < 2:
        reasons.append("requires two genuinely different discretization lineages")
    if good:
        spread = max(abs(r.omega - good[0].omega) for r in good)
        if spread > frequency_tol:
            reasons.append("independent solvers disagree")
    # A pointwise resolved pair is an exclusion at that point ONLY. Each branch
    # must itself have independent solver support; two CF calls do not suffice.
    partner = [r for r in evidence.partner_frequencies if r.usable(residual_tol)]
    partner_ok = (len(partner) == len(evidence.partner_frequencies) and len(partner) >= 2 and
                  len({r.discretization_lineage for r in partner}) >= 2)
    if not reasons and partner_ok:
        partner_spread = max(abs(r.omega-partner[0].omega) for r in partner)
        gap = abs(good[0].omega-partner[0].omega)
        uncertainty = max(r.uncertainty for r in good) + max(r.uncertainty for r in partner)
        if partner_spread <= frequency_tol and gap > max(gap_tol, 10*uncertainty):
            return {"verdict": Verdict.DISTINCT_AT_TESTED_POINT.value, "gap": float(gap),
                    "scope": "TESTED_POINT_ONLY", "continuum_certificate": False,
                    "note": "does not exclude collisions between grid points"}
    if not evidence.coalescence_located:
        reasons.append("coalescence parameters not independently located")
    matrix = evidence.matrix or {}
    if matrix.get("scope") != "FINITE_ANALYTIC_MATRIX_PENCIL" or not matrix.get("order_two_supported"):
        reasons.append("finite-pencil EP2 Jordan/geometric-multiplicity evidence missing")
    if not evidence.loops_converged or evidence.first_loop != [1, 0] or evidence.second_loop != [0, 1]:
        reasons.append("refined one-loop swap and two-loop restoration missing")
    if not (evidence.puiseux or {}).get("supported", False):
        reasons.append("shrinking-window Puiseux evidence missing")
    overlap = evidence.eigenvector_overlap
    if (overlap is None or not np.isfinite(overlap) or not 1-1e-5 <= overlap <= 1+1e-12 or
            not evidence.eigenvector_trend_verified):
        reasons.append("eigenvector collapse trend missing")
    return {"verdict": (Verdict.INCONCLUSIVE if reasons else Verdict.EP2_NUMERICALLY_SUPPORTED).value,
            "reasons": reasons, "scope": "FINITE_NUMERICAL_PROBLEM", "continuum_certificate": False}

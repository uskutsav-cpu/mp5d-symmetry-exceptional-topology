"""Convergence ladders whose missing rungs and numerical failures stay visible."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
from typing import Callable, Iterable

import numpy as np


@dataclass(frozen=True)
class Resolution:
    cf_depth: int = 200
    angular_n: int = 40
    radial_n: int = 160
    precision_dps: int = 16
    contour_length: float = 60.0
    scaling_angle_deg: float = 60.0
    continuation_step: float = 0.01

    def __post_init__(self):
        if min(self.cf_depth, self.angular_n, self.radial_n) < 4 or self.precision_dps < 15:
            raise ValueError("invalid precision or discretization size")
        if not 0 < self.scaling_angle_deg < 90 or self.contour_length <= 0 or self.continuation_step <= 0:
            raise ValueError("invalid contour or continuation settings")


@dataclass(frozen=True)
class Measurement:
    frequencies: tuple[complex, ...]
    residuals: tuple[float, ...]
    solver: str
    lineage: str
    effective_resolution: Resolution
    converged: bool
    uncertainty: tuple[float, ...] | None = None

    def to_dict(self):
        d = asdict(self)
        d["frequencies"] = [[z.real, z.imag] for z in self.frequencies]
        return d


def run_ladder(solve: Callable[[Resolution], Measurement], settings: Iterable[Resolution], *,
               frequency_tol: float = 1e-7, residual_tol: float = 1e-8,
               identity_radius: float = 0.1) -> dict:
    """Every rung independently invokes the solver; effective settings must match.

    The final three rungs must agree per fixed branch label and resolve any
    multi-branch gaps. The reported uncertainty is a heuristic based on drift,
    not a mathematically validated continuum error bound.
    """
    levels = list(settings)
    if len(levels) < 3 or len(set(levels)) != len(levels):
        raise ValueError("need at least three distinct resolutions")
    records: list[dict] = []
    good: list[Measurement] = []
    failures: list[str] = []
    for level in levels:
        try:
            result = solve(level)
            if result.effective_resolution != level:
                raise ValueError("solver did not actually use the requested settings")
            w, r = np.asarray(result.frequencies), np.asarray(result.residuals)
            if len(w) == 0 or r.shape != w.shape or not np.all(np.isfinite(w)) or not np.all(np.isfinite(r)):
                raise ValueError("invalid numerical outputs")
            if not result.converged or np.any(r < 0) or np.max(r) > residual_tol:
                raise ValueError("unconverged solution or unacceptable residual")
            if good:
                if len(w) != len(good[-1].frequencies):
                    raise ValueError("branch count changed")
                if result.lineage != good[-1].lineage or result.solver != good[-1].solver:
                    raise ValueError("a convergence ladder cannot switch solver implementations")
                if np.max(abs(w-np.asarray(good[-1].frequencies))) > identity_radius:
                    raise ValueError("possible branch jump across resolution ladder")
            record = result.to_dict()
            record["status"] = "MEASURED"
            records.append(record)
            good.append(result)
        except (ArithmeticError, ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            failures.append(str(exc))
            records.append({"requested_resolution": asdict(level), "status": "FAILED", "error": str(exc)})
    passed = not failures and len(good) == len(levels)
    max_drift = None
    if len(good) >= 3:
        tail = np.array([r.frequencies for r in good[-3:]])
        max_drift = float(np.max(abs(np.diff(tail, axis=0))))
        passed = passed and max_drift <= frequency_tol
        if len(tail[-1]) > 1:
            gaps = [abs(a-b) for i, a in enumerate(tail[-1]) for b in tail[-1][i+1:]]
            if min(gaps) <= max(10*max_drift, 10*frequency_tol):
                failures.append("distinct branches not resolved above ladder drift")
                passed = False
    return {"status": "PASS" if passed else "INCONCLUSIVE", "records": records,
            "failures": failures, "max_final_drift": max_drift,
            "scope": "NUMERICAL_CONVERGENCE_NOT_CONTINUUM_CERTIFICATION"}


def axis_ladders(base: Resolution, axes: dict[str, list]) -> dict[str, list[Resolution]]:
    """One-variable-at-a-time ladders; no axis is silently dropped."""
    out = {}
    for name, values in axes.items():
        if name not in asdict(base):
            raise ValueError(f"unknown resolution axis {name}")
        out[name] = [Resolution(**(asdict(base) | {name: v})) for v in values]
    return out


def robustness_labels(m1: int, m2: int, top_ell: int) -> list[dict]:
    """N=4,5 and next two allowed ell branches, explicitly out of formal scope."""
    if top_ell < abs(m1)+abs(m2) or (top_ell-abs(m1)-abs(m2)) % 2:
        raise ValueError("top ell is incompatible with the sector")
    return [{"m1": m1, "m2": m2, "ell": ell, "overtone": overtone,
             "scope": "OUT_OF_DOMAIN_ROBUSTNESS_PROBE"}
            for ell, overtone in product((top_ell, top_ell+2, top_ell+4), (4, 5))]

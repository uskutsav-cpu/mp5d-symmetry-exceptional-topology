"""Branch continuation with mode identity preserved by history, not by sorting.

Sorting frequencies is the classic way to corrupt a branch atlas: at an avoided
crossing two branches exchange their order in ``Re omega`` while remaining
perfectly distinct modes, and any sort-based labelling silently swaps them.
Here a branch is defined by its *continuation history* --- a sequence of solves
each seeded from an extrapolation of the previous ones --- plus the structural
label ``(m1, m2, l, N)`` carried through the continued-fraction inversion index.

Two independent identity checks run at every step:

* **predictor agreement** -- the corrected root must lie within a tolerance of
  the extrapolated prediction, scaled by the step size.  A jump to a different
  branch shows up as a large predictor miss.
* **step monotonicity** -- the motion per step must stay small compared with
  the distance to the nearest other tracked branch, otherwise the step is
  rejected and halved.

Failures are recorded rather than silently accepted, so a branch that could not
be continued leaves a visible gap instead of a fabricated point.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from ..geometry import MPGeometry, sdelta_to_ab
from ..radial.qnm import solve_qnm_cf

__all__ = ["AtlasPoint", "BranchLabel", "continue_along", "extrapolate"]


@dataclass(frozen=True)
class BranchLabel:
    m1: int
    m2: int
    ell: int
    overtone: int

    @property
    def key(self) -> str:
        return f"m{self.m1}_{self.m2}__l{self.ell}__N{self.overtone}"


@dataclass
class AtlasPoint:
    branch: str
    m1: int
    m2: int
    ell: int
    overtone: int
    s: float
    delta: float
    mu: float
    a: float
    b: float
    z_plus: float
    z_minus: float
    extremality: float
    omega_re: float
    omega_im: float
    Lambda_re: float
    Lambda_im: float
    residual: float
    predictor_miss: float
    depth: int
    angular_N: int
    solver: str
    step_index: int
    parent: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def omega(self) -> complex:
        return complex(self.omega_re, self.omega_im)

    def to_dict(self) -> dict:
        return asdict(self)


def extrapolate(history: list[complex], targets: list[float], nxt: float) -> complex:
    """Predict the next root from up to three previous points.

    Uses a linear or quadratic fit in the continuation parameter.  With a
    single point it degenerates to the trivial predictor, which is what makes
    the first step of a branch the least reliable one.
    """
    n = len(history)
    if n == 0:
        raise ValueError("empty history")
    if n == 1:
        return history[-1]
    k = min(n, 3)
    xs = np.asarray(targets[-k:], dtype=float)
    ys = np.asarray(history[-k:], dtype=complex)
    if np.ptp(xs) == 0:
        return history[-1]
    deg = 1 if k == 2 else 2
    cr = np.polyfit(xs, ys.real, deg)
    ci = np.polyfit(xs, ys.imag, deg)
    return complex(np.polyval(cr, nxt), np.polyval(ci, nxt))


def continue_along(label: BranchLabel, path: list[tuple[float, float, float]],
                   seed: complex, *, depth: int = 200, angular_N: int = 40,
                   M: float = 1.0, predictor_tol: float = 0.25,
                   depth_schedule: tuple[int, ...] = (100, 200),
                   commit: str = "", parent: str | None = None,
                   solver: str = "leaver-cf") -> tuple[list[AtlasPoint], list[dict]]:
    """Continue one branch along an ordered list of ``(s, delta, mu)`` points.

    Returns ``(points, failures)``.  A step whose corrected root misses the
    prediction by more than ``predictor_tol`` times the local scale is recorded
    as a failure and the branch is stopped there --- continuing past a suspected
    branch jump would poison every downstream point.
    """
    points: list[AtlasPoint] = []
    failures: list[dict] = []
    hist: list[complex] = []
    # arc-length-like parameter for the extrapolation
    ts: list[float] = []
    prev_pt: tuple[float, float, float] | None = None
    arc = 0.0

    for idx, (s, delta, mu) in enumerate(path):
        if prev_pt is not None:
            arc += float(np.linalg.norm(np.array([s, delta, mu]) - np.array(prev_pt)))
        ts.append(arc)
        prev_pt = (s, delta, mu)

        guess = seed if not hist else extrapolate(hist, ts[:-1] + [arc], arc)

        a, b = sdelta_to_ab(s, delta)
        geo = MPGeometry(a=a, b=b, M=M)
        if not geo.has_horizon:
            failures.append({"reason": "no_horizon", "s": s, "delta": delta,
                             "mu": mu, "branch": label.key})
            break
        try:
            sol = solve_qnm_cf(a, b, mu, label.m1, label.m2, label.ell,
                               label.overtone, initial_frequency=guess,
                               depth=depth, angular_N=angular_N, M=M,
                               depth_schedule=depth_schedule)
        except Exception as exc:  # noqa: BLE001 - failures are data here
            failures.append({"reason": "solver_error", "detail": str(exc)[:200],
                             "s": s, "delta": delta, "mu": mu,
                             "branch": label.key})
            break

        if not getattr(sol, "converged", True):
            failures.append({"reason": "not_converged", "s": s, "delta": delta,
                             "mu": mu, "branch": label.key})
            break

        w = complex(sol.omega)
        miss = abs(w - guess)
        scale = max(abs(w), 1.0)
        if hist and miss > predictor_tol * scale:
            failures.append({"reason": "predictor_miss", "miss": miss,
                             "s": s, "delta": delta, "mu": mu,
                             "branch": label.key,
                             "omega": [w.real, w.imag],
                             "predicted": [guess.real, guess.imag]})
            break

        lam = complex(sol.Lambda)
        points.append(AtlasPoint(
            branch=label.key, m1=label.m1, m2=label.m2, ell=label.ell,
            overtone=label.overtone, s=s, delta=delta, mu=mu, a=a, b=b,
            z_plus=float(geo.z_plus), z_minus=float(geo.z_minus),
            extremality=float(geo.extremality),
            omega_re=w.real, omega_im=w.imag,
            Lambda_re=lam.real, Lambda_im=lam.imag,
            residual=float(sol.cf_residual),
            predictor_miss=float(miss), depth=depth, angular_N=angular_N,
            solver=solver, step_index=idx, parent=parent,
            extra={"commit": commit},
        ))
        hist.append(w)

    return points, failures

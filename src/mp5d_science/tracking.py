"""Joint branch continuation with bijective assignment and explicit failures.

Frequency ordering is never a branch label. A small frequency gap also never
licenses accepting the same solver root twice. Adaptive retries are recorded.
All errors here are numerical/identity errors, not evidence against an EP.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment


class TrackingFailure(RuntimeError):
    """A requested branch identity could not be established."""


@dataclass(frozen=True)
class Mode:
    omega: complex
    residual: float
    vector: np.ndarray | None = field(default=None, repr=False, compare=False)
    uncertainty: float = 0.0
    converged: bool = True

    def validate(self, residual_tol: float) -> None:
        if not np.isfinite(self.omega):
            raise TrackingFailure("non-finite frequency")
        if not np.isfinite(self.residual) or not 0 <= self.residual <= residual_tol:
            raise TrackingFailure("unacceptable residual")
        if not self.converged:
            raise TrackingFailure("solver did not converge")
        if not np.isfinite(self.uncertainty) or self.uncertainty < 0:
            raise TrackingFailure("invalid uncertainty")
        if self.vector is not None:
            v = np.asarray(self.vector, dtype=complex)
            if v.ndim != 1 or not np.all(np.isfinite(v)) or np.linalg.norm(v) == 0:
                raise TrackingFailure("invalid eigenvector")


@dataclass(frozen=True)
class TrackingPolicy:
    residual_tol: float = 1e-8
    duplicate_atol: float = 1e-10
    predictor_atol: float = 1e-5
    predictor_rtol: float = 0.05
    max_motion_fraction: float = 0.45
    assignment_margin: float = 1e-5
    vector_weight: float = 0.5
    min_step: float = 1e-8
    max_halvings: int = 24
    max_evaluations: int = 10000

    def __post_init__(self) -> None:
        if not 0 < self.max_motion_fraction < 0.5:
            raise ValueError("motion fraction must lie in (0, 0.5)")
        for name in (
            "residual_tol",
            "duplicate_atol",
            "predictor_atol",
            "predictor_rtol",
            "assignment_margin",
            "min_step",
        ):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"invalid {name}")
        if self.max_halvings < 0 or self.max_evaluations < 1 or self.vector_weight < 0:
            raise ValueError("invalid tracking budget/weight")


def predict(
    history: Sequence[Sequence[complex]], times: Sequence[float], target: float
) -> np.ndarray:
    """Up-to-quadratic interpolation using ONLY previously accepted times.

    In particular, passing an extra target coordinate alongside the old roots
    is an error. This catches the previous off-by-one predictor bug.
    """
    if len(history) != len(times) or not history:
        raise ValueError("history and accepted times must have the same positive length")
    h = np.asarray(history, dtype=complex)
    t = np.asarray(times, dtype=float)
    if h.ndim != 2 or not np.all(np.isfinite(h)) or not np.all(np.isfinite(t)):
        raise ValueError("invalid predictor data")
    k = min(3, len(t))
    xs, ys = t[-k:], h[-k:]
    if len(np.unique(xs)) != k:
        raise ValueError("accepted predictor times must be distinct")
    if k == 1:
        return ys[-1].copy()
    # Centered Lagrange interpolation avoids polyfit conditioning at large arc length.
    out = np.zeros(h.shape[1], dtype=complex)
    for j in range(k):
        weight = 1.0
        for i in range(k):
            if i != j:
                weight *= (target - xs[i]) / (xs[j] - xs[i])
        out += weight * ys[j]
    return out


def _distinct(modes: Sequence[Mode], policy: TrackingPolicy) -> None:
    for i, a in enumerate(modes):
        for b in modes[i + 1 :]:
            floor = max(policy.duplicate_atol, a.uncertainty + b.uncertainty)
            if abs(a.omega - b.omega) <= floor:
                raise TrackingFailure(
                    "duplicate or unresolved roots: coalescence is NOT established"
                )


def _gap(modes: Sequence[Mode]) -> float:
    return min(
        (abs(a.omega - b.omega) for i, a in enumerate(modes) for b in modes[i + 1 :]),
        default=float("inf"),
    )


def assign(
    previous: Sequence[Mode],
    candidates: Sequence[Mode],
    predictions: Sequence[complex],
    policy: TrackingPolicy | None = None,
    *,
    enforce_motion: bool = True,
) -> tuple[list[Mode], dict]:
    if policy is None:
        policy = TrackingPolicy()

    """Globally minimum-cost one-to-one assignment, with ambiguity rejection.

    Compare the best matching with all alternatives formed by banning one of
    its chosen edges; this finds the exact second-best assignment cost. A
    frequency-only near-crossing that is undersampled is rejected, not guessed.
    """
    n = len(previous)
    if n == 0 or len(predictions) != n or len(candidates) != n:
        raise TrackingFailure("candidate count must equal the tracked branch count")
    for x in [*previous, *candidates]:
        x.validate(policy.residual_tol)
    _distinct(previous, policy)
    _distinct(candidates, policy)
    preds = np.asarray(predictions, dtype=complex)
    if not np.all(np.isfinite(preds)):
        raise TrackingFailure("non-finite prediction")
    scale = max(1.0, max(abs(x.omega) for x in previous))
    cost = np.abs(preds[:, None] - np.array([x.omega for x in candidates])[None, :]) / scale
    for i, a in enumerate(previous):
        for j, b in enumerate(candidates):
            if a.vector is not None and b.vector is not None:
                va, vb = np.asarray(a.vector), np.asarray(b.vector)
                if va.shape != vb.shape:
                    raise TrackingFailure("vectors need a common discretization for comparison")
                overlap = abs(np.vdot(va, vb)) / (np.linalg.norm(va) * np.linalg.norm(vb))
                cost[i, j] += policy.vector_weight * max(0.0, 1.0 - float(overlap))
    rows, cols = linear_sum_assignment(cost)
    best = float(cost[rows, cols].sum())
    second = float("inf")
    if n > 1:
        for i, j in zip(rows, cols, strict=True):
            alternate = cost.copy()
            alternate[i, j] = np.inf
            rr, cc = linear_sum_assignment(alternate)
            second = min(second, float(alternate[rr, cc].sum()))
        if second - best <= policy.assignment_margin:
            raise TrackingFailure("ambiguous bijective branch assignment; refine the path")
    ordered = [candidates[int(j)] for j in cols]
    misses = np.abs(np.array([m.omega for m in ordered]) - preds)
    if np.max(misses) > policy.predictor_atol + policy.predictor_rtol * scale:
        raise TrackingFailure("branch jump: predictor tolerance exceeded")
    motion = max(abs(b.omega - a.omega) for a, b in zip(previous, ordered, strict=True))
    separation = min(_gap(previous), _gap(ordered))
    if enforce_motion and motion >= policy.max_motion_fraction * separation:
        raise TrackingFailure("motion too large relative to branch gap; halve the step")
    return ordered, {
        "permutation": cols.tolist(),
        "predictor_miss": float(np.max(misses)),
        "max_motion": float(motion),
        "min_gap": float(separation),
        "assignment_margin": None if n == 1 else second - best,
    }


@dataclass
class TrackResult:
    parameters: list[list[float]]
    modes: list[list[Mode]]
    accepted_times: list[float]
    retries: list[dict]
    evaluations: int
    requested_vertices_completed: int
    complete: bool
    failure: str | None = None


RootSolver = Callable[[np.ndarray, Sequence[complex]], Sequence[Mode]]


def continue_branches(
    solver: RootSolver,
    path: Sequence[Sequence[float]],
    seeds: Sequence[Mode],
    policy: TrackingPolicy | None = None,
) -> TrackResult:
    if policy is None:
        policy = TrackingPolicy()

    """Traverse every requested segment, bisecting ambiguous steps in-place.

    ``solver`` must locate each candidate independently of branch assignment.
    Internal accepted subdivisions remain in the returned history. No file or
    global state is mutated. Exhausted budgets return an incomplete result.
    """
    vertices = np.asarray(path, dtype=float)
    if vertices.ndim != 2 or len(vertices) < 1 or not np.all(np.isfinite(vertices)):
        raise ValueError("path must be a nonempty finite array")
    if not seeds:
        raise ValueError("at least one seed is required")
    for mode in seeds:
        mode.validate(policy.residual_tol)
    _distinct(seeds, policy)
    params = [vertices[0].tolist()]
    modes = [list(seeds)]
    times = [0.0]
    retries: list[dict] = []
    evaluations = 0
    completed = 1
    failure: str | None = None
    for vertex in vertices[1:]:
        current = np.array(params[-1])
        distance = float(np.linalg.norm(vertex - current))
        if distance == 0:
            completed += 1
            continue
        # Endpoints in this stack are processed in geometric order.
        pending = [(vertex.copy(), 0)]
        while pending:
            target, level = pending.pop()
            current = np.array(params[-1])
            step = float(np.linalg.norm(target - current))
            if evaluations >= policy.max_evaluations:
                failure = "evaluation budget exhausted"
                break
            next_time = times[-1] + step
            predictions = predict([[m.omega for m in row] for row in modes], times, next_time)
            evaluations += 1
            try:
                candidates = list(solver(target.copy(), predictions.copy()))
                ordered, _ = assign(modes[-1], candidates, predictions, policy)
            except (TrackingFailure, ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
                retries.append(
                    {
                        "target": target.tolist(),
                        "step": step,
                        "level": level,
                        "reason": str(exc),
                        "evaluation": evaluations,
                    }
                )
                if step <= policy.min_step or level >= policy.max_halvings:
                    failure = f"unresolved branch identity: {exc}"
                    break
                pending.append((target, level + 1))
                pending.append(((current + target) / 2, level + 1))
                continue
            params.append(target.tolist())
            modes.append(ordered)
            times.append(next_time)
        if failure is not None:
            break
        completed += 1
    return TrackResult(
        params, modes, times, retries, evaluations, completed, failure is None, failure
    )


def loop_permutation(result: TrackResult, closure_tol: float = 1e-7) -> list[int]:
    """Require geometric loop closure and match endpoints bijectively."""
    if not result.complete:
        raise TrackingFailure("an incomplete loop has no established monodromy")
    if np.linalg.norm(np.array(result.parameters[-1]) - result.parameters[0]) > closure_tol:
        raise TrackingFailure("parameter path is not closed")
    start = np.array([x.omega for x in result.modes[0]])
    end = np.array([x.omega for x in result.modes[-1]])
    costs = abs(end[:, None] - start[None, :])
    rr, cc = linear_sum_assignment(costs)
    if np.max(costs[rr, cc]) > closure_tol:
        raise TrackingFailure("spectral endpoints do not close")
    return cc.tolist()

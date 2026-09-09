"""Bounded, reproducible stencil refinement of a resolved branch gap.

Local convergence is not a global lower bound. The same stencil also probes
multiple directions; independent continuation histories must be supplied by
the physical evaluator rather than fabricated from neighboring scalar gaps.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass(frozen=True)
class GapObservation:
    gap: float
    uncertainty: float
    branch_ids: tuple[str, str]
    independently_validated: bool
    continuation_directions: tuple[str, ...] = ()

    def validate(self, required_directions: int) -> None:
        if (
            not np.isfinite(self.gap)
            or self.gap <= 0
            or not np.isfinite(self.uncertainty)
            or self.uncertainty < 0
            or self.gap <= 10 * self.uncertainty
        ):
            raise ValueError("gap not resolved above numerical uncertainty")
        if len(set(self.branch_ids)) != 2 or not self.independently_validated:
            raise ValueError("branch identity or independent solver validation absent")
        if len(set(self.continuation_directions)) < required_directions:
            raise ValueError("independent continuation directions missing")


def refine_minimum(
    evaluate: Callable[[tuple[float, ...]], GapObservation],
    initial: Sequence[float],
    bounds: Sequence[tuple[float, float]],
    *,
    initial_radius: Sequence[float],
    parameter_tol: float = 1e-4,
    value_tol: float = 1e-7,
    max_levels: int = 15,
    max_evaluations: int = 10000,
    required_directions: int = 2,
) -> dict:
    x = np.asarray(initial, dtype=float)
    box = np.asarray(bounds, dtype=float)
    radius = np.asarray(initial_radius, dtype=float)
    if box.shape != (len(x), 2) or radius.shape != x.shape:
        raise ValueError("dimensions differ")
    if not np.all(np.isfinite(np.concatenate([x, box.ravel(), radius]))):
        raise ValueError("non-finite refinement configuration")
    if (
        np.any(box[:, 0] >= box[:, 1])
        or np.any(x < box[:, 0])
        or np.any(x > box[:, 1])
        or np.any(radius <= 0)
    ):
        raise ValueError("invalid refinement box or initial point")
    if parameter_tol <= 0 or value_tol <= 0 or max_levels < 2 or max_evaluations < 1:
        raise ValueError("invalid refinement tolerances/budget")
    cache: dict[tuple[float, ...], GapObservation] = {}
    failures = []
    history = []
    stable = 0
    status = "INCONCLUSIVE"
    previous = None
    identities = None
    calls = 0
    for level in range(max_levels):
        candidates = []
        level_failed = False
        for signs in product((-1, 0, 1), repeat=len(x)):
            p = tuple(float(t) for t in np.clip(x + radius * np.array(signs), box[:, 0], box[:, 1]))
            if p not in cache:
                if calls >= max_evaluations:
                    failures.append({"level": level, "reason": "evaluation budget exhausted"})
                    level_failed = True
                    break
                calls += 1
                try:
                    observation = evaluate(p)
                    observation.validate(required_directions)
                    if identities is None:
                        identities = observation.branch_ids
                    if observation.branch_ids != identities:
                        raise ValueError("branch identities changed inside refinement stencil")
                    cache[p] = observation
                except (ValueError, ArithmeticError, RuntimeError) as exc:
                    failures.append({"level": level, "point": list(p), "reason": str(exc)})
                    level_failed = True
                    continue
            if p in cache:
                candidates.append((p, cache[p]))
        if not candidates:
            break
        best_point, best = min(candidates, key=lambda pair: (pair[1].gap, pair[0]))
        new_x = np.array(best_point)
        movement = float(np.linalg.norm(new_x - x, ord=np.inf))
        change = None if previous is None else abs(best.gap - previous)
        history.append(
            {
                "level": level,
                "point": list(best_point),
                "gap": best.gap,
                "uncertainty": best.uncertainty,
                "radius": radius.tolist(),
                "movement": movement,
                "gap_change": change,
                "complete_stencil": not level_failed,
            }
        )
        if (
            not level_failed
            and change is not None
            and movement <= parameter_tol
            and change <= value_tol
            and max(radius) <= parameter_tol
        ):
            stable += 1
        else:
            stable = 0
        x, previous = new_x, best.gap
        if stable >= 2 and not failures:
            status = "LOCALLY_STABLE_NUMERICAL_MINIMUM"
            break
        if calls >= max_evaluations:
            break
        radius /= 2
    return {
        "status": status,
        "point": x.tolist(),
        "gap": previous,
        "levels": history,
        "failures": failures,
        "evaluations": calls,
        "cached_points": len(cache),
        "scope": "LOCAL_REFINEMENT_ONLY",
        "global_exclusion_proved": False,
        "boundary_active": [
            bool(abs(v - lo) <= parameter_tol or abs(v - hi) <= parameter_tol)
            for v, (lo, hi) in zip(x, box, strict=True)
        ],
    }

"""Physical adapters for continuation, seven-axis ladders and boundary refinement.

Every route starts from an explicitly identified anchor, not an arbitrary
frequency that is silently assigned an overtone label. Both radial solvers
must agree at the anchor before any adaptive search is allowed to start.
Numerical error estimates remain heuristic; none is a continuum enclosure.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal

import numpy as np

from .adaptive import GapObservation, refine_minimum
from .convergence import Measurement, Resolution, axis_ladders, run_ladder
from .physical import solve_a, solve_c
from .radial_polynomial import RadialParameters
from .tracking import Mode, TrackingPolicy, continue_branches


@dataclass(frozen=True)
class PairAnchor:
    point: tuple[float, float, float]
    m1: int
    m2: int
    ell: int
    labels: tuple[str, str]
    overtones: tuple[int, int]
    frequencies: tuple[complex, complex]

    def __post_init__(self):
        if len(self.frequencies) != 2 or len(self.labels) != 2 or len(self.overtones) != 2:
            raise ValueError("exactly two frequencies, labels and overtones required")
        if len(set(self.labels)) != 2 or len(set(self.overtones)) != 2:
            raise ValueError("two distinct, explicitly identified branches required")
        if min(self.overtones) < 0 or len(self.point) != 3 or not np.all(np.isfinite(self.point)):
            raise ValueError("invalid anchor labels or parameters")
        if (
            not all(np.isfinite(z) for z in self.frequencies)
            or abs(self.frequencies[0] - self.frequencies[1]) < 1e-8
        ):
            raise ValueError("anchor roots are nonfinite or unresolved")
        parameters(self, self.point)


def parameters(anchor: PairAnchor, point) -> RadialParameters:
    s, d, mu = map(lambda x: Decimal(str(x)), point)
    return RadialParameters(str(s + d), str(s - d), str(mu), anchor.m1, anchor.m2, anchor.ell)


def coordinate_route(start, end, order=(0, 1, 2), step=0.01):
    if tuple(sorted(order)) != (0, 1, 2) or not np.isfinite(step) or step <= 0:
        raise ValueError("route requires all three axes and positive finite step")
    current = np.array(start, dtype=float)
    end = np.array(end, dtype=float)
    if current.shape != (3,) or end.shape != (3,) or not np.all(np.isfinite([current, end])):
        raise ValueError("finite three-dimensional endpoints required")
    result = [tuple(current)]
    for axis in order:
        delta = end[axis] - current[axis]
        n = int(np.ceil(abs(delta) / step))
        beginning = current[axis]
        for j in range(1, n + 1):
            current = current.copy()
            current[axis] = beginning + delta * j / n
            result.append(tuple(current))
    return result


class PairEvaluator:
    def __init__(self, anchor: PairAnchor, resolution: Resolution, agreement_tol=2e-6):
        if resolution.precision_dps < 25:
            raise ValueError("physical A/C validation requires A working precision >=25 dps")
        self.anchor = anchor
        self.resolution = resolution
        self.agreement_tol = agreement_tol
        self.calls = []
        self.preflight_passed = False
        self.anchor_modes = []

    def pair(self, point, seeds):
        p = parameters(self.anchor, point)
        r = self.resolution
        primary = []
        checks = []
        for seed, N in zip(seeds, self.anchor.overtones, strict=True):
            A = solve_a(
                p, seed, overtone=N, depth=r.cf_depth, angular_n=r.angular_n, dps=r.precision_dps
            )
            C = solve_c(
                p,
                seed,
                radial_n=r.radial_n,
                angular_n=r.angular_n,
                length=r.contour_length,
                angle_deg=r.scaling_angle_deg,
            )
            agreement = abs(A.omega - C.omega)
            self.calls.append(
                {
                    "point": list(map(float, point)),
                    "overtone_label": N,
                    "A": A.to_dict(),
                    "C": C.to_dict(),
                    "agreement": agreement,
                }
            )
            if not A.converged or not C.converged or agreement > self.agreement_tol:
                raise ValueError(
                    "physical pair lacks same-root agreement and cross-inversion consistency"
                )
            primary.append(A)
            checks.append(C)
        gap = abs(primary[0].omega - primary[1].omega)
        error = max(abs(a.omega - c.omega) for a, c in zip(primary, checks, strict=True))
        if gap <= max(10 * error, 1e-7):
            raise ValueError("pair not resolved above numerical cross-check error")
        return primary, checks

    def preflight(self):
        primary, checks = self.pair(self.anchor.point, self.anchor.frequencies)
        displacement = max(
            abs(a.omega - w) for a, w in zip(primary, self.anchor.frequencies, strict=True)
        )
        if displacement > self.agreement_tol:
            raise ValueError("independent solver found different roots than the named anchor")
        self.anchor_modes = [
            Mode(a.omega, a.residual, uncertainty=abs(a.omega - c.omega))
            for a, c in zip(primary, checks, strict=True)
        ]
        self.preflight_passed = True
        return primary, checks

    def follow(self, path):
        if not self.preflight_passed:
            raise ValueError("anchor preflight has not passed")
        r = self.resolution

        def solve(point, predictions):
            p = parameters(self.anchor, point)
            modes = []
            for w, N in zip(predictions, self.anchor.overtones, strict=True):
                out = solve_a(
                    p, w, overtone=N, depth=r.cf_depth, angular_n=r.angular_n, dps=r.precision_dps
                )
                modes.append(Mode(out.omega, out.residual, converged=out.converged))
            return modes

        seeds = self.anchor_modes
        result = continue_branches(solve, path, seeds, policy=TrackingPolicy(max_evaluations=20000))
        if not result.complete:
            raise ValueError("continuation route failed before its requested endpoint")
        return result

    def observation(self, point):
        if not self.preflight_passed:
            self.preflight()
        results = []
        names = []
        # Distinct geometric detours, even when only one coordinate changes or
        # the target equals the anchor. Different names alone are not evidence
        # of different continuation histories.
        for axis, order in ((0, (0, 1, 2)), (1, (2, 1, 0))):
            waypoint = np.array(self.anchor.point, dtype=float)
            lo, hi = ((0.0, 0.42), (-0.2, 0.2))[axis]
            waypoint[axis] += 0.005 if waypoint[axis] + 0.005 <= hi else -0.005
            if not lo <= waypoint[axis] <= hi:
                raise ValueError("cannot construct independent continuation detour")
            path = coordinate_route(
                self.anchor.point, waypoint, order, self.resolution.continuation_step
            )
            path += coordinate_route(waypoint, point, order, self.resolution.continuation_step)[1:]
            tracked = self.follow(path)
            seeds = tuple(m.omega for m in tracked.modes[-1])
            a, c = self.pair(point, seeds)
            results.append((a, c))
            names.append("geometric-detour-axis-" + str(axis))
        first = results[0][0]
        drift = max(abs(results[1][0][i].omega - first[i].omega) for i in range(2))
        error = max(
            drift,
            *[abs(a.omega - c.omega) for aa, cc in results for a, c in zip(aa, cc, strict=True)],
            1e-12,
        )
        if drift > self.agreement_tol:
            raise ValueError("different continuation directions disagree")
        return GapObservation(
            abs(first[0].omega - first[1].omega), error, self.anchor.labels, True, tuple(names)
        )


def physical_ladders(anchor: PairAnchor, base: Resolution, axes: dict):
    reports = {}
    for axis, levels in axis_ladders(base, axes).items():

        def evaluate(level):
            evaluator = PairEvaluator(anchor, level)
            A, C = evaluator.preflight()
            if axis == "continuation_step":
                # Actually use the requested step on an in-domain excursion and
                # return. Repeating the same endpoint solve is not a step ladder.
                center = np.array(anchor.point)
                neighbor = center.copy()
                neighbor[0] = max(0.0, center[0] - 0.02)
                path = coordinate_route(center, neighbor, step=level.continuation_step)
                path += coordinate_route(neighbor, center, step=level.continuation_step)[1:]
                if len(path) < 3:
                    raise ValueError("nontrivial in-domain step ladder excursion required")
                tracked = evaluator.follow(path)
                A, C = evaluator.pair(anchor.point, tuple(m.omega for m in tracked.modes[-1]))
            return Measurement(
                tuple(a.omega for a in A),
                tuple(max(a.residual, c.residual) for a, c in zip(A, C, strict=True)),
                "A-direct+C-bordered",
                "two-independent-radial-discretizations",
                level,
                True,
                tuple(abs(a.omega - c.omega) for a, c in zip(A, C, strict=True)),
            )

        reports[axis] = run_ladder(evaluate, levels, frequency_tol=1e-6)
    return {
        "status": "PASS"
        if reports and all(r["status"] == "PASS" for r in reports.values())
        else "INCONCLUSIVE",
        "axes": reports,
        "scope": "POINTWISE_NUMERICAL_LADDERS_NOT_CONTINUUM_CERTIFICATION",
    }


def physical_boundary(anchor: PairAnchor, resolution: Resolution, plan: dict):
    evaluator = PairEvaluator(anchor, resolution)
    # Failure here is recorded as a failure; no scalar-gap surrogate is fitted.
    evaluator.preflight()
    result = refine_minimum(
        evaluator.observation,
        anchor.point,
        [(0, 0.42), (-0.2, 0.2), (0, 1.8)],
        initial_radius=plan["initial_radius"],
        parameter_tol=plan["parameter_tolerance"],
        value_tol=plan["gap_tolerance"],
        max_evaluations=plan["max_evaluations"],
    )
    result["physical_evaluations"] = evaluator.calls
    return result

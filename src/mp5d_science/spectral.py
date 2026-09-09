"""Numerical contours with an explicit, caller-audited analyticity contract.

Sampling does NOT prove analyticity or a rigorous root count. Meromorphic
continued fractions cannot be wrapped as analytic without a mathematical
regularization proof that controls all poles and branch cuts in the domain.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


class AnalyticityError(ValueError):
    pass


class ContourFailure(ArithmeticError):
    pass


@dataclass(frozen=True)
class AnalyticSpectralFunction:
    function: Callable[[complex], complex]
    center: complex
    radius: float
    justification: str
    representation: str
    known_singularities: tuple[complex, ...] = ()

    def __post_init__(self) -> None:
        if not np.isfinite(self.center) or not np.isfinite(self.radius) or self.radius <= 0:
            raise ValueError("invalid analytic domain")
        if not self.justification.strip():
            raise AnalyticityError("an explicit mathematical justification is required")
        if self.representation not in {
            "polynomial",
            "finite_analytic_determinant",
            "regularized_evans_function",
            "explicit_entire_function",
        }:
            raise AnalyticityError("uncontrolled meromorphic or nonanalytic representation")
        if any(
            not np.isfinite(s) or abs(s - self.center) <= self.radius
            for s in self.known_singularities
        ):
            raise AnalyticityError("declared analytic domain contains a known singularity")

    def __call__(self, z: complex) -> complex:
        if abs(z - self.center) > self.radius * (1 + 1e-13):
            raise AnalyticityError("evaluation outside the audited analytic domain")
        value = complex(self.function(z))
        if not np.isfinite(value):
            raise ContourFailure("non-finite spectral function")
        return value

    def require_disc(self, center: complex, radius: float) -> None:
        if not np.isfinite(radius) or radius <= 0 or not np.isfinite(center):
            raise ValueError("invalid contour")
        if abs(center - self.center) + radius > self.radius * (1 + 1e-13):
            raise AnalyticityError("contour interior is not covered by the analytic domain")


def polynomial_spectral(
    coefficients: list[complex], center: complex, radius: float
) -> AnalyticSpectralFunction:
    c = np.asarray(coefficients, dtype=complex)
    if c.ndim != 1 or len(c) < 2 or not np.all(np.isfinite(c)) or c[0] == 0:
        raise ValueError("finite descending polynomial coefficients required")
    return AnalyticSpectralFunction(
        lambda z: complex(np.polyval(c, z)),
        center,
        radius,
        "A finite polynomial is entire.",
        "polynomial",
    )


def finite_determinant(
    matrix: Callable[[complex], np.ndarray],
    center: complex,
    radius: float,
    entrywise_analyticity_justification: str,
    fixed_row_scale: np.ndarray | None = None,
) -> AnalyticSpectralFunction:
    """A determinant of an entrywise analytic matrix, with FIXED row scaling.

    Absolute-value row scaling depending on frequency is prohibited: it destroys
    analyticity. The matrix provider's entrywise analyticity remains an explicit
    proof obligation; a selected angular branch need not be globally analytic.
    """
    sample = np.asarray(matrix(center), dtype=complex)
    if sample.ndim != 2 or sample.shape[0] != sample.shape[1]:
        raise ValueError("matrix must be square")
    scale = (
        np.ones(sample.shape[0])
        if fixed_row_scale is None
        else np.asarray(fixed_row_scale, dtype=float)
    )
    if scale.shape != (sample.shape[0],) or not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError("row scales must be fixed positive finite values")

    def evaluate(z: complex) -> complex:
        mat = np.asarray(matrix(z), dtype=complex)
        if mat.shape != sample.shape or not np.all(np.isfinite(mat)):
            raise ContourFailure("matrix shape changed or entries became nonfinite")
        return complex(np.linalg.det(mat / scale[:, None]))

    return AnalyticSpectralFunction(
        evaluate, center, radius, entrywise_analyticity_justification, "finite_analytic_determinant"
    )


def _samples(f: AnalyticSpectralFunction, center: complex, radius: float, n: int):
    t = 2 * np.pi * np.arange(n) / n
    z = center + radius * np.exp(1j * t)
    values = np.array([f(complex(p)) for p in z])
    amplitudes = abs(values)
    if np.min(amplitudes) <= max(np.max(amplitudes) * 1e-13, np.finfo(float).tiny):
        raise ContourFailure("contour approaches a zero too closely or has excessive dynamic range")
    # Multiplying normalized phases avoids overflow/underflow in raw ratios.
    phases = values / amplitudes
    increments = np.angle(np.roll(phases, -1) * phases.conj())
    return z, values, increments


def count_zeros(
    f: AnalyticSpectralFunction,
    center: complex,
    radius: float,
    *,
    initial_n: int = 64,
    max_n: int = 16384,
    integer_tol: float = 1e-7,
    max_phase_step: float = np.pi / 3,
) -> dict:
    """Require three successive phase-resolved winding counts to agree.

    This is a convergence diagnostic, not an interval argument-principle proof.
    A deliberately false analyticity declaration remains a caller error.
    """
    if not isinstance(f, AnalyticSpectralFunction):
        raise AnalyticityError("raw callables/continued fractions are not accepted")
    f.require_disc(center, radius)
    if initial_n < 16 or max_n < 4 * initial_n:
        raise ValueError("budget must allow at least three contour resolutions")
    if not 0 < max_phase_step < np.pi or integer_tol <= 0:
        raise ValueError("invalid contour tolerances")
    records = []
    consecutive = 0
    prior = None
    n = initial_n
    while n <= max_n:
        _, _, increments = _samples(f, center, radius, n)
        winding = float(np.sum(increments) / (2 * np.pi))
        rounded = round(winding)
        resolved = (
            abs(winding - rounded) <= integer_tol
            and np.max(abs(increments)) < max_phase_step
            and rounded >= 0
        )
        records.append(
            {
                "n": n,
                "winding": winding,
                "max_phase_step": float(np.max(abs(increments))),
                "resolved": bool(resolved),
            }
        )
        consecutive = consecutive + 1 if resolved and prior == rounded else (1 if resolved else 0)
        if consecutive >= 3:
            return {
                "count": rounded,
                "converged": True,
                "ladder": records,
                "scope": "NUMERICAL_ARGUMENT_PRINCIPLE",
                "rigorous": False,
            }
        prior = rounded if resolved else None
        n *= 2
    raise ContourFailure(f"root count did not stabilize: {records}")


def roots_in_disc(
    f: AnalyticSpectralFunction,
    center: complex,
    radius: float,
    *,
    initial_n: int = 128,
    max_roots: int = 8,
    moment_tol: float = 1e-6,
) -> tuple[list[complex], dict]:
    """Converged contour moments in centered/scaled coordinates, then root check."""
    report = count_zeros(f, center, radius, initial_n=initial_n)
    count = report["count"]
    if count == 0:
        return [], report
    if count > max_roots:
        raise ContourFailure("contour encloses more roots than the configured bound")
    n = report["ladder"][-1]["n"]
    estimates = []
    for size in (n, 2 * n):
        z, values, _ = _samples(f, center, radius, size)
        t = 2 * np.pi * np.arange(size) / size
        frequencies = np.fft.fftfreq(size, 1 / size)
        derivative_t = np.fft.ifft(1j * frequencies * np.fft.fft(values))
        q = derivative_t / (1j * values)
        moments = [complex(np.mean(np.exp(1j * k * t) * q)) for k in range(count + 1)]
        if abs(moments[0] - count) > moment_tol:
            raise ContourFailure("complex moment count disagrees with winding count")
        elementary = [1 + 0j]
        for k in range(1, count + 1):
            elementary.append(
                sum((-1) ** (j - 1) * elementary[k - j] * moments[j] for j in range(1, k + 1)) / k
            )
        local_roots = np.roots([(-1) ** j * elementary[j] for j in range(count + 1)])
        estimates.append([complex(center + radius * r) for r in local_roots])
    from scipy.optimize import linear_sum_assignment

    costs = abs(np.array(estimates[0])[:, None] - np.array(estimates[1])[None, :])
    rr, cc = linear_sum_assignment(costs)
    discrepancy = float(np.max(costs[rr, cc]))
    if discrepancy > moment_tol * max(1, radius):
        raise ContourFailure("contour moment roots do not converge")
    roots = estimates[-1]
    _, boundary, _ = _samples(f, center, radius, 2 * n)
    residual_scale = max(abs(boundary))
    for root in roots:
        if abs(root - center) >= radius or abs(f(root)) > moment_tol * residual_scale:
            raise ContourFailure("moment reconstruction failed independent root residual check")
    report.update({"moment_resolution": 2 * n, "root_ladder_difference": discrepancy})
    return roots, report

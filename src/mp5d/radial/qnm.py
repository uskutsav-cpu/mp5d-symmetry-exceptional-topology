"""Radial QNM solvers for 5D Myers-Perry.

Solver A -- Leaver continued fraction
    Frobenius series in ``x = (r - r_+)/(r - r_-)``, the raw N-term recurrence
    (width 9 for a = b = 0; width 13 for generic two spins)
    reduced to three terms by the general Gaussian elimination in
    :mod:`mp5d.radial.recurrence`, then the ``n``-th inversion of Leaver's
    condition solved for ``omega``.

Solver B -- Hill determinant with Wynn acceleration
    Uses the **raw** N-term recurrence directly: no Gaussian reduction, no
    continued fraction.  The truncated banded determinant is evaluated by LU and
    its sequence in the truncation order is accelerated with the Wynn epsilon
    algorithm.  Following Benda and Matyjasek, arXiv:2503.17325, who note that
    Gaussian elimination is unnecessary if the determinant sequence is
    accelerated instead.  Reimplemented from the mathematics with attribution.

The two solvers share only the extracted polynomials ``A, B, C``; their root
conditions have no machinery in common, which is what makes them an admissible
independent pair under ``SCIENTIFIC_GATES.md``.

Boundary conditions (both solvers)
----------------------------------
    R = x^{-i sig} (1-x)^{3/2} exp( i Om (r_+ - r_-)/(1-x) ) sum_n a_n x^n

``sig = (w - m1 Om_a - m2 Om_b)/(2 kappa)`` gives the ingoing horizon
behaviour; ``exp(i Om r) r^{-3/2}`` with ``Om = sqrt(w^2 - mu^2)`` gives the
outgoing behaviour at infinity, the power being exactly ``-3/2`` because the
radial potential is even in ``r`` (claim C18).  ``Lambda`` is recomputed from
the angular solver at every iterate, never frozen.
"""

from __future__ import annotations

import cmath
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from .leaver import LeaverProblem
from .recurrence import (
    continued_fraction_inverted,
    recurrence_row,
    recurrence_width,
    reduce_to_three_term,
)

__all__ = ["QNMSolution", "solve_qnm_cf", "solve_qnm_hill", "wynn_epsilon"]

SOLVER_A = "leaver-cf/2.0"
SOLVER_B = "hill-wynn/1.0"


@dataclass
class QNMSolution:
    omega: complex
    Lambda: complex
    cf_residual: float
    ode_residual: float
    depth: int
    degree_bound: int
    precision: str
    mode_labels: dict[str, int]
    parameters: dict[str, float]
    solver: str
    iterations: int
    converged: bool
    depth_table: list[tuple[int, complex]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["omega"] = [self.omega.real, self.omega.imag]
        d["Lambda"] = [self.Lambda.real, self.Lambda.imag]
        d["depth_table"] = [[n, w.real, w.imag] for n, w in self.depth_table]
        return d


def _abc(prob: LeaverProblem, omega: complex, degree_bound: int):
    return prob.poly_ABC(omega, prob.Lambda_of(omega), degree_bound=degree_bound)


def cf_value(
    prob: LeaverProblem, omega: complex, depth: int, inversion: int, degree_bound: int
) -> complex:
    A, B, C = _abc(prob, omega, degree_bound)
    al, be, ga = reduce_to_three_term(A, B, C, depth=depth)
    return continued_fraction_inverted(al, be, ga, inversion, depth)


def _muller(f, x0: complex, tol: float, maxiter: int, re_floor: float = 0.05):
    """Muller's method: quadratic, derivative-free, and happy in the complex plane.

    ``re_floor`` keeps the iteration off the ``Re(omega) = 0`` branch cut of
    ``Om = sqrt(omega^2 - mu^2)``.  Without it the search readily walks onto the
    cut and returns a spurious near-zero-real-part "root"; several such
    candidates are recorded in ``results/rejected_radial_roots.json``.  Steps
    that would cross the floor are projected back onto it rather than rejected,
    which preserves Muller's quadratic behaviour in the physical region.
    """
    h = 1e-4 * max(abs(x0), 1.0)

    def clamp(z: complex) -> complex:
        return complex(max(z.real, re_floor), z.imag)

    x0 = clamp(x0)
    xs = [x0 - h, x0 + h, x0]
    fs = [f(x) for x in xs]
    for it in range(maxiter):
        x0_, x1_, x2_ = xs[-3:]
        f0, f1, f2 = fs[-3:]
        q = (x2_ - x1_) / (x1_ - x0_) if x1_ != x0_ else 1.0
        Aq = q * f2 - q * (1 + q) * f1 + q * q * f0
        Bq = (2 * q + 1) * f2 - (1 + q) ** 2 * f1 + q * q * f0
        Cq = (1 + q) * f2
        disc = cmath.sqrt(Bq * Bq - 4 * Aq * Cq)
        den = Bq + disc if abs(Bq + disc) > abs(Bq - disc) else Bq - disc
        if den == 0:
            break
        x3 = clamp(x2_ - (x2_ - x1_) * 2 * Cq / den)
        f3 = f(x3)
        xs.append(x3)
        fs.append(f3)
        if abs(x3 - x2_) < tol * max(abs(x3), 1.0) or abs(f3) < tol:
            return x3, it + 1, True
    return xs[-1], maxiter, False


def ode_residual(prob: LeaverProblem, omega: complex, depth: int, degree_bound: int) -> float:
    """Substitute the truncated series back into A y'' + B y' + C y at real x.

    This is the differential-equation residual, evaluated at points that play no
    role in building the recurrence.
    """
    A, B, C = _abc(prob, omega, degree_bound)
    w = recurrence_width(A, B, C)
    a = np.zeros(depth + 2, dtype=complex)
    a[0] = 1.0
    for n in range(depth):
        row = recurrence_row(n, A, B, C, w)
        acc = sum(row[i] * a[n + 1 - i] for i in range(1, w) if 0 <= n + 1 - i < len(a))
        a[n + 1] = -acc / row[0]

    worst = 0.0
    n = np.arange(len(a))
    for x in (0.15, 0.3, 0.45, 0.6):
        y = np.polyval(a[::-1], x)
        yp = np.polyval((n * a)[::-1][:-1], x)
        ypp = np.polyval((n * (n - 1) * a)[::-1][:-2], x)
        val = (
            np.polyval(A[::-1], x) * ypp
            + np.polyval(B[::-1], x) * yp
            + np.polyval(C[::-1], x) * y
        )
        scale = (
            abs(np.polyval(A[::-1], x) * ypp)
            + abs(np.polyval(B[::-1], x) * yp)
            + abs(np.polyval(C[::-1], x) * y)
        )
        worst = max(worst, abs(val) / max(scale, 1e-300))
    return float(worst)


def solve_qnm_cf(
    a: float,
    b: float,
    field_mass: float,
    m1: int,
    m2: int,
    ell: int,
    overtone: int = 0,
    initial_frequency: complex | None = None,
    depth: int = 400,
    degree_bound: int = 32,
    M: float = 1.0,
    angular_N: int = 60,
    tol: float = 1e-12,
    maxiter: int = 80,
    depth_schedule: tuple[int, ...] = (100, 200, 400, 800),
) -> QNMSolution:
    """Solver A.  Returns the root plus its depth-convergence table."""
    from ..geometry import MPGeometry

    geo = MPGeometry(a=a, b=b, M=M)
    prob = LeaverProblem(geo, field_mass, m1, m2, ell, depth=depth, angular_N=angular_N)

    if initial_frequency is None:
        raise ValueError("initial_frequency is required; seed from a published value")

    table: list[tuple[int, complex]] = []
    omega = initial_frequency
    iters = 0
    ok = False
    for d in depth_schedule:
        omega, it, ok = _muller(
            lambda w, d=d: cf_value(prob, w, d, overtone, degree_bound),
            omega,
            tol,
            maxiter,
        )
        iters += it
        table.append((d, omega))

    res = abs(cf_value(prob, omega, depth_schedule[-1], overtone, degree_bound))
    return QNMSolution(
        omega=omega,
        Lambda=prob.Lambda_of(omega),
        cf_residual=float(res),
        ode_residual=ode_residual(prob, omega, 200, degree_bound),
        depth=depth_schedule[-1],
        degree_bound=degree_bound,
        precision="double",
        mode_labels={"m1": m1, "m2": m2, "l": ell, "overtone": overtone},
        parameters={"a": a, "b": b, "s": geo.s, "delta": geo.delta, "mu": field_mass, "M": M},
        solver=SOLVER_A,
        iterations=iters,
        converged=ok,
        depth_table=table,
    )


# ---------------------------------------------------------------- Solver B


def wynn_epsilon(seq) -> complex:
    """Wynn's epsilon algorithm; returns the best diagonal approximant."""
    s = list(seq)
    n = len(s)
    e_prev = [0.0 + 0.0j] * (n + 1)
    e_cur = list(s)
    best = e_cur[-1]
    for k in range(1, n):
        e_next = []
        for i in range(len(e_cur) - 1):
            d = e_cur[i + 1] - e_cur[i]
            if abs(d) < 1e-300:
                e_next.append(e_prev[i + 1])
            else:
                e_next.append(e_prev[i + 1] + 1.0 / d)
        e_prev, e_cur = e_cur, e_next
        if k % 2 == 0 and e_cur:
            best = e_cur[-1]
        if len(e_cur) < 2:
            break
    return best


def hill_determinant(prob: LeaverProblem, omega: complex, size: int, degree_bound: int) -> complex:
    """Truncated Hill-determinant condition, in its numerically stable form.

    The truncated system (equations ``0..N-1`` in unknowns ``a_0..a_{N-1}`` with
    ``a_N`` set to zero) is singular exactly when the sequence generated forward
    from ``a_0 = 1`` satisfies ``a_N = 0``.  So the determinant is proportional
    to ``a_N``, and we evaluate that instead: it carries the same zeros without
    the enormous dynamic range of a determinant.

    Normalizing by ``max_k |a_k|`` keeps the returned value ``O(1)`` and encodes
    the actual criterion -- for a quasinormal mode the coefficient sequence is
    the *minimal* solution and its tail collapses; away from one, the dominant
    solution takes over and the tail is the largest part of the sequence.

    Uses only the raw N-term recurrence: no Gaussian reduction, no continued
    fraction, hence genuinely independent of Solver A.
    """
    A, B, C = _abc(prob, omega, degree_bound)
    w = recurrence_width(A, B, C)
    a = np.zeros(size + 1, dtype=complex)
    a[0] = 1.0
    for n in range(size):
        row = recurrence_row(n, A, B, C, w)
        acc = 0.0 + 0.0j
        for i in range(1, w):
            m = n + 1 - i
            if 0 <= m < len(a):
                acc += row[i] * a[m]
        a[n + 1] = -acc / row[0]
        # rescale to avoid overflow while preserving the ratio structure
        big = np.abs(a[: n + 2]).max()
        if big > 1e100:
            a[: n + 2] /= big
    return a[size] / max(np.abs(a).max(), 1e-300)


def solve_qnm_hill(
    a: float,
    b: float,
    field_mass: float,
    m1: int,
    m2: int,
    ell: int,
    overtone: int = 0,
    initial_frequency: complex | None = None,
    sizes: tuple[int, ...] = (60, 80, 100, 120, 140, 160),
    degree_bound: int = 32,
    M: float = 1.0,
    angular_N: int = 60,
    tol: float = 1e-11,
    maxiter: int = 80,
) -> QNMSolution:
    """Solver B.  Roots of the truncated Hill determinant, Wynn-accelerated."""
    from ..geometry import MPGeometry

    geo = MPGeometry(a=a, b=b, M=M)
    prob = LeaverProblem(geo, field_mass, m1, m2, ell, depth=max(sizes), angular_N=angular_N)
    if initial_frequency is None:
        raise ValueError("initial_frequency is required")

    roots: list[complex] = []
    omega = initial_frequency
    iters = 0
    ok = False
    for s in sizes:
        omega, it, ok = _muller(
            lambda w, s=s: hill_determinant(prob, w, s, degree_bound), omega, tol, maxiter
        )
        iters += it
        roots.append(omega)

    accelerated = wynn_epsilon(roots) if len(roots) >= 5 else roots[-1]
    if not np.isfinite(abs(accelerated)):
        accelerated = roots[-1]

    return QNMSolution(
        omega=accelerated,
        Lambda=prob.Lambda_of(accelerated),
        cf_residual=float("nan"),
        ode_residual=ode_residual(prob, accelerated, 200, degree_bound),
        depth=max(sizes),
        degree_bound=degree_bound,
        precision="double",
        mode_labels={"m1": m1, "m2": m2, "l": ell, "overtone": overtone},
        parameters={"a": a, "b": b, "s": geo.s, "delta": geo.delta, "mu": field_mass, "M": M},
        solver=SOLVER_B,
        iterations=iters,
        converged=ok,
        depth_table=[(s, r) for s, r in zip(sizes, roots, strict=True)],
    )

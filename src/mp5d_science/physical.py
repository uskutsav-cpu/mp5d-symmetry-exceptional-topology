"""New, explicitly versioned numerical solvers using the repository's equation.

A-direct uses algebraic (non-FFT) radial coefficients and true arbitrary
precision. C-bordered uses recurrence-free contour collocation and a fixed
bordered-system residual, rather than misidentifying a singular value as a
holomorphic residual. The formulations share the angular discretization and
physical equation, but not the radial discretization or root condition.

Neither solver establishes continuum convergence or an EP by itself.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass

import mpmath as mp
import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy import linalg, optimize

from .radial_polynomial import (
    COEFFICIENT_VERSION,
    RadialParameters,
    coefficient_formula,
    horizons,
    mul,
    power,
)


@dataclass
class PhysicalSolution:
    omega: complex
    residual: float
    converged: bool
    solver: str
    lineage: str
    settings: dict
    omega_decimal: tuple[str, str]
    angular_decimal: tuple[str, str]
    diagnostics: dict

    def to_dict(self):
        return {
            "omega": [self.omega.real, self.omega.imag],
            "residual": self.residual,
            "converged": self.converged,
            "solver": self.solver,
            "lineage": self.lineage,
            "settings": self.settings,
            "omega_decimal": list(self.omega_decimal),
            "angular_decimal": list(self.angular_decimal),
            "diagnostics": self.diagnostics,
        }


def angular_matrix(params: RadialParameters, omega: complex, n: int = 40) -> np.ndarray:
    if n <= params.angular_index or n < 4:
        raise ValueError("angular truncation too small for branch")
    alpha, beta = abs(params.m2), abs(params.m1)
    index = np.arange(n, dtype=float)
    z = 2 * index + alpha + beta
    diag = np.divide(beta * beta - alpha * alpha, z * (z + 2), out=np.zeros(n), where=z != 0)
    j = np.arange(n - 1, dtype=float)
    zoff = 2 * j + alpha + beta
    off = np.sqrt(
        4
        * (j + 1)
        * (j + alpha + beta + 1)
        * (j + 1 + alpha)
        * (j + 1 + beta)
        / ((zoff + 1) * (zoff + 2) ** 2 * (zoff + 3))
    )
    c2 = (omega * omega - float(params.mu) ** 2) * (float(params.a) ** 2 - float(params.b) ** 2)
    ell = 2 * index + alpha + beta
    return (
        np.diag(ell * (ell + 2) - c2 * (1 - diag) / 2)
        + np.diag(c2 * off / 2, 1)
        + np.diag(c2 * off / 2, -1)
    )


def angular_double(
    params: RadialParameters, omega: complex, n: int = 40, homotopy_steps: int = 12
) -> complex:
    """Track the selected angular branch by eigenvector overlap from c2=0."""
    if homotopy_steps < 2:
        raise ValueError("angular homotopy requires at least two steps")
    original = angular_matrix(params, omega, n)
    base = np.diag(
        [
            float(
                (2 * i + abs(params.m1) + abs(params.m2))
                * (2 * i + abs(params.m1) + abs(params.m2) + 2)
            )
            for i in range(n)
        ]
    )
    vector = np.eye(n, dtype=complex)[:, params.angular_index]
    val = complex(params.ell * (params.ell + 2))
    if float(params.a) ** 2 != float(params.b) ** 2:
        for step in range(1, homotopy_steps + 1):
            vals, vecs = linalg.eig(
                base + (original - base) * step / homotopy_steps, check_finite=False
            )
            overlaps = abs(vector.conj() @ vecs) / np.linalg.norm(vecs, axis=0)
            chosen = int(np.argmax(overlaps))
            # Ambiguous angular identity requires more homotopy resolution.
            sorted_overlap = np.sort(overlaps)
            if sorted_overlap[-1] - sorted_overlap[-2] < 1e-5:
                raise ValueError("angular branch assignment is ambiguous")
            vector, val = vecs[:, chosen], complex(vals[chosen])
    return val - (omega * omega - float(params.mu) ** 2) * float(params.b) ** 2


def angular_highprec(params: RadialParameters, omega, n: int):
    """Finite Jacobi determinant root with every coefficient evaluated at mp.dps.

    A double-precision continuation supplies only a starting guess. The final
    eigenvalue and characteristic residual are recomputed at working precision.
    """
    a, b, mu = map(mp.mpf, (params.a, params.b, params.mu))
    c2 = (omega * omega - mu * mu) * (a * a - b * b)
    if c2 == 0:
        return mp.mpf(params.ell * (params.ell + 2)) - (omega * omega - mu * mu) * b * b
    if n <= params.angular_index:
        raise ValueError("angular truncation excludes requested branch")
    alpha, beta = abs(params.m2), abs(params.m1)
    diag, off2 = [], []
    for j in range(n):
        z = 2 * j + alpha + beta
        xdiag = mp.mpf(beta * beta - alpha * alpha) / (z * (z + 2)) if z else mp.mpf(0)
        ell = 2 * j + alpha + beta
        diag.append(ell * (ell + 2) - c2 * (1 - xdiag) / 2)
        if j < n - 1:
            sq = mp.mpf(4 * (j + 1) * (j + alpha + beta + 1) * (j + 1 + alpha) * (j + 1 + beta)) / (
                (z + 1) * (z + 2) ** 2 * (z + 3)
            )
            off2.append(c2 * c2 * sq / 4)

    def determinant(value):
        old, current = (
            mp.mpc(1),
            (diag[0] - value) / max(1, abs(params.m1) + abs(params.m2) + 1) ** 2,
        )
        previous_scale = mp.mpf(max(1, abs(params.m1) + abs(params.m2) + 1) ** 2)
        for j in range(1, n):
            local_scale = mp.mpf((2 * j + alpha + beta + 1) ** 2)
            following = (
                (diag[j] - value) * current - off2[j - 1] * old / previous_scale
            ) / local_scale
            old, current, previous_scale = current, following, local_scale
        return current

    guess = mp.mpc(angular_double(params, complex(omega), n)) + (omega * omega - mu * mu) * b * b
    root = mp.findroot(
        determinant,
        (guess, guess + mp.mpf("0.0001")),
        solver="secant",
        tol=mp.power(10, -mp.mp.dps + 8),
        maxsteps=40,
        verify=True,
    )
    return root - (omega * omega - mu * mu) * b * b


def recurrence_cf(A, B, C, depth: int, inversion: int):
    """Finite 9-term recurrence and exact-width Gaussian reduction, in scalar arithmetic."""
    if depth < max(4, inversion + 2) or inversion < 0:
        raise ValueError("invalid continued-fraction depth/inversion")
    width = max(len(A) - 2, len(B) - 1, len(C)) + 1
    alpha, beta, gamma = [0] * depth, [0] * depth, [0] * depth

    def at(p, i):
        return p[i] if 0 <= i < len(p) else 0

    for n in range(depth):
        row = []
        for i in range(width):
            j = n + 1 - i
            row.append(at(A, i + 1) * j * (j - 1) + at(B, i) * j + at(C, i - 1) if j >= 0 else 0)
        for i in range(width - 1, 2, -1):
            if row[i] == 0:
                continue
            j = n + 2 - i
            if gamma[j] == 0:
                raise ZeroDivisionError("zero Gaussian pivot; do not silently regularize")
            f = row[i] / gamma[j]
            row[i - 2] -= f * alpha[j]
            row[i - 1] -= f * beta[j]
        alpha[n], beta[n], gamma[n] = row[:3]
    down = 0
    for j in range(depth - 1, inversion, -1):
        down = alpha[j - 1] * gamma[j] / (beta[j] - down)
    up = 0
    for j in range(inversion):
        up = alpha[j] * gamma[j + 1] / (beta[j] - up)
    return beta[inversion] - down - up


def solve_a(
    params: RadialParameters,
    seed: complex | tuple[str, str],
    *,
    overtone: int = 0,
    depth: int = 200,
    angular_n: int = 40,
    dps: int = 40,
    max_steps: int = 40,
) -> PhysicalSolution:
    """Solve a specified finite CF truncation; a separate ladder checks convergence."""
    if dps < 16:
        raise ValueError("working precision must be at least 16 decimal digits")
    with mp.workdps(dps):
        a, b, mu, M = map(mp.mpf, (params.a, params.b, params.mu, params.M))
        p, m = horizons(a, b, M, mp.sqrt)
        calls = 0

        def condition(omega):
            nonlocal calls
            calls += 1
            if omega.real <= 0 or omega.imag >= 0:
                raise ValueError("iteration left the configured damped positive-frequency sheet")
            lam = angular_highprec(params, omega, angular_n)
            k = mp.sqrt(omega * omega - mu * mu)
            if k.real < 0:
                k = -k
            A, B, C, _ = coefficient_formula(
                a, b, mu, params.m1, params.m2, omega, lam, p, m, k, mp.j
            )
            # Preserve the audited recurrence selection: removing a common infinity
            # factor changes this finite Gaussian-reduction root condition.
            # The exact degree-9/8/7 formula is lifted back to degree-13/12/11.
            t4 = power([1, -1], 4)
            return recurrence_cf(mul(A, t4), mul(B, t4), mul(C, t4), depth, overtone)

        x = mp.mpc(*seed) if isinstance(seed, tuple) else mp.mpc(seed)
        tolerance = mp.power(10, -dps + 8)
        # This perturbation forces a real solve instead of returning a seed solely
        # because the old solver's absolute residual happens to be small.
        root = mp.findroot(
            condition,
            (x - mp.mpf("0.00013"), x + mp.mpf("0.00017")),
            solver="secant",
            tol=tolerance,
            maxsteps=max_steps,
            verify=True,
        )
        residual = abs(condition(root))
        lam = angular_highprec(params, root, angular_n)
        k = mp.sqrt(root * root - mu * mu)
        if k.real < 0:
            k = -k
        A, B, C, _ = coefficient_formula(a, b, mu, params.m1, params.m2, root, lam, p, m, k, mp.j)
        lifted = [mul(v, power([1, -1], 4)) for v in (A, B, C)]
        cross = {
            str(i): abs(recurrence_cf(*lifted, depth, i))
            for i in range(min(depth - 2, max(4, overtone + 2)))
        }
        cross_tol = max(mp.mpf("1e-8"), mp.sqrt(tolerance))
        cross_ok = all(value < cross_tol for value in cross.values())
        # A single inverted CF can have an inversion-specific removable or
        # pivot-induced zero. It is not accepted merely because Newton converged.
        return PhysicalSolution(
            complex(root),
            float(residual),
            bool(residual < mp.sqrt(tolerance) and cross_ok),
            "A-direct-algebra/1.0",
            "frobenius-recurrence",
            {"cf_depth": depth, "angular_n": angular_n, "precision_dps": dps},
            (mp.nstr(root.real, dps), mp.nstr(root.imag, dps)),
            (mp.nstr(lam.real, dps), mp.nstr(lam.imag, dps)),
            {
                "function_evaluations": calls,
                "coefficient_version": COEFFICIENT_VERSION,
                "cross_inversion_residuals": {i: mp.nstr(v, dps) for i, v in cross.items()},
                "cross_inversion_consistent": cross_ok,
                "root_validity": "SUPPORTED_FINITE_CF_ROOT"
                if cross_ok
                else "INVERSION_SPECIFIC_ROOT_REQUIRES_REJECTION",
                "continuum_convergence_established": False,
            },
        )


def physical_ode(params: RadialParameters, omega: complex, lam: complex, r):
    """Untransformed radial equation coefficients, independent of polynomial code."""
    a, b, mu, M = map(float, (params.a, params.b, params.mu, params.M))
    r2 = r * r
    delta = r2 + a * a + b * b - M + a * a * b * b / r2
    delta_prime = 2 * r - 2 * a * a * b * b / (r2 * r)
    W = (
        (r2 + a * a) * (r2 + b * b) * omega
        - params.m1 * a * (r2 + b * b)
        - params.m2 * b * (r2 + a * a)
    )
    G = a * b * omega - a * params.m2 - b * params.m1
    V = (
        W * W / (r2 * r2 * delta)
        - G * G / r2
        - (a * a + b * b) * omega * omega
        + 2 * omega * (a * params.m1 + b * params.m2)
        - mu * mu * r2
        - lam
    )
    return delta, delta / r + delta_prime, V


class ContourProblem:
    """Recurrence-free Chebyshev collocation with a fixed, analytic border.

    The finite contour is valid only where Im(k*exp(i theta)) > 0. This
    formulation explicitly refuses quasiresonant points that need a near-90
    degree contour. No claim of solving the massive hyperboloidal problem.
    """

    def __init__(
        self,
        params: RadialParameters,
        seed: complex,
        *,
        radial_n: int = 160,
        angular_n: int = 40,
        length: float = 60,
        angle_deg: float = 60,
    ):
        if radial_n < 12 or length <= 0 or not 0 < angle_deg < 89:
            raise ValueError("invalid contour configuration")
        self.params, self.n, self.angular_n = params, radial_n, angular_n
        self.length, self.theta = float(length), np.deg2rad(angle_deg)
        a, b, M = map(float, (params.a, params.b, params.M))
        self.p, self.m = horizons(a, b, M, np.sqrt)
        nodes = np.cos(np.pi * (2 * np.arange(radial_n - 1) + 1) / (2 * (radial_n - 1)))
        self.rho = length * (nodes + 1) / 2
        self.V = cheb.chebvander(nodes, radial_n - 1)
        derivative = np.zeros((radial_n, radial_n))
        for j in range(radial_n):
            coeff = np.zeros(j + 1)
            coeff[-1] = 1
            dc = cheb.chebder(coeff)
            derivative[: len(dc), j] = dc
        self.D = self.V @ derivative * (2 / length)
        self.D2 = self.V @ (derivative @ derivative) * (2 / length) ** 2
        self.derivative_coefficients = derivative
        raw = self.raw_matrix(seed)
        self.row_scale = np.max(abs(raw), axis=1)
        if np.any(self.row_scale == 0):
            raise ValueError("zero collocation row")
        U, _, Vh = linalg.svd(raw / self.row_scale[:, None], check_finite=False)
        # Fixed border vectors: unlike re-evaluated SVD vectors they preserve
        # complex analyticity locally, until the bordered matrix becomes singular.
        self.border_left, self.border_right = U[:, -1], Vh[-1].conj()
        self.evaluations = 0

    def horizon_exponent(self, omega):
        a, b = float(self.params.a), float(self.params.b)
        p, m = self.p, self.m
        Wp = (
            (p * p + a * a) * (p * p + b * b) * omega
            - self.params.m1 * a * (p * p + b * b)
            - self.params.m2 * b * (p * p + a * a)
        )
        return -1j * Wp / (2 * p * (p * p - m * m))

    def check_radiation(self, omega):
        k = cmath.sqrt(omega * omega - float(self.params.mu) ** 2)
        if k.real < 0:
            k = -k
        rate = (k * np.exp(1j * self.theta)).imag
        if omega.real <= 0 or omega.imag >= 0 or rate <= 0:
            raise ValueError("radiation condition is not selected on this contour")
        if rate * self.length < 12:
            raise ValueError("far-end outgoing suppression too weak for this contour length")
        return rate

    def raw_matrix(self, omega):
        self.check_radiation(omega)
        lam = angular_double(self.params, omega, self.angular_n)
        phase = np.exp(1j * self.theta)
        r = self.p + self.rho * phase
        A, B, V = physical_ode(self.params, omega, lam, r)
        exponent = self.horizon_exponent(omega)
        C2 = A
        C1 = 2 * A * exponent / self.rho + B * phase
        C0 = (
            A * exponent * (exponent - 1) / self.rho**2
            + B * phase * exponent / self.rho
            + V * phase**2
        )
        interior = C2[:, None] * self.D2 + C1[:, None] * self.D + C0[:, None] * self.V
        return np.vstack([interior, np.ones(self.n, dtype=complex)])

    def matrix(self, omega):
        return self.raw_matrix(omega) / self.row_scale[:, None]

    def bordered_residual(self, omega):
        self.evaluations += 1
        mat = self.matrix(omega)
        border = np.zeros((self.n + 1, self.n + 1), dtype=complex)
        border[:-1, :-1] = mat
        border[:-1, -1] = self.border_left
        border[-1, :-1] = self.border_right.conj()
        rhs = np.zeros(self.n + 1, dtype=complex)
        rhs[-1] = 1
        solved = linalg.solve(border, rhs, assume_a="gen", check_finite=False)
        if not np.all(np.isfinite(solved)):
            raise ArithmeticError("singular or non-finite bordered solve")
        return complex(solved[-1])

    def offgrid_residual(self, omega, coefficients):
        # Independent collocation points, no recurrence and no reuse of ODE rows.
        nodes = np.linspace(-0.98, 0.98, 2 * self.n + 3)
        rho = self.length * (nodes + 1) / 2
        values = cheb.chebval(nodes, coefficients)
        first = cheb.chebval(nodes, cheb.chebder(coefficients)) * (2 / self.length)
        second = cheb.chebval(nodes, cheb.chebder(coefficients, 2)) * (2 / self.length) ** 2
        phase = np.exp(1j * self.theta)
        lam = angular_double(self.params, omega, self.angular_n)
        A, B, V = physical_ode(self.params, omega, lam, self.p + rho * phase)
        p = self.horizon_exponent(omega)
        terms = np.array(
            [
                A * second,
                (2 * A * p / rho + B * phase) * first,
                (A * p * (p - 1) / rho**2 + B * phase * p / rho + V * phase**2) * values,
            ]
        )
        scale = np.sum(abs(terms), axis=0)
        # Do not report relative noise at points where the entire mode is
        # exponentially suppressed. Also report the global normwise residual.
        normwise = float(np.linalg.norm(np.sum(terms, axis=0)) / max(np.linalg.norm(scale), 1e-300))
        mask = scale > np.max(scale) * 1e-9
        relative = float(np.max(abs(np.sum(terms, axis=0))[mask] / scale[mask]))
        return normwise, relative


def solve_c(
    params: RadialParameters,
    seed: complex,
    *,
    radial_n: int = 160,
    angular_n: int = 40,
    length: float = 60,
    angle_deg: float = 60,
    tol: float = 1e-10,
    max_evaluations: int = 100,
) -> PhysicalSolution:
    problem = ContourProblem(
        params, seed, radial_n=radial_n, angular_n=angular_n, length=length, angle_deg=angle_deg
    )

    def f(x):
        residual = problem.bordered_residual(complex(*x))
        return [residual.real, residual.imag]

    # A genuinely perturbed initial point; no equality-by-reusing-the-A-root.
    x0 = [seed.real + 1.7e-5, seed.imag - 1.3e-5]
    root = optimize.root(f, x0, method="hybr", options={"xtol": tol, "maxfev": max_evaluations})
    omega = complex(*root.x)
    residual = abs(problem.bordered_residual(omega))
    mat = problem.matrix(omega)
    _, singular, Vh = linalg.svd(mat, check_finite=False)
    normwise, pointwise = problem.offgrid_residual(omega, Vh[-1].conj())
    lam = angular_double(params, omega, angular_n)
    converged = bool(residual < tol and singular[-1] / singular[0] < tol and normwise < 1e-6)
    return PhysicalSolution(
        omega,
        float(residual),
        converged,
        "C-bordered-collocation/1.0",
        "recurrence-free-complex-contour",
        {
            "radial_n": radial_n,
            "angular_n": angular_n,
            "precision_dps": 16,
            "contour_length": length,
            "scaling_angle_deg": angle_deg,
        },
        (format(omega.real, ".17g"), format(omega.imag, ".17g")),
        (format(lam.real, ".17g"), format(lam.imag, ".17g")),
        {
            "function_evaluations": problem.evaluations,
            "optimizer_success": bool(root.success),
            "optimizer_message": str(root.message),
            "normalized_sigma_min": float(singular[-1] / singular[0]),
            "singular_value_gap": float(singular[-2] / max(singular[-1], np.finfo(float).tiny)),
            "offgrid_ode_residual_normwise": normwise,
            "offgrid_ode_residual_pointwise": pointwise,
            "outgoing_suppression_exponent": problem.check_radiation(omega) * length,
            "continuum_convergence_established": False,
        },
    )


def check_existing_cf_root(
    params: RadialParameters,
    omega: complex,
    *,
    depth=200,
    angular_n=40,
    overtone=0,
    dps=45,
    tolerance=1e-6,
) -> dict:
    """Guard a legacy double-precision root without changing its frequency.

    This is a root-validity diagnostic, not a rigorous normalized error bound.
    The looser tolerance accommodates the legacy floating-point input; a
    submission result must still pass the full high-precision A/C ladders.
    """
    if tolerance <= 0 or depth < 4 or angular_n < 4 or dps < 25:
        raise ValueError("invalid finite root-validation configuration")
    try:
        with mp.workdps(dps):
            w = mp.mpc(str(omega.real), str(omega.imag))
            if not mp.isfinite(w) or w.real <= 0 or w.imag >= 0:
                raise ValueError("frequency outside configured damped sheet")
            a, b, mu, M = map(mp.mpf, (params.a, params.b, params.mu, params.M))
            p, m = horizons(a, b, M, mp.sqrt)
            k = mp.sqrt(w * w - mu * mu)
            if k.real < 0:
                k = -k
            lam = angular_highprec(params, w, angular_n)
            A, B, C, _ = coefficient_formula(a, b, mu, params.m1, params.m2, w, lam, p, m, k, mp.j)
            ABC = [mul(poly, power([1, -1], 4)) for poly in (A, B, C)]
            residuals = {
                str(i): float(abs(recurrence_cf(*ABC, depth, i)))
                for i in range(min(depth - 2, max(4, overtone + 2)))
            }
            valid = all(np.isfinite(x) and x < tolerance for x in residuals.values())
        return {
            "valid": bool(valid),
            "cross_inversion_residuals": residuals,
            "tolerance": tolerance,
            "scope": "FINITE_CF_ROOT_CONSISTENCY_NOT_CONTINUUM_CERTIFICATION",
        }
    except (ArithmeticError, ValueError, np.linalg.LinAlgError) as exc:
        return {
            "valid": False,
            "reason": str(exc),
            "scope": "FINITE_CF_ROOT_CONSISTENCY_NOT_CONTINUUM_CERTIFICATION",
        }

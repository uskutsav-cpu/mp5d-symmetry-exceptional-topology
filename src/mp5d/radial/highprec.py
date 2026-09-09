"""Arbitrary-precision path: recurrence coefficients computed at working precision.

The double-precision path (:mod:`mp5d.radial.leaver`, :mod:`mp5d.radial.qnm`)
extracts the polynomial coefficients ``A, B, C`` with a NumPy FFT, which floors
every downstream residual near ``1e-12``.  Running only the *root solver* in
high precision would not help: the coefficients themselves would still be
double.  So this module recomputes the whole chain in ``mpmath``:

1. the ODE coefficient functions ``P2, P1, P0`` and the clearing factor, in
   ``mpmath.mpc``;
2. the polynomial coefficients, by an arbitrary-precision **DFT** on the unit
   circle with a half-sample phase offset (same contour as the double path, but
   evaluated as an explicit ``O(N^2)`` transform rather than an FFT);
3. the angular eigenvalue, from the Jacobi three-term recurrence via the same
   continued fraction -- *not* from a double-precision LAPACK eigensolve;
4. the Gaussian reduction, the continued fraction and Muller's method.

Everything is therefore consistent at the requested precision.
"""

from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from ..geometry import MPGeometry

__all__ = [
    "HighPrecisionProblem",
    "solve_qnm_mp",
    "MPSolution",
    "dft_polynomial_coefficients",
    "asymptotic_ratio",
]


def dft_polynomial_coefficients(fn, degree_bound: int, tol=None):
    """Arbitrary-precision recovery of polynomial coefficients on the unit circle.

    Mirrors :func:`mp5d.radial.leaver.polynomial_coefficients` exactly, including
    the half-sample phase offset that keeps every sample away from ``u = 1``, but
    computes the transform directly in ``mpmath``.
    """
    N = 1
    while N < 2 * (degree_bound + 1):
        N *= 2
    pts = [mp.expjpi(2 * mp.mpf(k + mp.mpf(1) / 2) / N) for k in range(N)]
    vals = [fn(u) for u in pts]

    coef = []
    for j in range(N):
        acc = mp.mpc(0)
        for k in range(N):
            acc += vals[k] * mp.expjpi(-2 * mp.mpf(j) * k / N)
        coef.append(acc / N * mp.expjpi(-mp.mpf(j) / N))

    scale = max(abs(c) for c in coef)
    if scale == 0:
        return [mp.mpc(0)] * (degree_bound + 1)
    if tol is None:
        tol = mp.mpf(10) ** (-(mp.mp.dps - 6))
    tail = max((abs(c) for c in coef[degree_bound + 1 :]), default=mp.mpf(0)) / scale
    if tail > tol:
        raise ValueError(
            f"coefficient tail {mp.nstr(tail, 5)} exceeds {mp.nstr(tol, 5)}: "
            "clearing factor incomplete or degree_bound too small"
        )
    return coef[: degree_bound + 1]


def _reduce_three_term_mp(A, B, C, depth: int):
    """Gaussian reduction, mpmath arithmetic.  Mirrors the double-precision one."""

    def g(arr, j):
        return arr[j] if 0 <= j < len(arr) else mp.mpc(0)

    # honest degrees
    def trimlen(arr):
        big = max((abs(x) for x in arr), default=mp.mpf(0))
        if big == 0:
            return 1
        thr = big * mp.mpf(10) ** (-(mp.mp.dps - 4))
        last = 0
        for i, x in enumerate(arr):
            if abs(x) > thr:
                last = i
        return last + 1

    width = max(trimlen(A) - 2, trimlen(B) - 1, trimlen(C))
    width = max(width, 1) + 1

    alpha = [mp.mpc(0)] * depth
    beta = [mp.mpc(0)] * depth
    gamma = [mp.mpc(0)] * depth

    for n in range(depth):
        row = []
        for i in range(width):
            m = n + 1 - i
            row.append(g(A, i + 1) * m * (m - 1) + g(B, i) * m + g(C, i - 1))
        for i in range(width):
            if n + 1 - i < 0:
                row[i] = mp.mpc(0)

        i = width - 1
        while i >= 3:
            if row[i] != 0:
                m = n + 1 - i
                j = m + 1
                if j < 0 or m < 0:
                    row[i] = mp.mpc(0)
                    i -= 1
                    continue
                piv = gamma[j]
                if piv == 0:
                    raise ZeroDivisionError(f"zero pivot gamma[{j}] at n={n}")
                f = row[i] / piv
                row[i - 2] -= f * alpha[j]
                row[i - 1] -= f * beta[j]
                row[i] = mp.mpc(0)
            i -= 1

        alpha[n], beta[n] = row[0], row[1]
        gamma[n] = row[2] if width > 2 else mp.mpc(0)
    return alpha, beta, gamma


def asymptotic_ratio(n, c, order: int = 1):
    """Minimal-solution ratio ``R_n = a_{n+1}/a_n`` at large ``n``.

    Structure (derived, see ``docs/RECURRENCE_ASYMPTOTICS.md``): the recurrence's
    characteristic equation is ``A(1/R) = 0`` and ``x = 1`` is a multiple root of
    ``A``, so ``R = 1`` is degenerate and the corrections come in powers of
    ``n^{-1/2}``:

        R_n = 1 + u1 n^{-1/2} + u2 n^{-1} + ...,   u1 = -sqrt(-2c)

    with ``c = i Omega (r_+ - r_-)`` the coefficient of the peeled exponential
    ``exp(c/(1-x))``.  The sign is fixed by minimality: the decaying solution is
    the one with ``Re(u1) < 0`` for the branch of ``sqrt`` used here.

    ``order = 0`` returns 1 (no tail).  Only ``order <= 1`` is supported: the
    higher coefficients have NOT been derived, and asking for them raises rather
    than silently returning a wrong tail.
    """
    if order <= 0:
        return mp.mpc(1)
    if order > 1:
        raise NotImplementedError(
            "only tail_order <= 1 is derived; u2 and beyond are not available"
        )
    u1 = -mp.sqrt(-2 * c)
    if u1.real > 0:
        u1 = -u1
    return 1 + u1 / mp.sqrt(n)


def _cf_inverted_mp(alpha, beta, gamma, inversion: int, depth: int, c=None, tail_order: int = 0):
    """Backward evaluation.  ``tail_order > 0`` closes the truncation with the
    asymptotic minimal-solution ratio instead of the crude ``frac = 0``.

    Terminal condition: the CF variable satisfies ``f_{k+1} = -alpha_k r_k``
    where ``r_k = a_{k+1}/a_k``, so seeding ``f_depth`` with the asymptotic
    ``r_{depth-1}`` is the exact analogue of Nollert's remainder estimate.
    """
    k = int(inversion)
    if tail_order > 0 and c is not None:
        R = asymptotic_ratio(depth - 1, c, tail_order)
        frac = -alpha[depth - 1] * R
    else:
        frac = mp.mpc(0)
    for j in range(depth - 1, k, -1):
        den = beta[j] - frac
        frac = alpha[j - 1] * gamma[j] / den
    down = frac
    up = mp.mpc(0)
    for j in range(0, k):
        den = beta[j] - up
        up = alpha[j] * gamma[j + 1] / den
    return beta[k] - down - up


def _muller_mp(f, x0, tol, maxiter: int, re_floor):
    h = mp.mpf("1e-4") * max(abs(x0), mp.mpf(1))

    def clamp(z):
        return mp.mpc(max(z.real, re_floor), z.imag)

    x0 = clamp(x0)
    xs = [x0 - h, x0 + h, x0]
    fs = [f(x) for x in xs]
    for it in range(maxiter):
        x0_, x1_, x2_ = xs[-3:]
        f0, f1, f2 = fs[-3:]
        q = (x2_ - x1_) / (x1_ - x0_)
        Aq = q * f2 - q * (1 + q) * f1 + q * q * f0
        Bq = (2 * q + 1) * f2 - (1 + q) ** 2 * f1 + q * q * f0
        Cq = (1 + q) * f2
        disc = mp.sqrt(Bq * Bq - 4 * Aq * Cq)
        den = Bq + disc if abs(Bq + disc) > abs(Bq - disc) else Bq - disc
        if den == 0:
            break
        x3 = clamp(x2_ - (x2_ - x1_) * 2 * Cq / den)
        f3 = f(x3)
        xs.append(x3)
        fs.append(f3)
        if abs(x3 - x2_) < tol * max(abs(x3), mp.mpf(1)) or abs(f3) < tol:
            return x3, it + 1, True
    return xs[-1], maxiter, False


@dataclass
class MPSolution:
    omega: complex
    Lambda: complex
    cf_residual: float
    dps: int
    depth: int
    converged: bool
    iterations: int
    #: full-precision decimal strings -- ``omega`` alone truncates to double and
    #: would hide exactly the digits this module exists to produce.
    omega_str: str = ""
    lambda_str: str = ""

    def omega_mp(self):
        with mp.workdps(self.dps):
            return mp.mpmathify(self.omega_str) if self.omega_str else mp.mpc(self.omega)


class HighPrecisionProblem:
    """Arbitrary-precision MP5D radial problem for one ``(m1, m2, l)`` sector."""

    def __init__(
        self,
        geo: MPGeometry,
        mu,
        m1: int,
        m2: int,
        ell: int,
        dps: int = 50,
        angular_depth: int = 200,
    ):
        self.geo = geo
        self.m1, self.m2, self.ell = m1, m2, ell
        self.dps = dps
        self.angular_depth = angular_depth
        rest = ell - abs(m1) - abs(m2)
        if rest < 0 or rest % 2 != 0:
            raise ValueError(f"l={ell} incompatible with (m1,m2)=({m1},{m2})")
        self.k_ang = rest // 2
        with mp.workdps(dps + 10):
            self.mu = mp.mpf(mu)
            self.a, self.b, self.M = mp.mpf(geo.a), mp.mpf(geo.b), mp.mpf(geo.M)
            # Horizons must be recomputed HERE at working precision.  Reading
            # geo.z_plus would import a double-precision root and silently cap
            # the whole chain near 16 digits, which is exactly the leak this
            # module exists to remove.
            c = self.M - self.a**2 - self.b**2
            disc = c**2 - 4 * self.a**2 * self.b**2
            sq = mp.sqrt(disc) if disc > 0 else mp.mpf(0)
            self.zp = (c + sq) / 2
            self.zm = (c - sq) / 2
            if self.zm < 0:
                self.zm = mp.mpf(0)
            self.rp = mp.sqrt(self.zp)
            self.rm = mp.sqrt(self.zm)

    # -- angular eigenvalue, high precision, via the Jacobi 3-term CF --------
    def Ahat(self, omega):
        """``Ahat`` from the Jacobi recurrence, solved by continued fraction.

        Deliberately *not* a LAPACK eigensolve: that would silently reintroduce
        double precision into the high-precision chain.
        """
        al, be = abs(self.m2), abs(self.m1)
        L = al + be
        c2 = (omega**2 - self.mu**2) * (self.a**2 - self.b**2)
        target = mp.mpf(self.ell) * (self.ell + 2)
        if c2 == 0:
            return mp.mpc(target)

        def cf(Ah):
            A = [mp.mpc(0), mp.mpc(1), mp.mpc(-1)]
            B = [mp.mpc(al + 1), mp.mpc(-(al + be + 2))]
            C = [(Ah - L * (L + 2)) / 4, c2 / 4]
            a_, b_, g_ = _reduce_three_term_mp(A, B, C, self.angular_depth)
            return _cf_inverted_mp(a_, b_, g_, self.k_ang, self.angular_depth)

        tol = mp.mpf(10) ** (-(self.dps - 8))
        root, _, _ = _muller_mp(cf, mp.mpc(target), tol, 60, mp.mpf("-1e300"))
        return root

    def Lambda(self, omega):
        return self.Ahat(omega) - (omega**2 - self.mu**2) * self.b**2

    # -- radial ODE coefficients in the Leaver variable ---------------------
    def _P(self, u, omega, Lam):
        a, b, M, mu = self.a, self.b, self.M, self.mu
        m1, m2 = self.m1, self.m2
        rp, rm = self.rp, self.rm
        d = rp - rm

        r = (rp - rm * u) / (1 - u)
        r2 = r * r
        q = a**2 + b**2 - M
        w = a**2 * b**2

        Delta = r2 + q + w / r2
        dDelta = 2 * r - 2 * w / (r2 * r)
        W = (r2 + a**2) * (r2 + b**2) * omega - m1 * a * (r2 + b**2) - m2 * b * (r2 + a**2)
        G = a * b * omega - a * m2 - b * m1
        V = (
            W**2 / (r2 * r2 * Delta)
            - G**2 / r2
            - (a**2 + b**2) * omega**2
            + 2 * omega * (a * m1 + b * m2)
            - mu**2 * r2
            - Lam
        )

        dudr = (1 - u) ** 2 / d
        d2udr2 = -2 * (1 - u) ** 3 / d**2
        A_u = Delta * dudr**2
        B_u = Delta * d2udr2 + (Delta / r + dDelta) * dudr

        Om = mp.sqrt(omega**2 - mu**2)
        if Om.real < 0:
            Om = -Om
        zp = self.zp
        Om_a = a / (zp + a**2)
        Om_b = b / (zp + b**2)
        kappa = rp * (zp - self.zm) / ((zp + a**2) * (zp + b**2))
        sig = (omega - m1 * Om_a - m2 * Om_b) / (2 * kappa)
        c = mp.mpc(0, 1) * Om * d

        LF = -mp.mpc(0, 1) * sig / u - mp.mpf("1.5") / (1 - u) + c / (1 - u) ** 2
        LFp = mp.mpc(0, 1) * sig / u**2 - mp.mpf("1.5") / (1 - u) ** 2 + 2 * c / (1 - u) ** 3
        FppF = LFp + LF * LF
        return (A_u, 2 * A_u * LF + B_u, A_u * FppF + B_u * LF + V)

    def _clearing(self, u):
        a, b, M = self.a, self.b, self.M
        rp, rm = self.rp, self.rm
        Nr = rp - rm * u
        r = Nr / (1 - u)
        r2 = r * r
        q = a**2 + b**2 - M
        w = a**2 * b**2
        Delta_num = (r2 * r2 + q * r2 + w) * (1 - u) ** 4
        return u**2 * (1 - u) ** 4 * Delta_num * Nr**4

    def poly_ABC(self, omega, Lam, degree_bound: int = 32):
        out = []
        for idx in range(3):

            def f(u, idx=idx):
                return self._P(u, omega, Lam)[idx] * self._clearing(u)

            out.append(dft_polynomial_coefficients(f, degree_bound))

        # strip the common power of u introduced by the clearing factor
        def val(arr):
            big = max(abs(x) for x in arr)
            thr = big * mp.mpf(10) ** (-(self.dps - 4))
            for i, x in enumerate(arr):
                if abs(x) > thr:
                    return i
            return len(arr)

        k = min(val(a) for a in out)
        return [a[k:] for a in out] if k else out

    def exp_coefficient(self, omega):
        """``c = i Omega (r_+ - r_-)`` -- the peeled exponential's coefficient."""
        Om = mp.sqrt(omega**2 - self.mu**2)
        if Om.real < 0:
            Om = -Om
        return mp.mpc(0, 1) * Om * (self.rp - self.rm)

    def cf_value(
        self, omega, depth: int, inversion: int, degree_bound: int = 32, tail_order: int = 0
    ):
        Lam = self.Lambda(omega)
        A, B, C = self.poly_ABC(omega, Lam, degree_bound)
        al, be, ga = _reduce_three_term_mp(A, B, C, depth)
        c = self.exp_coefficient(omega) if tail_order > 0 else None
        return _cf_inverted_mp(al, be, ga, inversion, depth, c, tail_order)


def solve_qnm_mp(
    a,
    b,
    mu,
    m1,
    m2,
    ell,
    overtone=0,
    initial_frequency=None,
    dps: int = 50,
    depth: int = 300,
    degree_bound: int = 32,
    M: float = 1.0,
    maxiter: int = 60,
    tail_order: int = 0,
) -> MPSolution:
    """Solve one QNM entirely at ``dps`` decimal digits."""
    if initial_frequency is None:
        raise ValueError("initial_frequency is required")
    geo = MPGeometry(a=a, b=b, M=M)
    with mp.workdps(dps):
        prob = HighPrecisionProblem(geo, mu, m1, m2, ell, dps=dps)
        tol = mp.mpf(10) ** (-(dps - 8))
        w0 = mp.mpc(initial_frequency.real, initial_frequency.imag)
        root, iters, ok = _muller_mp(
            lambda w: prob.cf_value(w, depth, overtone, degree_bound, tail_order),
            w0,
            tol,
            maxiter,
            mp.mpf("0.05"),
        )
        res = abs(prob.cf_value(root, depth, overtone, degree_bound, tail_order))
        lam = prob.Lambda(root)
        return MPSolution(
            omega=complex(root),
            Lambda=complex(lam),
            cf_residual=float(res),
            dps=dps,
            depth=depth,
            converged=ok,
            iterations=iters,
            omega_str=mp.nstr(root, dps, strip_zeros=False),
            lambda_str=mp.nstr(lam, dps, strip_zeros=False),
        )

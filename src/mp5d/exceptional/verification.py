"""EP verification: root tracking, monodromy, Puiseux fitting, Jordan chains.

Every routine here is written against an *abstract* spectral condition -- a
callable ``F(omega, p)`` analytic in ``omega`` -- so that the whole machinery can
be validated on synthetic problems whose exceptional structure is known in
closed form before it is pointed at the Myers-Perry spectrum.  That ordering
matters: a detector that has never been shown to fire on a real EP cannot
support a negative result.

The three tests implemented are the ones that actually distinguish an EP2 from
an avoided crossing:

* **zero count** -- two zeros inside a contour that shrinks with the candidate
  (argument principle, no root finder involved);
* **monodromy** -- transporting the roots around a loop in parameter space
  exchanges them, and a second loop restores them;
* **Puiseux** -- the splitting grows like ``t^{1/2}``, not like ``t``.

An avoided crossing passes none of these.  A numerical artifact passes none of
them reproducibly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "MonodromyResult",
    "PuiseuxFit",
    "roots_in_disc",
    "track_roots_on_loop",
    "monodromy",
    "fit_puiseux",
    "jordan_chain_matrix",
]


def roots_in_disc(f, center: complex, radius: float, n: int = 256,
                  max_roots: int = 4) -> list[complex]:
    """Locate all zeros of ``f`` in a disc using contour moments.

    Uses the Kravanja/Van Barel approach: the moments
    ``s_k = (1/2 pi i) oint z^k f'/f dz`` are the power sums of the enclosed
    roots, and Newton's identities convert power sums to a monic polynomial
    whose roots are the enclosed zeros.  ``f'`` is obtained from the same
    circle by Cauchy differentiation, so only ``f`` samples are needed.

    This is derivative-free in the user's sense (no analytic ``f'`` required)
    and does not depend on any starting guess, which makes it independent
    evidence relative to a Newton search.
    """
    t = 2.0 * np.pi * np.arange(n) / n
    pts = center + radius * np.exp(1j * t)
    vals = np.array([complex(f(p)) for p in pts], dtype=complex)
    if np.any(vals == 0):
        raise FloatingPointError("f vanishes on the contour")

    # f' on the contour by spectral differentiation of the periodic samples.
    # f(center + R e^{it}) as a function of t has derivative i R e^{it} f'.
    spec = np.fft.fft(vals)
    freqs = np.fft.fftfreq(n, d=1.0 / n)
    dvals_dt = np.fft.ifft(1j * freqs * spec)
    fprime = dvals_dt / (1j * radius * np.exp(1j * t))

    ratio = fprime / vals
    dz_dt = 1j * radius * np.exp(1j * t)
    # s_k = (1/2 pi i) * integral over t of (z^k f'/f) dz/dt dt
    moments = []
    for k in range(2 * max_roots + 1):
        integrand = (pts ** k) * ratio * dz_dt
        moments.append(np.sum(integrand) / n / (1j))  # (2pi/n) * .../(2 pi i)

    count = moments[0].real
    m = round(count)
    if abs(count - m) > 0.1:
        raise ValueError(f"non-integer zero count {count:.4f} in disc")
    if m <= 0:
        return []
    if m > max_roots:
        raise ValueError(f"{m} roots exceeds max_roots={max_roots}")

    # Newton's identities: power sums p_k -> elementary symmetric e_k
    p = [complex(moments[k]) for k in range(m + 1)]
    e = [1.0 + 0j]
    for k in range(1, m + 1):
        acc = 0.0 + 0j
        for i in range(1, k + 1):
            acc += ((-1) ** (i - 1)) * e[k - i] * p[i]
        e.append(acc / k)
    # monic polynomial z^m - e1 z^{m-1} + e2 z^{m-2} - ...
    coeffs = [((-1) ** k) * e[k] for k in range(m + 1)]
    return [complex(z) for z in np.roots(coeffs)]


@dataclass
class MonodromyResult:
    n_points: int
    start: list[complex]
    end: list[complex]
    permutation: list[int]
    is_swap: bool
    max_step: float
    closed_error: float
    detail: dict = field(default_factory=dict)


def track_roots_on_loop(root_fn, loop_points, n_roots: int = 2) -> list[list[complex]]:
    """Follow ``n_roots`` roots along a closed parameter loop.

    ``root_fn(p)`` must return the list of roots near the candidate for
    parameter ``p``.  Roots are matched between consecutive parameter values by
    nearest-neighbour assignment, which is valid provided the loop is sampled
    finely enough that the motion per step is small compared with the root
    separation -- checked and reported by :func:`monodromy`.
    """
    history: list[list[complex]] = []
    prev: list[complex] | None = None
    for p in loop_points:
        rts = list(root_fn(p))
        if len(rts) < n_roots:
            raise ValueError(f"only {len(rts)} roots found at parameter {p}")
        rts = rts[:n_roots] if prev is None else rts
        if prev is None:
            ordered = sorted(rts, key=lambda z: (z.real, z.imag))[:n_roots]
        else:
            # greedy nearest-neighbour matching against previous ordering
            pool = list(rts)
            ordered = []
            for q in prev:
                j = min(range(len(pool)), key=lambda i: abs(pool[i] - q))
                ordered.append(pool.pop(j))
        history.append(ordered)
        prev = ordered
    return history


def monodromy(root_fn, center_params, radius: float, n_points: int = 48,
              n_roots: int = 2, plane=(0, 1)) -> MonodromyResult:
    """Transport roots around a circle in a 2-plane of parameter space.

    ``center_params`` is a real vector; ``plane`` names the two components that
    are varied.  Returns the permutation induced on the tracked roots.  For a
    genuine EP2 enclosed by the loop the permutation is a transposition.
    """
    c = np.asarray(center_params, dtype=float)
    i, j = plane
    loop = []
    for k in range(n_points + 1):
        th = 2.0 * np.pi * k / n_points
        p = c.copy()
        p[i] += radius * np.cos(th)
        p[j] += radius * np.sin(th)
        loop.append(p)

    hist = track_roots_on_loop(root_fn, loop, n_roots=n_roots)
    start, end = hist[0], hist[-1]

    steps = [
        max(abs(hist[k + 1][r] - hist[k][r]) for r in range(n_roots))
        for k in range(len(hist) - 1)
    ]
    max_step = max(steps) if steps else 0.0

    perm = []
    for e in end:
        perm.append(min(range(n_roots), key=lambda r: abs(start[r] - e)))
    closed = max(abs(end[k] - start[perm[k]]) for k in range(n_roots))
    is_swap = (n_roots == 2 and perm == [1, 0])

    return MonodromyResult(
        n_points=n_points, start=start, end=end, permutation=perm,
        is_swap=is_swap, max_step=float(max_step), closed_error=float(closed),
        detail={"radius": radius, "plane": list(plane)},
    )


@dataclass
class PuiseuxFit:
    exponent: float
    residual_sqrt: float
    residual_linear: float
    prefers_sqrt: bool
    coefficient: complex
    n_points: int


def fit_puiseux(ts, splittings) -> PuiseuxFit:
    """Fit ``|omega_+ - omega_-| ~ C t^p`` and test ``p = 1/2`` against ``p = 1``.

    At an EP2 the two branches separate like ``t^{1/2}``; at an avoided crossing
    the separation tends to a nonzero constant, and at an ordinary (transversal)
    crossing it is linear.  The exponent is obtained by least squares in log-log
    space, and both hypotheses are scored so the comparison is explicit rather
    than eyeballed.
    """
    t = np.asarray(ts, dtype=float)
    d = np.asarray([abs(x) for x in splittings], dtype=float)
    ok = (t > 0) & (d > 0)
    t, d = t[ok], d[ok]
    if len(t) < 3:
        raise ValueError("need at least 3 positive samples")
    A = np.vstack([np.log(t), np.ones_like(t)]).T
    sol, *_ = np.linalg.lstsq(A, np.log(d), rcond=None)
    p, logC = sol

    def resid(power: float) -> float:
        c = np.exp(np.mean(np.log(d) - power * np.log(t)))
        return float(np.sqrt(np.mean((np.log(d) - (power * np.log(t) + np.log(c))) ** 2)))

    r_half, r_one = resid(0.5), resid(1.0)
    return PuiseuxFit(
        exponent=float(p), residual_sqrt=r_half, residual_linear=r_one,
        prefers_sqrt=bool(r_half < r_one), coefficient=complex(np.exp(logC)),
        n_points=int(len(t)),
    )


def jordan_chain_matrix(T, dT, tol: float = 1e-8) -> dict:
    """Jordan-chain test for a matrix-valued spectral problem ``T(omega)``.

    Given ``T = T(omega_*)`` and ``dT = T'(omega_*)``, solve

        T u0 = 0,                 (right null vector)
        T u1 = -dT u0             (generalized vector, solvability required)

    An EP2 requires the *second* equation to be solvable, i.e. ``dT u0`` must be
    orthogonal to the left null vector ``v``: ``v^H dT u0 = 0`` is precisely the
    condition ``dD/domega = 0`` on the determinant.  The returned
    ``solvability`` is ``|v^H dT u0|`` normalized; small means defective.
    """
    T = np.asarray(T, dtype=complex)
    dT = np.asarray(dT, dtype=complex)
    U, S, Vh = np.linalg.svd(T)
    u0 = Vh[-1].conj()          # right null vector
    v = U[:, -1]                # left null vector
    smin = float(S[-1])
    gap = float(S[-2]) if len(S) > 1 else np.inf

    num = complex(v.conj() @ (dT @ u0))
    scale = float(np.linalg.norm(dT @ u0)) or 1.0
    solvability = abs(num) / scale

    rhs = -(dT @ u0)
    u1, *_ = np.linalg.lstsq(T, rhs, rcond=None)
    chain_residual = float(np.linalg.norm(T @ u1 - rhs) / (np.linalg.norm(rhs) or 1.0))

    return {
        "sigma_min": smin,
        "sigma_gap": gap,
        "geometric_multiplicity_1": bool(smin < tol * max(gap, 1.0)),
        "solvability": float(solvability),
        "defective": bool(solvability < 1e-3),
        "chain_residual": chain_residual,
        "u0": u0,
        "u1": u1,
        "v": v,
    }

"""Scale-free exceptional-point diagnostics via Cauchy integrals.

Why not ``|dF/domega|``
-----------------------
The spectral condition ``F`` is defined only up to an arbitrary nonvanishing
analytic prefactor: Solver A's continued fraction can be rescaled, inverted at a
different depth, or multiplied by any ``g(omega) != 0`` without changing its
zero set.  Therefore ``|dF/domega|`` at a root carries **no invariant meaning** --
it can be made as small as one likes by rescaling ``F``.  A bound built from it
is not a bound on the physics.

The invariant object is the *root separation*.  Near a simple root ``w1`` of
``F = c (w - w1)(w - w2) g(w)``, the Taylor coefficients about ``w1`` obey

    a1 = F'(w1)      = c (w1 - w2) g(w1)
    a2 = F''(w1)/2   = c g(w1) + O(w1 - w2)

so

    a1 / a2  ->  (w1 - w2)                                        (*)

and both ``c`` and ``g`` cancel.  ``|a1/a2|`` is a genuine estimate of the
distance to the partner root, invariant under rescaling of ``F``, and it
vanishes exactly at an EP2.  This module uses (*) as the primary diagnostic.

The Taylor coefficients themselves come from a Cauchy integral on a circle
rather than from finite differences.  ``F`` is analytic in ``omega`` away from
its poles, so the FFT of ``F`` sampled on a circle gives every coefficient at
once, with error controlled by the circle radius rather than by subtractive
cancellation.  The same samples give the argument-principle zero count for
free, which is the one EP test that does not depend on any root finder.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "TaylorData",
    "cauchy_taylor",
    "root_separation",
    "zero_count",
    "winding_zero_count",
]


@dataclass(frozen=True)
class TaylorData:
    """Taylor coefficients of ``F`` about ``center`` on a circle of ``radius``."""

    center: complex
    radius: float
    coeffs: np.ndarray  # a[k] = F^(k)(center)/k!
    n_samples: int

    @property
    def a0(self) -> complex:
        return complex(self.coeffs[0])

    @property
    def a1(self) -> complex:
        return complex(self.coeffs[1])

    @property
    def a2(self) -> complex:
        return complex(self.coeffs[2])


def cauchy_taylor(f, center: complex, radius: float, n: int = 64,
                  n_coeffs: int = 8) -> TaylorData:
    """Taylor coefficients of ``f`` about ``center`` from samples on a circle.

    Uses ``a[k] = (1/n) sum_j f(center + radius e^{i t_j}) e^{-i k t_j} / radius^k``,
    the trapezoidal (=spectrally accurate) discretization of the Cauchy integral.

    ``radius`` must be small enough that ``f`` is analytic on and inside the
    circle, but large enough that ``radius^k`` does not underflow the sampling
    error.  ``n_coeffs << n`` keeps the high-``k`` amplification bounded.
    """
    if radius <= 0:
        raise ValueError("radius must be positive")
    if n_coeffs >= n // 2:
        raise ValueError("n_coeffs must be well below n/2 to stay accurate")
    t = 2.0 * np.pi * np.arange(n) / n
    pts = center + radius * np.exp(1j * t)
    vals = np.array([complex(f(p)) for p in pts], dtype=complex)
    if not np.all(np.isfinite(vals)):
        raise FloatingPointError("f returned a non-finite value on the contour")
    spec = np.fft.fft(vals) / n
    ks = np.arange(n_coeffs)
    coeffs = spec[:n_coeffs] / (radius ** ks)
    return TaylorData(center=center, radius=radius, coeffs=coeffs, n_samples=n)


def root_separation(f, root: complex, radius: float, n: int = 64) -> dict:
    """Scale-free estimate of the distance from ``root`` to its nearest partner.

    Returns a dict with ``separation`` = ``|a1/a2|``, the residual ``|a0|``
    (should be ~0 at a converged root), and the raw coefficients.

    ``separation -> 0`` is the signature of an EP2.  A *large* separation is a
    certificate that no EP2 sits at this root, and unlike ``|F'|`` it cannot be
    faked by rescaling ``F``.
    """
    td = cauchy_taylor(f, root, radius, n=n, n_coeffs=6)
    a0, a1, a2 = td.a0, td.a1, td.a2
    sep = abs(a1) / abs(a2) if a2 != 0 else np.inf
    scale = max(abs(a1), abs(a2) * radius)
    return {
        "separation": float(sep),
        "residual": float(abs(a0)),
        "residual_relative": float(abs(a0) / scale) if scale > 0 else np.inf,
        "a0": complex(a0),
        "a1": complex(a1),
        "a2": complex(a2),
        "radius": float(radius),
    }


def winding_zero_count(f, center: complex, radius: float, n: int = 256) -> float:
    """Number of zeros minus poles of ``f`` inside the circle, by winding number.

    Computed as the total change in ``arg f`` around the contour divided by
    ``2 pi``.  This needs no derivative of ``f`` and no root finder, so it is
    genuinely independent evidence.  The return value is a float; it should be
    within a few times ``1e-2`` of an integer for a well-resolved contour.
    """
    t = 2.0 * np.pi * np.arange(n + 1) / n
    pts = center + radius * np.exp(1j * t)
    vals = np.array([complex(f(p)) for p in pts], dtype=complex)
    if np.any(vals == 0):
        raise FloatingPointError("f vanishes on the contour; move or resize it")
    ang = np.unwrap(np.angle(vals))
    return float((ang[-1] - ang[0]) / (2.0 * np.pi))


def zero_count(f, center: complex, radius: float, n: int = 256,
               tol: float = 0.08) -> int:
    """Integer zero count inside the contour, or raise if not resolved."""
    w = winding_zero_count(f, center, radius, n=n)
    k = round(w)
    if abs(w - k) > tol:
        raise ValueError(
            f"winding number {w:.4f} is not close to an integer; "
            "refine the contour (increase n or move the circle)"
        )
    return int(k)

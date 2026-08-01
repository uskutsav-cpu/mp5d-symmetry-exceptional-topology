"""Five-dimensional Myers-Perry geometry with two rotation parameters.

Conventions are fixed and documented in ``docs/CONVENTIONS.md``.  Every formula
here is exercised by ``tests/unit/test_geometry.py`` against closed-form limits.

Metric (Boyer-Lindquist-like coordinates ``(t, r, theta, phi, psi)``,
``theta in [0, pi/2]``, ``phi, psi in [0, 2 pi)``), signature ``(-++++)``::

    ds^2 = -dt^2 + (M / Sigma) (dt - a sin^2 t dphi - b cos^2 t dpsi)^2
           + (Sigma / Delta) dr^2 + Sigma dtheta^2
           + (r^2 + a^2) sin^2 t dphi^2 + (r^2 + b^2) cos^2 t dpsi^2

with

    Sigma = r^2 + a^2 cos^2 t + b^2 sin^2 t
    Pi    = (r^2 + a^2)(r^2 + b^2)
    Delta = (Pi - M r^2) / r^2

``M`` is the Myers-Perry mass parameter (dimension length^2), *not* the scalar
mass.  The scalar mass is written ``mu`` throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = ["MPGeometry", "sdelta_to_ab", "ab_to_sdelta"]


def sdelta_to_ab(s: float, delta: float) -> tuple[float, float]:
    """Adapted rotation coordinates ``(s, delta)`` -> ``(a, b)``."""
    return s + delta, s - delta


def ab_to_sdelta(a: float, b: float) -> tuple[float, float]:
    """``(a, b)`` -> adapted rotation coordinates ``(s, delta)``."""
    return 0.5 * (a + b), 0.5 * (a - b)


@dataclass(frozen=True)
class MPGeometry:
    """Asymptotically flat 5D Myers-Perry black hole.

    Parameters
    ----------
    M : mass parameter (default 1, which fixes the length unit).
    a, b : rotation parameters.
    """

    a: float
    b: float
    M: float = 1.0

    # -- horizon structure ------------------------------------------------
    @property
    def discriminant(self) -> float:
        """``(M - a^2 - b^2)^2 - 4 a^2 b^2``; non-negative iff a horizon exists."""
        return (self.M - self.a**2 - self.b**2) ** 2 - 4 * self.a**2 * self.b**2

    @property
    def has_horizon(self) -> bool:
        """True on the physical (non-naked-singularity) region."""
        return self.discriminant >= 0.0 and (self.M - self.a**2 - self.b**2) >= 0.0

    @property
    def extremality(self) -> float:
        """``M - (|a| + |b|)^2``.  Zero at extremality, positive for sub-extremal."""
        return self.M - (abs(self.a) + abs(self.b)) ** 2

    @property
    def _sqrt_disc(self) -> float:
        # At exactly extremal parameters the discriminant is analytically zero but
        # can evaluate to a tiny negative float; clamping keeps the extremal limit
        # reachable instead of raising.
        return math.sqrt(max(self.discriminant, 0.0))

    @property
    def z_plus(self) -> float:
        """Outer horizon radius squared ``r_+^2``."""
        return 0.5 * (self.M - self.a**2 - self.b**2 + self._sqrt_disc)

    @property
    def z_minus(self) -> float:
        """Inner horizon radius squared ``r_-^2``.  Note ``z_+ z_- = a^2 b^2``."""
        return 0.5 * (self.M - self.a**2 - self.b**2 - self._sqrt_disc)

    @property
    def r_plus(self) -> float:
        return math.sqrt(self.z_plus)

    @property
    def r_minus(self) -> float:
        return math.sqrt(max(self.z_minus, 0.0))

    # -- horizon thermodynamics -------------------------------------------
    @property
    def Omega_a(self) -> float:
        """Horizon angular velocity conjugate to ``phi`` (partner of ``m1``)."""
        return self.a / (self.z_plus + self.a**2)

    @property
    def Omega_b(self) -> float:
        """Horizon angular velocity conjugate to ``psi`` (partner of ``m2``)."""
        return self.b / (self.z_plus + self.b**2)

    @property
    def kappa(self) -> float:
        """Surface gravity ``kappa = r_+ (r_+^2 - r_-^2) / [(r_+^2+a^2)(r_+^2+b^2)]``.

        Reduces to ``1 / (2 r_+)`` for Schwarzschild-Tangherlini, i.e.
        ``T_H = (D-3) / (4 pi r_+) = 1 / (2 pi r_+)`` in ``D = 5``.
        """
        zp = self.z_plus
        return self.r_plus * (zp - self.z_minus) / ((zp + self.a**2) * (zp + self.b**2))

    @property
    def T_H(self) -> float:
        return self.kappa / (2.0 * math.pi)

    # -- adapted coordinates ----------------------------------------------
    @property
    def s(self) -> float:
        """Average rotation ``(a + b) / 2``."""
        return 0.5 * (self.a + self.b)

    @property
    def delta(self) -> float:
        """Spin asymmetry ``(a - b) / 2``."""
        return 0.5 * (self.a - self.b)

    @classmethod
    def from_sdelta(cls, s: float, delta: float, M: float = 1.0) -> "MPGeometry":
        a, b = sdelta_to_ab(s, delta)
        return cls(a=a, b=b, M=M)

    def exchanged(self) -> "MPGeometry":
        """The ``a <-> b`` image, i.e. ``delta -> -delta`` at fixed ``s``."""
        return MPGeometry(a=self.b, b=self.a, M=self.M)

    # -- radial structure --------------------------------------------------
    def P(self, z: complex) -> complex:
        """``P(z) = (z + a^2)(z + b^2) - M z = (z - z_+)(z - z_-)`` with ``z = r^2``."""
        return (z + self.a**2) * (z + self.b**2) - self.M * z

    def Delta(self, r: complex) -> complex:
        z = r * r
        return self.P(z) / z

"""Direct algebraic radial coefficients: no FFT, samples, or degree guessing.

Derived from the radial equation in the audited repository at 72abf0bf.
Ascending coefficients have rigorous structural degree bounds (9, 8, 7).
The known common factor u**3*(1-u)**4 is removed analytically. Polynomial
long division is by a MONIC exact factor, so no numerical degree/pivot test is
used. The returned quotient, not a rounded-away remainder, defines the formula.

The companion symbolic audit proves the remainder is
    p**2 * (W(p)**2 - 4*p**2*(p**2-m**2)**2*sigma**2) * (1-u)**2,
which vanishes by the horizon-exponent identity. This is an exact finite
coefficient identity, NOT a convergence proof for the infinite recurrence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import mpmath as mp

DEGREE_BOUNDS = (9, 8, 7)
COEFFICIENT_VERSION = "mp5d-direct-algebra/1.0"


def add(*polynomials):
    n = max(map(len, polynomials))
    out = [0]*n
    for p in polynomials:
        for j, a in enumerate(p):
            out[j] += a
    return out


def scale(p, c):
    return [c*a for a in p]


def mul(*polynomials):
    out = [1]
    for p in polynomials:
        q = [0]*(len(out)+len(p)-1)
        for i, a in enumerate(out):
            for j, b in enumerate(p):
                q[i+j] += a*b
        out = q
    return out


def power(p, n):
    if not isinstance(n, int) or n < 0:
        raise ValueError("polynomial exponent must be a nonnegative integer")
    out = [1]
    for _ in range(n):
        out = mul(out, p)
    return out


def derivative(p):
    return [j*p[j] for j in range(1, len(p))] or [0]


def evaluate(p, x):
    value = 0
    for coefficient in reversed(p):
        value = value*x + coefficient
    return value


def divide_monic(p, divisor):
    """Long division in a scalar ring; no trimming or near-zero comparisons."""
    if len(divisor) < 2 or divisor[-1] != 1:
        raise ValueError("exact monic divisor required")
    rem = list(p)
    degree = len(divisor)-1
    if len(rem) <= degree:
        return [0], rem
    quotient = [0]*(len(rem)-degree)
    for i in range(len(quotient)-1, -1, -1):
        q = rem[i+degree]
        quotient[i] = q
        for j, coefficient in enumerate(divisor):
            rem[i+j] -= q*coefficient
    return quotient, rem[:degree]


@dataclass(frozen=True)
class RadialParameters:
    a: str
    b: str
    mu: str
    m1: int
    m2: int
    ell: int
    M: str = "1"

    def __post_init__(self):
        if not all(isinstance(v, str) for v in (self.a, self.b, self.mu, self.M)):
            raise TypeError("use decimal strings to preserve inputs under precision ladders")
        if any(not isinstance(n, int) for n in (self.m1, self.m2, self.ell)):
            raise TypeError("mode labels must be integers")
        if self.ell < abs(self.m1)+abs(self.m2) or (self.ell-abs(self.m1)-abs(self.m2)) % 2:
            raise ValueError("ell is incompatible with azimuthal labels")
        a, b, mu, M = map(mp.mpf, (self.a, self.b, self.mu, self.M))
        if not all(mp.isfinite(v) for v in (a, b, mu, M)) or mu < 0 or M <= 0:
            raise ValueError("invalid physical parameters")
        if (abs(a)+abs(b))**2 >= M:
            raise ValueError("strictly subextremal geometry required")

    @property
    def angular_index(self):
        return (self.ell-abs(self.m1)-abs(self.m2))//2


def horizons(a, b, M, sqrt):
    h = M-a*a-b*b
    p = sqrt((h+sqrt(h*h-4*a*a*b*b))/2)
    # Stable at a=0 or b=0; avoids cancellation in h-sqrt(discriminant).
    m = sqrt(a*a*b*b)/p
    return p, m


def coefficient_formula(a: Any, b: Any, mu: Any, m1: int, m2: int,
                        omega: Any, lam: Any, p: Any, m: Any,
                        outgoing_k: Any, imaginary_unit: Any = 1j):
    """Scalar-generic algebra: usable with mpmath, exact SymPy, or Arb/Acb.

    Inputs must satisfy outgoing_k**2=omega**2-mu**2 and the physical horizon
    identities. The caller chooses a consistent analytic square-root sheet.
    No floating-point cast occurs anywhere in this coefficient construction.
    """
    u, t, R = [0, 1], [1, -1], [p, -m]
    d = p-m
    q, v = -(p*p+m*m), p*p*m*m
    L = mul([2*p, -(p+m)], [p+m, -2*m])
    D = scale(mul(u, L), d*d)
    W = add(scale(mul(add(power(R, 2), scale(power(t, 2), a*a)),
                          add(power(R, 2), scale(power(t, 2), b*b))), omega),
            scale(mul(add(power(R, 2), scale(power(t, 2), b*b)), power(t, 2)), -m1*a),
            scale(mul(add(power(R, 2), scale(power(t, 2), a*a)), power(t, 2)), -m2*b))
    Wp = (p*p+a*a)*(p*p+b*b)*omega-m1*a*(p*p+b*b)-m2*b*(p*p+a*a)
    sigma = Wp/(2*p*(p*p-m*m))
    G = a*b*omega-a*m2-b*m1
    H = add(scale(power(t, 2), -imaginary_unit*sigma),
            scale(mul(u, t), -(p-p+3)/2), scale(u, imaginary_unit*outgoing_k*d))
    J = add(power(H, 2), mul(u, power(t, 2), derivative(H)),
            scale(mul(H, add(power(t, 2), scale(mul(u, t), -2))), -1))
    Bn = add(scale(mul(D, t, R), -2), scale(mul(power(R, 4), t), 3*d),
             scale(mul(power(t, 3), power(R, 2)), d*q), scale(power(t, 5), -d*v))
    Cc = -(a*a+b*b)*omega*omega+2*omega*(a*m1+b*m2)-lam
    A = scale(mul(u, power(L, 2), power(t, 2), power(R, 2)), d*d)
    B = add(scale(mul(power(L, 2), power(R, 2), H), 2*d*d), mul(L, R, Bn))
    numerator = add(scale(mul(power(L, 2), power(R, 2), J), d*d), mul(L, R, Bn, H),
                    mul(power(R, 2), power(W, 2)), scale(mul(D, power(R, 6)), -mu*mu),
                    scale(mul(power(t, 4), D, power(R, 2)), -G*G),
                    scale(mul(power(t, 2), D, power(R, 4)), Cc))
    C, remainder = divide_monic(numerator, [0, 1, -2, 1])
    if tuple(len(v)-1 for v in (A, B, C)) != DEGREE_BOUNDS:
        raise AssertionError("internal structural degree changed")
    return A, B, C, remainder


def coefficients(params: RadialParameters, omega, lam, *, dps: int = 50,
                 outgoing_k=None):
    """Arbitrary-precision coefficient evaluation; precision is actually used."""
    if dps < 15:
        raise ValueError("at least 15 decimal digits required")
    with mp.workdps(dps):
        a, b, mu, M = map(mp.mpf, (params.a, params.b, params.mu, params.M))
        w, L = mp.mpc(omega), mp.mpc(lam)
        p, m = horizons(a, b, M, mp.sqrt)
        k = mp.sqrt(w*w-mu*mu) if outgoing_k is None else mp.mpc(outgoing_k)
        if k.real < 0:
            k = -k
        if abs(k*k-(w*w-mu*mu)) > mp.power(10, -dps+7)*max(1, abs(w*w)):
            raise ValueError("outgoing momentum does not satisfy the dispersion relation")
        return coefficient_formula(a, b, mu, params.m1, params.m2, w, L, p, m, k, mp.j)


def ball_coefficients(params: RadialParameters, omega_real: str, omega_imag: str,
                      lambda_real: str, lambda_imag: str, *, bits: int = 192):
    """Outward-rounded coefficient balls for the finite coefficient formula.

    Requires python-flint. Lambda must be enclosed separately; this function
    never claims to enclose an unknown angular eigenvalue by rounding a guess.
    Returned balls do not certify a root of the infinite QNM problem.
    """
    try:
        from flint import acb, arb, ctx
    except ImportError as exc:
        raise RuntimeError("python-flint is required; no floating fallback is allowed") from exc
    old = ctx.prec
    try:
        ctx.prec = bits
        a, b, mu, M = (arb(x) for x in (params.a, params.b, params.mu, params.M))
        p, m = horizons(a, b, M, lambda x: x.sqrt())
        omega = acb(arb(omega_real), arb(omega_imag))
        lam = acb(arb(lambda_real), arb(lambda_imag))
        k = (omega*omega-mu*mu).sqrt()
        if not k.real > 0:
            raise ValueError("outgoing square-root sheet is not separated from the cut")
        A, B, C, _ = coefficient_formula(a, b, mu, params.m1, params.m2, omega, lam,
                                         p, m, k, acb(0, 1))
        return {"A": A, "B": B, "C": C, "bits": bits,
                "scope": "FINITE_COEFFICIENT_FORMULA_ONLY", "degree_bounds": DEGREE_BOUNDS}
    finally:
        ctx.prec = old

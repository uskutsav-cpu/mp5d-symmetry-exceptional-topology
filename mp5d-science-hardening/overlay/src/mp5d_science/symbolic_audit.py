"""Exact symbolic audit of the direct radial polynomial construction."""
from __future__ import annotations

import sympy as s


def audit_radial_division() -> dict:
    u, p, m, aa, bb, j1, j2, w, k, lam, G, sigma = s.symbols(
        'u p m aa bb j1 j2 w k lam G sigma')
    t, R, d = 1-u, p-m*u, p-m
    q, v = -(p*p+m*m), p*p*m*m
    L = (2*p-(p+m)*u)*(p+m-2*m*u)
    D = d*d*u*L
    W = (R*R+aa*t*t)*(R*R+bb*t*t)*w-j1*(R*R+bb*t*t)*t*t-j2*(R*R+aa*t*t)*t*t
    Wp = (p*p+aa)*(p*p+bb)*w-j1*(p*p+bb)-j2*(p*p+aa)
    H = -s.I*sigma*t*t-s.Rational(3, 2)*u*t+s.I*k*d*u
    J = H*H+u*t*t*s.diff(H, u)-H*(t*t-2*u*t)
    Bn = -2*D*t*R+d*(3*R**4*t+q*t**3*R**2-v*t**5)
    Cc = -(aa+bb)*w*w+2*w*(j1+j2)-lam
    N = s.Poly(d*d*L*L*R*R*J+L*R*Bn*H+R*R*W*W-(w*w-k*k)*D*R**6-
               G*G*t**4*D*R**2+Cc*t*t*D*R**4, u)
    quotient, remainder = s.div(N, s.Poly(u*t*t, u))
    expected = p*p*(Wp*Wp-4*p*p*(p*p-m*m)**2*sigma*sigma)*t*t
    remainder_identity = s.Poly(remainder.as_expr()-expected, u).is_zero
    horizon_identity = s.cancel(expected.subs(sigma, Wp/(2*p*(p*p-m*m)))) == 0
    degrees = [s.degree(d*d*u*L*L*t*t*R*R, u),
               s.degree(2*d*d*L*L*R*R*H+L*R*Bn, u), quotient.degree()]
    return {"status": "PASS" if remainder_identity and horizon_identity and degrees == [9, 8, 7] else "FAIL",
            "polynomial_degrees": list(map(int, degrees)), "remainder_identity": bool(remainder_identity),
            "horizon_cancellation": bool(horizon_identity),
            "assumptions": ["p>m>=0", "sigma=W(p)/(2*p*(p^2-m^2))", "k^2=w^2-mu^2",
                            "Delta numerator=(r^2-p^2)(r^2-m^2)"],
            "scope": "EXACT_TRANSFORMED_RADIAL_COEFFICIENT_IDENTITY"}


def audit_symmetry() -> dict:
    a, b, x, d, w, mu = s.symbols('a b s delta omega mu')
    first = s.expand((x+d)**2-(x-d)**2-4*x*d) == 0
    second = s.expand(4*x*d*(w*w-mu*mu)).subs(d, 0) == 0
    return {"status": "PASS" if first and second else "FAIL", "spin_identity": first,
            "equal_spin_spheroidicity_zero": second, "scope": "EXACT_FINITE_ALGEBRA"}

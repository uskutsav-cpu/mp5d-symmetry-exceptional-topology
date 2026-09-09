"""Rigorous, exact-rational interval Krawczyk checks for FINITE polynomials.

No binary floating-point operation participates in the proof. This fallback is
slower than Arb but makes small finite algebraic certificates reproducible even
without python-flint. It does not turn a finite radial determinant into a
continuum quasinormal-mode theorem.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from typing import Sequence


def rational(x):
    if isinstance(x, float):
        raise TypeError("binary floats are not accepted as exact decimal inputs; use strings")
    return x if isinstance(x, Q) else Q(x)


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self):
        object.__setattr__(self, 'lo', rational(self.lo))
        object.__setattr__(self, 'hi', rational(self.hi))
        if self.lo > self.hi:
            raise ValueError("empty interval")

    @classmethod
    def point(cls, x):
        return cls(rational(x), rational(x))

    @classmethod
    def around(cls, midpoint, radius):
        midpoint, radius = rational(midpoint), rational(radius)
        if radius <= 0:
            raise ValueError("positive exact radius required")
        return cls(midpoint-radius, midpoint+radius)

    def __add__(self, other):
        if not isinstance(other, Interval):other=Interval.point(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__=__add__

    def __neg__(self):
        return Interval(-self.hi,-self.lo)

    def __sub__(self,other):
        if not isinstance(other,Interval):other=Interval.point(other)
        return self+(-other)

    def __rsub__(self,other):
        return Interval.point(other)-self

    def __mul__(self,other):
        if not isinstance(other,Interval):other=Interval.point(other)
        products=[self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi]
        return Interval(min(products),max(products))

    __rmul__=__mul__

    def reciprocal(self):
        if self.lo<=0<=self.hi:
            raise ZeroDivisionError("interval contains zero")
        return Interval(1/self.hi,1/self.lo)

    def __truediv__(self,other):
        if not isinstance(other,Interval):other=Interval.point(other)
        return self*other.reciprocal()

    def abs_upper(self):
        return max(abs(self.lo),abs(self.hi))

    def inside(self,other):
        return other.lo<self.lo and self.hi<other.hi

    def to_dict(self):
        return {'lo':str(self.lo),'hi':str(self.hi)}


def cadd(a,b):
    return a[0]+b[0],a[1]+b[1]


def cmul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def exact_complex(x):
    if not isinstance(x,(list,tuple)) or len(x)!=2:
        raise TypeError("complex coefficients must be [exact real, exact imaginary]")
    return rational(x[0]),rational(x[1])


def polynomial_value(coefficients,z):
    out=(0,0)
    for coefficient in coefficients:
        out=cadd(cmul(out,z),coefficient)
    return out


def polynomial_derivative(coefficients):
    degree=len(coefficients)-1
    return [(a*(degree-j),b*(degree-j)) for j,(a,b) in enumerate(coefficients[:-1])]


def krawczyk_polynomial(coefficients: Sequence[Sequence], center: Sequence, radius,
                       *, problem_id: str) -> dict:
    """Existence and uniqueness inside a rectangular complex box, when PASS.

    Coefficients are descending exact rational complex pairs. The polynomial
    coefficients themselves must already represent the exact finite problem;
    this routine cannot certify uncertainty in rounded input coefficients.
    """
    if not problem_id or len(coefficients)<2:
        raise ValueError("a named nonconstant finite polynomial is required")
    coefficients=[exact_complex(c) for c in coefficients]
    if coefficients[0]==(0,0):
        raise ValueError("leading coefficient is zero")
    center=exact_complex(center)
    r=rational(radius)
    box=[Interval.around(x,r) for x in center]
    f=polynomial_value(coefficients,center)
    derivative=polynomial_derivative(coefficients)
    a,b=polynomial_value(derivative,center)
    denominator=a*a+b*b
    if denominator==0:
        return {'status':'INCONCLUSIVE','reason':'singular point preconditioner',
                'scope':'FINITE_EXACT_POLYNOMIAL','problem_id':problem_id,
                'continuum_certificate':False}
    R=[[a/denominator,b/denominator],[-b/denominator,a/denominator]]
    da,db=polynomial_value(derivative,(box[0],box[1]))
    J=[[da,-db],[db,da]]
    E=[[Interval.point(int(i==j))-sum((R[i][k]*J[k][j] for k in range(2)),Interval.point(0))
        for j in range(2)] for i in range(2)]
    corrected=[center[i]-sum(R[i][j]*f[j] for j in range(2)) for i in range(2)]
    delta=[box[i]-center[i] for i in range(2)]
    image=[Interval.point(corrected[i])+sum((E[i][j]*delta[j] for j in range(2)),Interval.point(0))
           for i in range(2)]
    contraction=max(sum(E[i][j].abs_upper() for j in range(2)) for i in range(2))
    inclusion=all(image[i].inside(box[i]) for i in range(2))
    passed=inclusion and contraction<1
    return {'status':'CERTIFIED_UNIQUE_FINITE_ROOT' if passed else 'INCONCLUSIVE',
            'problem_id':problem_id,'scope':'FINITE_EXACT_POLYNOMIAL',
            'arithmetic':'EXACT_RATIONAL_INTERVALS','polynomial_degree':len(coefficients)-1,
            'box':[x.to_dict() for x in box],'krawczyk_image':[x.to_dict() for x in image],
            'strict_inclusion':inclusion,'contraction_upper':str(contraction),
            'continuum_certificate':False}


def radial_finite_polynomial(ell: int, size: int) -> dict:
    """Exact finite Schwarzschild-Tangherlini radial Hill polynomial, M=1, mu=0.

    Its determinant is computed from the RAW upper-Hessenberg recurrence by a
    division-free determinant recurrence. No Gaussian-reduced CF is used.
    This narrowly scoped certificate must not be advertised as generic spinning
    radial certification, nor as a certificate for the infinite boundary value
    problem. Horizon p=1, inner horizon m=0, k=omega, Lambda=ell(ell+2).
    """
    import sympy as s
    from .radial_polynomial import coefficient_formula
    if ell<0 or not 2<=size<=64:
        raise ValueError("ell>=0 and finite size 2..64 required")
    w=s.Symbol('w')
    A,B,C,_=coefficient_formula(s.Integer(0),s.Integer(0),s.Integer(0),0,0,w,
                                s.Integer(ell*(ell+2)),s.Integer(1),s.Integer(0),w,s.I)
    def get(p,i):return p[i] if 0<=i<len(p) else s.Integer(0)
    def entry(n,j):
        return s.Poly(s.expand(get(A,n-j+2)*j*(j-1)+get(B,n-j+1)*j+get(C,n-j)),w,domain=s.QQ_I)
    one=s.Poly(1,w,domain=s.QQ_I)
    det=[one]
    for n in range(1,size+1):
        total=s.Poly(0,w,domain=s.QQ_I)
        product=one
        for length in range(1,min(n,8)+1):
            j=n-length
            if length>1:product*=entry(j,j+1)
            total+=(-1)**(length-1)*entry(n-1,j)*product*det[j]
        det.append(total)
    poly=det[-1]
    # Divide by a nonzero constant to reduce exact rational coefficient sizes.
    poly=poly.monic()
    coefficients=[]
    for coefficient in poly.all_coeffs():
        re,im=s.expand_complex(coefficient).as_real_imag()
        coefficients.append([str(re),str(im)])
    return {'problem_id':f'static-massless-radial-l{ell}-K{size}',
            'coefficients':coefficients,'size':size,'ell':ell,'degree':poly.degree(),
            'scope':'FINITE_RAW_RADIAL_DETERMINANT','continuum_certificate':False}

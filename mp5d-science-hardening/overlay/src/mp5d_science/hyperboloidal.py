"""Hyperboloidal quadratic-pencil TESTBED, not a massive MP5D Solver F.

For V(x)=V0*sech(x)^2, y=tanh(x), tau=t-log(cosh(x)), the ansatz
psi=exp(-i*w*tau)*u(y) gives
 (1-y*y)u'' + 2*(i*w-1)*y*u' + (w*w+i*w-V0)u = 0.
The degenerate endpoint equations impose regularity. This testbed has exact
frequencies +/-sqrt(V0-1/4)-i*(n+1/2). Its success does NOT establish the
Gevrey/radiation selection required at asymptotically flat massive infinity.
"""
from __future__ import annotations
import numpy as np
from scipy import linalg


def chebyshev_lobatto(n: int):
    if not isinstance(n,int) or n<6:
        raise ValueError('need at least six polynomial intervals')
    y=np.cos(np.pi*np.arange(n+1)/n)
    weights=(-1.)**np.arange(n+1)
    weights[[0,-1]]*=2
    delta=y[:,None]-y[None,:]
    D=(weights[:,None]/weights[None,:])/(delta+np.eye(n+1))
    D-=np.diag(np.sum(D,axis=1))
    return y,D


def quadratic_spectrum(P0,P1,P2):
    matrices=[np.asarray(p,dtype=complex) for p in (P0,P1,P2)]
    n=matrices[0].shape[0]
    if any(p.shape!=(n,n) or not np.all(np.isfinite(p)) for p in matrices):
        raise ValueError('finite square coefficient matrices of the same size required')
    P0,P1,P2=matrices
    A=np.block([[np.zeros((n,n)),np.eye(n)],[-P0,-P1]])
    B=linalg.block_diag(np.eye(n),P2)
    values,vectors=linalg.eig(A,B,check_finite=False)
    records=[]
    for w,vector in zip(values,vectors.T,strict=True):
        if not np.isfinite(w):continue
        u=vector[:n]
        terms=[P0@u,w*(P1@u),w*w*(P2@u)]
        scale=sum(np.linalg.norm(t) for t in terms)
        residual=np.linalg.norm(sum(terms))/max(scale,np.finfo(float).tiny)
        records.append({'omega':[float(w.real),float(w.imag)],'residual':float(residual)})
    return records


def poschl_teller_spectrum(V0: float=1.,n: int=32):
    if not np.isfinite(V0) or V0<=.25:
        raise ValueError('this testbed requires V0>1/4')
    y,D=chebyshev_lobatto(n)
    I=np.eye(n+1)
    P0=np.diag(1-y*y)@(D@D)-2*np.diag(y)@D-V0*I
    P1=1j*(2*np.diag(y)@D+I)
    return {'status':'FINITE_TESTBED_SPECTRUM','model':'Poschl-Teller',
        'V0':V0,'resolution':n,'records':quadratic_spectrum(P0,P1,I),
        'mp5d_solver_F_validated':False,'massive_quasiresonant_coverage':False}


def testbed_benchmark(V0=1.,resolutions=(20,28,36),max_overtone=3):
    from scipy.optimize import linear_sum_assignment
    rows=[]
    for n in resolutions:
        result=poschl_teller_spectrum(V0,n)
        values=np.array([complex(*r['omega']) for r in result['records']])
        targets=np.array([np.sqrt(V0-.25)-1j*(k+.5) for k in range(max_overtone+1)])
        row,col=linear_sum_assignment(abs(targets[:,None]-values[None,:]))
        errors=abs(targets[row]-values[col])
        rows.append({'resolution':n,'max_error':float(np.max(errors)),
                     'matched':[result['records'][i] for i in col]})
    return {'status':'TESTBED_PASS' if all(r['max_error']<1e-5 for r in rows) else 'INCONCLUSIVE',
            'scope':'EXACTLY_SOLVABLE_TESTBED_ONLY','records':rows,'mp5d_solver_F_validated':False}

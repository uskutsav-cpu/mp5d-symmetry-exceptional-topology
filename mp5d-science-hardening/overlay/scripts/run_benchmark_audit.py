#!/usr/bin/env python3
"""Published static benchmarks, finite certification and a hyperboloidal testbed."""
from __future__ import annotations
import argparse
import cmath
import math
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import mpmath as mp
from mp5d_science.physical import solve_a,solve_c
from mp5d_science.radial_polynomial import RadialParameters
from mp5d_science.provenance import envelope,read_json,atomic_json
from mp5d_science.symbolic_audit import audit_radial_division
from mp5d_science.hyperboloidal import testbed_benchmark
from mp5d_science.certification import radial_finite_polynomial,krawczyk_polynomial


def finite_certificate(size=12):
    exact=radial_finite_polynomial(0,size)
    with mp.workdps(65):
        coeff=[mp.mpc(mp.mpf(a),mp.mpf(b)) for a,b in exact['coefficients']]
        root=mp.findroot(lambda z:mp.polyval(coeff,z),(.53-.38j,.531-.381j),tol=mp.mpf('1e-58'),maxsteps=80)
        center=[mp.nstr(root.real,55),mp.nstr(root.imag,55)]
    certificate=krawczyk_polynomial(exact['coefficients'],center,'1e-30',problem_id=exact['problem_id'])
    return {'exact_problem':exact,'center':center,'certificate':certificate}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--development',action='store_true')
    parser.add_argument('--max-overtone',type=int,default=4)
    args=parser.parse_args()
    if not 0<=args.max_overtone<=4:parser.error('published manifest contains overtones 0..4')
    manifest=read_json(ROOT/'data/regressions/matyjasek_2021.json')
    rows=[]
    for record in manifest['records']:
        if record['overtone']>args.max_overtone:continue
        with mp.workdps(40):
            value=mp.mpc(*record['omega_tilde'])/(2*mp.pi)
            seed=complex(value)
        # Compatible sector for odd/even ell; static spectrum depends only on ell.
        params=RadialParameters('0','0','0',record['ell']%2,0,record['ell'])
        measurements=[]
        for depth in (160,240,320):
            try:
                sol=solve_a(params,seed,overtone=record['overtone'],depth=depth,angular_n=16,dps=45)
                measurements.append(sol.to_dict()|{'absolute_error_from_published':abs(sol.omega-seed)})
            except Exception as exc:measurements.append({'status':'NUMERICAL_FAILURE','error':str(exc),'cf_depth':depth})
        minimum=math.degrees(cmath.phase(complex(seed.real,-seed.imag)))
        angle=min(88.5,max(60.,(minimum+89.)/2))
        rate=(seed*cmath.exp(1j*math.radians(angle))).imag
        length=max(60.,16/rate)
        for n in (160,240,320):
            try:
                sol=solve_c(params,seed,radial_n=n,angular_n=16,length=length,angle_deg=angle)
                measurements.append(sol.to_dict()|{'absolute_error_from_published':abs(sol.omega-seed)})
            except Exception as exc:measurements.append({'status':'NUMERICAL_FAILURE','error':str(exc),'radial_n':n})
        passed=all(m.get('converged') and m.get('absolute_error_from_published',math.inf)<1e-5 for m in measurements)
        rows.append({'id':record['id'],'normalized_published_frequency':[seed.real,seed.imag],
                     'status':'NUMERICAL_BENCHMARK_PASS' if passed else 'CONVERGENCE_NOT_ESTABLISHED',
                     'records':measurements})
        print(record['id'],rows[-1]['status'],flush=True)
    payload={'status':'DIAGNOSTIC_AUDIT','benchmark_source':manifest['source'],'benchmarks':rows,
             'symbolic_radial':audit_radial_division(), 'finite_certificate':finite_certificate(),
             'hyperboloidal_testbed':testbed_benchmark(),
             'scope':'STATIC_BENCHMARKS_FINITE_CERTIFICATION_TESTBED_NOT_SUBMISSION_PASS'}
    ready=(len(rows)==15 and all(r['status']=='NUMERICAL_BENCHMARK_PASS' for r in rows)
           and payload['symbolic_radial']['status']=='PASS'
           and payload['finite_certificate']['certificate']['status']=='CERTIFIED_UNIQUE_FINITE_ROOT'
           and payload['hyperboloidal_testbed']['status']=='TESTBED_PASS')
    payload['status']='SUBMISSION_BENCHMARK_PASS' if ready else 'SCIENTIFIC_REVALIDATION_REQUIRED'
    atomic_json(args.output,envelope('benchmark-and-algebra-audit',payload,ROOT,allow_development=args.development))
    return 0 if ready else 2
if __name__=='__main__':sys.exit(main())

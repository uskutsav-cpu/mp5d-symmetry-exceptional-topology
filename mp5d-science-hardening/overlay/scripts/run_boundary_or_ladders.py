#!/usr/bin/env python3
"""Physical boundary/ladders runner. Invalid anchors cause an explicit BLOCKED result."""
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mp5d_science.research import PairAnchor,physical_boundary,physical_ladders
from mp5d_science.convergence import Resolution
from mp5d_science.provenance import envelope,read_json,atomic_json

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task',choices=['boundary','ladders'])
    parser.add_argument('--anchor-file',type=Path,help='explicit single anchor record; never infer overtone labels')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--development',action='store_true')
    args=parser.parse_args()
    fixture=read_json(ROOT/'data/regressions/ten_interactions.json')
    record=read_json(args.anchor_file) if args.anchor_file else fixture['records'][0]
    plan=read_json(ROOT/'config/science/research_plan.json')
    anchor=PairAnchor(tuple(float(record[k]) for k in ('s','delta','mu')),record['m1'],record['m2'],record['ell'],
        tuple(record['branches']),tuple(record['overtone_labels']),tuple(complex(*z) for z in record['saved_frequencies']))
    base=Resolution(cf_depth=240,angular_n=24,radial_n=220,precision_dps=45,contour_length=40,scaling_angle_deg=65)
    try:
        payload=(physical_boundary(anchor,base,plan['adaptive_boundary']) if args.task=='boundary'
                 else physical_ladders(anchor,base,plan['axes']))
    except Exception as exc:
        payload={'status':'BLOCKED','reason':str(exc),'exception_type':type(exc).__name__,
            'initial_anchor':record,'global_exclusion_proved':False,
            'rule':'Do not optimize or recertify a gap whose branch identity failed preflight.'}
    kind='boundary-refinement' if args.task=='boundary' else 'dangerous-point-convergence'
    atomic_json(args.output,envelope(kind,payload,ROOT,allow_development=args.development))
    print(payload['status'])
    return 0 if payload['status'] in {'PASS','LOCALLY_STABLE_NUMERICAL_MINIMUM'} else 2
if __name__=='__main__':sys.exit(main())

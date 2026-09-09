#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mp5d_science.provenance import read_json
from mp5d_science.claims import check_claim_registry,withdrawn_paths

def main():
    parser=argparse.ArgumentParser(description='Check active claims for withdrawn diagnostics and missing evidence.')
    parser.add_argument('--registry',type=Path,default=ROOT/'config/science/claims.json')
    parser.add_argument('--active-result',type=Path,action='append',default=[])
    args=parser.parse_args();errors=[]
    try:
        errors.extend(check_claim_registry(read_json(args.registry)))
        for path in args.active_result:
            errors.extend(f'{path}: withdrawn diagnostic at {x}' for x in withdrawn_paths(read_json(path)))
    except (OSError,ValueError,TypeError) as exc:errors.append(str(exc))
    for error in errors:print(error,file=sys.stderr)
    return 1 if errors else 0
if __name__=='__main__':sys.exit(main())

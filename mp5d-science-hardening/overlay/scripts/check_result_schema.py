#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from mp5d_science.provenance import read_json,check_result

def main():
    parser=argparse.ArgumentParser(description='Validate envelopes without promoting scientific status.')
    parser.add_argument('files',type=Path,nargs='+')
    args=parser.parse_args();errors=[]
    for path in args.files:
        try:errors.extend(f'{path}: {e}' for e in check_result(read_json(path)))
        except (OSError,ValueError,TypeError) as exc:errors.append(f'{path}: {exc}')
    for error in errors:print(error,file=sys.stderr)
    if not errors:print('Result envelopes are structurally valid; this is not a science PASS.')
    return 1 if errors else 0
if __name__=='__main__':sys.exit(main())

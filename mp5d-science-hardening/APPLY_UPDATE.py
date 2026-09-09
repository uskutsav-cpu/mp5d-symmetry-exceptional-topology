#!/usr/bin/env python3
"""Apply the reviewed overlay to a CLEAN checkout of the exact audited revision.

Never fetch, reset, switch branches, commit, push, delete user changes or change
credentials. Preflight validates every change before writing anything. Existing
files are restored if a write fails. Run --dry-run first and inspect git diff.
"""
from __future__ import annotations
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import os
import subprocess
import stat
import sys
import tempfile

BASE='72abf0bf68d41d07de0fbd5056c40b7a713a6cb5'
REPOSITORY='uskutsav-cpu/mp5d-symmetry-exceptional-topology'
HERE=Path(__file__).resolve().parent


def git(root,*args):
    return subprocess.check_output(['git','-C',str(root),*args],text=True,stderr=subprocess.PIPE).strip()


def replace_once(source,old,new):
    if source.count(old)!=1:raise ValueError('expected exact source fragment was not found exactly once: '+old[:80])
    return source.replace(old,new,1)


def replace_function_body(source,name,body):
    """AST-located body replacement retains the public signature for explicit errors."""
    tree=ast.parse(source)
    nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name==name]
    if len(nodes)!=1:raise ValueError('expected one top-level function '+name)
    node=nodes[0]
    lines=source.splitlines(keepends=True)
    lines[node.body[0].lineno-1:node.end_lineno]=['    '+line+'\n' for line in body.splitlines()]
    updated=''.join(lines)
    ast.parse(updated)
    return updated


def withdrawal_body(name):
    return (f'"""Retired unsafe legacy API: {name}. See docs/SCIENCE_HARDENING.md."""\n'
        f'raise RuntimeError("WITHDRAWN legacy {name}: use mp5d_science explicit evidence APIs; "\n'
        '                   "an unconverged/meromorphic calculation is not EP evidence")')


def overlay_files():
    files={}
    for path in sorted((HERE/'overlay').rglob('*')):
        if not path.is_file():continue
        relative=path.relative_to(HERE/'overlay')
        if '__pycache__' in relative.parts or '.pytest_cache' in relative.parts or path.suffix=='.pyc':continue
        files[relative.as_posix()]=path.read_bytes()
    manifest_path=HERE/'PAYLOAD_SHA256.json'
    if not manifest_path.exists():raise ValueError('bundle integrity manifest missing')
    manifest=json.loads(manifest_path.read_text())
    actual={name:sha256(content).hexdigest() for name,content in files.items()}
    if actual!=manifest['overlay_files']:raise ValueError('bundle payload hash mismatch')
    return files


def build_changes(root):
    if git(root,'rev-parse','HEAD')!=BASE:
        raise ValueError('checkout must be the exact audited research/absolute-final commit '+BASE)
    if git(root,'status','--porcelain','--untracked-files=all'):
        raise ValueError('checkout is not clean; preserve/commit your own changes first')
    remote=git(root,'remote','get-url','origin').removesuffix('.git').rstrip('/')
    if not remote.endswith('/'+REPOSITORY) and not remote.endswith(':'+REPOSITORY):
        raise ValueError('origin is not '+REPOSITORY)
    branch=git(root,'branch','--show-current')
    if not branch or branch in {'main','master','research/absolute-final'}:
        raise ValueError('create a new working branch before applying; do not patch the baseline branch')
    changes=overlay_files()
    collisions=[path for path in changes if (root/path).exists()]
    if collisions:raise ValueError('new overlay paths already exist: '+', '.join(collisions))
    def existing(name):
        current=(root/name).read_text()
        original=subprocess.check_output(['git','-C',str(root),'show',f'{BASE}:{name}']).decode()
        if current!=original:raise ValueError('base file differs: '+name)
        return current
    name='src/mp5d/continuation/branch.py'
    changes[name]=replace_once(existing(name),'extrapolate(hist, ts[:-1] + [arc], arc)',
                              'extrapolate(hist, ts[:-1], arc)').encode()
    name='src/mp5d/exceptional/diagnostics.py'
    content=existing(name)
    for fn in ('root_separation','zero_count'):
        content=replace_function_body(content,fn,withdrawal_body(fn))
    changes[name]=content.encode()
    name='src/mp5d/exceptional/verification.py'
    content=existing(name)
    for fn in ('roots_in_disc','track_roots_on_loop','jordan_chain_matrix'):
        content=replace_function_body(content,fn,withdrawal_body(fn))
    changes[name]=content.encode()
    name='src/mp5d/radial/qnm.py'
    content=existing(name)
    tree=ast.parse(content)
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='solve_qnm_cf')
    lines=content.splitlines(keepends=True)
    fragment=''.join(lines[node.lineno-1:node.end_lineno])
    fragment=replace_once(fragment,'    return QNMSolution(',"""    from mp5d_science.physical import check_existing_cf_root
    from mp5d_science.radial_polynomial import RadialParameters
    _root_guard = check_existing_cf_root(
        RadialParameters(str(a), str(b), str(field_mass), m1, m2, ell, str(M)),
        complex(omega), depth=depth_schedule[-1], angular_n=angular_N, overtone=overtone)
    return QNMSolution(""")
    fragment=replace_once(fragment,'converged=ok,','converged=bool(ok and _root_guard["valid"]),')
    lines[node.lineno-1:node.end_lineno]=[fragment]
    changes[name]=''.join(lines).encode()
    name='src/mp5d/radial/solver_c.py'
    content=existing(name)
    body="""\"\"\"Corrected locally analytic bordered-system solver; iterations reports evaluations.\"\"\"
from mp5d_science.physical import solve_c
from mp5d_science.radial_polynomial import RadialParameters
params = RadialParameters(str(a), str(b), str(mu), m1, m2, ell, str(M))
out = solve_c(params, initial_frequency, radial_n=resolution, angular_n=angular_N,
              length=L, angle_deg=float(np.degrees(theta)), tol=tol,
              max_evaluations=maxiter)
return SolverCResult(
    omega=out.omega, Lambda=complex(float(out.angular_decimal[0]), float(out.angular_decimal[1])),
    smallest_singular_value=out.diagnostics['normalized_sigma_min'],
    singular_value_gap=out.diagnostics['singular_value_gap'], residual=out.residual,
    theta=theta, L=L, resolution=resolution, converged=out.converged,
    iterations=out.diagnostics['function_evaluations'], solver=out.solver)
"""
    changes[name]=replace_function_body(content,'solve_qnm_c',body).encode()
    name='tests/unit/test_ep_machinery_synthetic.py'
    existing(name)
    changes[name]=b'''"""Legacy unsafe API retirement. Positive/negative EP controls now live in
 tests/submission/test_evidence.py, test_spectral.py and test_tracking.py.
They include polynomial EP2/EP3, semisimple crossings, monodromy and Puiseux.
This replaces the former universal Taylor-ratio/rescaling assumptions.
"""
import pytest
from mp5d.exceptional import diagnostics, verification

@pytest.mark.parametrize("function,args,kwargs", [
    (diagnostics.root_separation,(lambda z:z,0.),{"radius":1.}),
    (diagnostics.zero_count,(lambda z:z,0.),{"radius":1.}),
    (verification.roots_in_disc,(lambda z:z,0.),{"radius":1.}),
    (verification.jordan_chain_matrix,([[0,1],[0,0]],[[1,0],[0,1]]),{}),
])
def test_retired_unsafe_api_cannot_gate_science(function,args,kwargs):
    with pytest.raises(RuntimeError, match="WITHDRAWN"):
        function(*args,**kwargs)
'''
    name='pyproject.toml'
    changes[name]=replace_once(existing(name),'packages = ["src/mp5d"]',
                              'packages = ["src/mp5d", "src/mp5d_science"]').encode()
    name='README.md'
    banner='''> **Scientific revalidation required — do not freeze or submit the legacy bound.**
> The hardening audit reproduced inversion-specific roots among all ten saved
> closest interactions. Corrected recurrence-free checks do not independently
> confirm the saved N2 branches. All numerical results described below remain
> historical, not newly certified. Current claim authority is
> `config/science/claims.json`; see `docs/SCIENCE_HARDENING.md`.
> The hyperboloidal module is a testbed, not a validated massive MP5D Solver F.

'''
    changes[name]=(banner+existing(name)).encode()
    name='results/status.json'
    historical=json.loads(existing(name))
    changes[name]=(json.dumps({
        'status': 'SCIENTIFIC_REVALIDATION_REQUIRED',
        'authority': 'HISTORICAL_ONLY',
        'primary_result': 'No refreshed exclusion bound is established by this installation.',
        'reason': 'All ten historical N2 roots failed the new cross-inversion and independent-solver audit.',
        'claim_registry': 'config/science/claims.json',
        'historical_status': historical,
        'installation_note': 'Historical data and original provenance are preserved below; installation generated no new scientific results.'
    }, indent=2, allow_nan=False)+'\n').encode()
    # Replace the permissive manual workflow instead of leaving a misleading green route.
    name='.github/workflows/research-validation.yml'
    existing(name)
    changes[name]=b'''name: Research validation (hardening)
on:
  workflow_dispatch:
permissions:
  contents: read
jobs:
  validate:
    runs-on: ubuntu-latest
    timeout-minutes: 360
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0
      - uses: actions/setup-python@v5
        with:
          python-version: '3.13.5'
      - run: bash scripts/reproduce_submission_science.sh "${RUNNER_TEMP}/mp5d-science"
'''
    for name,data in changes.items():
        target=(root/name).resolve()
        if not target.is_relative_to(root):raise ValueError('unsafe overlay path '+name)
        if target.is_symlink() or any(parent.is_symlink() for parent in (root/name).parents if parent!=root.parent):
            raise ValueError('refusing symlinked destination '+name)
        if name.endswith('.py'):ast.parse(data.decode(),filename=name)
    return changes


def atomic_write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else 0o644
    if path.suffix == '.sh':
        mode |= 0o111
    fd,temp=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent)
    os.fchmod(fd, mode)
    try:
        with os.fdopen(fd,'wb') as handle:handle.write(data)
        os.replace(temp,path)
    finally:
        if os.path.exists(temp):os.unlink(temp)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkout',type=Path)
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    root=args.checkout.resolve()
    try:
        changes=build_changes(root)
        for name in sorted(changes):print(('MODIFY ' if (root/name).exists() else 'ADD    ')+name)
        if args.dry_run:
            print(f'DRY RUN: {len(changes)} validated file changes; nothing written.');return 0
        previous={name:(root/name).read_bytes() if (root/name).exists() else None for name in changes}
        written=[]
        try:
            for name,data in changes.items():
                atomic_write(root/name,data);written.append(name)
        except BaseException:
            for name in reversed(written):
                if previous[name] is None:(root/name).unlink(missing_ok=True)
                else:atomic_write(root/name,previous[name])
            raise
        print('Applied. Review git diff; no commit, push or science-freeze tag was created.')
        return 0
    except (ValueError,OSError,subprocess.CalledProcessError,SyntaxError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);return 2
if __name__=='__main__':sys.exit(main())

"""Fail-closed provenance for code commits and post-commit research artifacts.

An artifact cannot meaningfully embed the hash of the Git commit that contains
that very artifact (a circular hash dependency). Freeze clean source first;
produce external/ignored artifacts afterwards with that exact source commit.
Never replace historical generating commits with a more convenient HEAD.
"""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata
from pathlib import Path
import json
import math
import platform
import re
import subprocess
import sys
import tempfile
import os

SCHEMA_VERSION = 1
COMMIT_RE = re.compile(r'^[0-9a-f]{40}$')


def canonical_json(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value) -> str:
    return sha256(canonical_json(value)).hexdigest()


def file_hash(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(['git', '-C', str(root), *args], text=True,
                                   stderr=subprocess.PIPE).strip()


def environment() -> dict:
    packages={}
    for name in ('numpy','scipy','sympy','mpmath','pytest','python-flint'):
        try: packages[name]=metadata.version(name)
        except metadata.PackageNotFoundError: packages[name]=None
    value={'python':platform.python_version(), 'implementation':platform.python_implementation(),
           'platform':platform.platform(), 'packages':packages}
    return value | {'sha256':digest(value)}


def source_snapshot(root: Path) -> dict:
    """Content digest includes source, tests, configs and benchmark fixtures, not outputs."""
    root=Path(root)
    files={}
    for directory in ('src','scripts','config','tests','data/regressions','lean','docs','.github'):
        base=root/directory
        if not base.exists(): continue
        for path in sorted(base.rglob('*')):
            if path.is_file() and not any(p.startswith('.') or p=='__pycache__' for p in path.relative_to(base).parts):
                files[path.relative_to(root).as_posix()]=file_hash(path)
    for name in ('requirements.lock', 'requirements-science.lock', 'requirements-arb.lock',
                 'pyproject.toml', 'pytest-science.ini', 'README.md'):
        path=root/name
        if path.exists():files[name]=file_hash(path)
    return {'files':files,'sha256':digest(files)}


def provenance(root: Path, *, allow_development: bool=False) -> dict:
    root=Path(root).resolve()
    snap=source_snapshot(root)
    try:
        commit=git(root,'rev-parse','HEAD')
        dirty=bool(git(root,'status','--porcelain','--untracked-files=all'))
    except (subprocess.CalledProcessError,FileNotFoundError):
        commit=None;dirty=True
    if not allow_development and (not commit or dirty):
        raise ValueError('authoritative runs require a clean committed checkout')
    return {'source_commit':commit,'source_tree_sha256':snap['sha256'],
            'source_files':snap['files'],'dirty_checkout':dirty,
            'authority':'DEVELOPMENT_ONLY' if dirty or not commit else 'COMMITTED_SOURCE_RUN',
            'environment':environment(), 'created_utc':datetime.now(timezone.utc).isoformat()}


def envelope(kind: str, payload: dict, root: Path, *, allow_development=False) -> dict:
    if not kind: raise ValueError('named result kind required')
    record={'schema_version':SCHEMA_VERSION,'kind':kind,
            'provenance':provenance(root,allow_development=allow_development),'payload':payload}
    return record | {'integrity_sha256':digest(record)}


def atomic_json(path: Path, value) -> None:
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    # Serializing first ensures invalid numerical data cannot truncate an old artifact.
    data=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'
    fd,temporary=tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as handle:
            handle.write(data);handle.flush();os.fsync(handle.fileno())
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def check_result(record: dict, *, expected_commit: str|None=None,
                 expected_source_digest: str|None=None, authoritative=False) -> list[str]:
    errors=[]
    if not isinstance(record,dict):return ['result must be an object']
    if record.get('schema_version')!=SCHEMA_VERSION:errors.append('unsupported schema_version')
    if not isinstance(record.get('kind'),str) or not record.get('kind'):errors.append('missing result kind')
    for key in ('provenance','payload'):
        if not isinstance(record.get(key),dict):errors.append('missing object '+key)
    checksum=record.get('integrity_sha256')
    try:
        actual=digest({k:v for k,v in record.items() if k!='integrity_sha256'})
        if checksum!=actual:errors.append('integrity hash mismatch')
    except (ValueError,TypeError):errors.append('non-JSON or non-finite result')
    prov=record.get('provenance',{})
    if not isinstance(prov,dict):return errors
    env=prov.get('environment',{})
    if not isinstance(env,dict) or env.get('sha256')!=digest({k:v for k,v in env.items() if k!='sha256'}):
        errors.append('environment hash missing or invalid')
    source_files=prov.get('source_files')
    if not isinstance(source_files,dict) or not source_files or prov.get('source_tree_sha256')!=digest(source_files):
        errors.append('source snapshot missing or invalid')
    if expected_source_digest and prov.get('source_tree_sha256')!=expected_source_digest:
        errors.append('result source digest differs from current checkout')
    if authoritative:
        if not COMMIT_RE.fullmatch(prov.get('source_commit') or ''):errors.append('exact source commit missing')
        if prov.get('dirty_checkout') is not False:errors.append('dirty checkout cannot be authoritative')
        if prov.get('authority')!='COMMITTED_SOURCE_RUN':errors.append('development result is not authoritative')
        if not expected_commit:errors.append('science-freeze commit must be specified externally')
    if expected_commit and prov.get('source_commit')!=expected_commit:
        errors.append('generating commit is not the exact requested science-freeze commit')
    return errors


def read_json(path: Path) -> dict:
    def invalid_constant(x): raise ValueError('non-finite JSON literal '+x)
    return json.loads(Path(path).read_text(),parse_constant=invalid_constant)

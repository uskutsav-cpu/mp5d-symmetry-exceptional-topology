"""Immutable run identity, exact task coverage, and auditable resumable receipts."""
from __future__ import annotations

import os
import subprocess
import time
import uuid
from pathlib import Path

import hashlib
import importlib.metadata
import json
import platform
import sys
import tempfile
import re
import traceback


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def digest(value):
    return hashlib.sha256(canonical_json(value)).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = canonical_json(value) + b"\n"
    fd, temporary = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(encoded)
            f.flush()
            os.fsync(f.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def environment():
    packages = sorted((d.metadata.get("Name", "unknown"), d.version)
                      for d in importlib.metadata.distributions())
    return {"python": sys.version, "platform": platform.platform(),
            "packages": [list(item) for item in packages],
            "threads": {k: os.environ.get(k) for k in
                        ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}}


AXES = {
    'cf_depth': [160, 240, 320, 480], 'angular_n': [24, 40, 64],
    'precision_dps': [30, 50, 80], 'radial_n': [160, 220, 280, 360],
    'contour_length': [40, 60, 90], 'scaling_angle_deg': [60, 65, 70],
    'continuation_step': [0.01, 0.005, 0.0025],
}
ARTIFACTS = {
    'benchmarks.json': ('benchmark-and-algebra-audit', 'SUBMISSION_BENCHMARK_PASS'),
    'candidates.json': ('candidate-regression-audit', 'POINTWISE_CHECKS_PASS'),
    'atlas.json': ('fresh-branch-atlas', 'VALIDATED'),
    'convergence.json': ('dangerous-point-convergence', 'PASS'),
    'boundary.json': ('boundary-refinement', 'LOCALLY_STABLE_NUMERICAL_MINIMUM'),
    'robustness.json': ('higher-mode-probe', 'VALIDATED_OUT_OF_DOMAIN_PROBE'),
    'bibliography.json': ('bibliography-convention-audit', 'ALL_SUBMISSION_INPUTS_VERIFIED'),
    'lean.json': ('lean-kernel-build', 'BUILD_PASS'),
}


class Blocked(RuntimeError):
    """A scientific or provenance prerequisite has not been established."""


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(['git', '-C', str(root), *args], text=True,
                                   stderr=subprocess.PIPE).strip()


def validate_plan(plan: dict):
    if plan.get('axes') != AXES:
        raise Blocked('The merged seven-axis ladders must match exactly; no dropped or substituted rung')
    h = plan['higher_modes']
    if h['overtones'] != [4, 5] or h['ell_offsets'] != [0, 2, 4]:
        raise Blocked('Higher-mode coverage must include N=4,5 and ell offsets 0,2,4')
    if plan['formal_domain_proposal']['overtones'] != [0, 1, 2, 3]:
        raise Blocked('Formal overtone domain changed; an explicit new campaign version is required')
    if plan['adaptive_boundary']['independent_continuation_directions'] != ['forward', 'reverse']:
        raise Blocked('Both genuinely distinct approaches are required')


def snapshot(root: Path, strict: bool) -> dict:
    """Hash every tracked input, including bibliography and all published tables."""
    root = root.resolve()
    try:
        if Path(git(root, 'rev-parse', '--show-toplevel')).resolve() != root:
            raise Blocked('Source directory must be the worktree root')
        commit = git(root, 'rev-parse', 'HEAD')
        dirty = bool(git(root, 'status', '--porcelain', '--untracked-files=all'))
        if git(root, 'submodule', 'status'):
            raise Blocked('Unimplemented submodule provenance; refusing authoritative source')
        names = git(root, 'ls-files', '-z', '--cached', '--others', '--exclude-standard').split('\0')
    except (subprocess.CalledProcessError, FileNotFoundError):
        if strict:
            raise Blocked('Authoritative execution needs an existing clean committed checkout') from None
        commit, dirty = None, True
        names = [p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()
                 and not any(x in {'.git', '.venv', '__pycache__', '.pytest_cache', '.lake', '.ruff_cache'}
                             for x in p.relative_to(root).parts)]
    if strict and dirty:
        raise Blocked('Source is dirty or has untracked inputs: commit it before reproduction')
    hashes = {}
    for name in sorted(n for n in names if n):
        path = root / name
        if path.is_symlink() or not path.is_file():
            raise Blocked('Missing or nonregular source input: ' + name)
        hashes[name] = file_hash(path)
    if not hashes:
        raise Blocked('Empty source snapshot')
    return {'source_commit': commit, 'source_files': hashes,
            'source_tree_sha256': digest(hashes), 'dirty_checkout': dirty}


def coverage(expected: list[str], tasks: dict) -> dict:
    if not expected or len(expected) != len(set(expected)):
        raise Blocked('Expected task list must be nonempty and unique')
    missing = sorted(set(expected) - tasks.keys())
    unexpected = sorted(tasks.keys() - set(expected))
    failed = sorted(k for k in expected if k in tasks and tasks[k].get('status') != 'PASS')
    return {'expected': expected, 'completed': sorted(tasks), 'missing': missing,
            'unexpected': unexpected, 'failed': failed, 'complete': not (missing or unexpected or failed)}


def require_coverage(payload: dict):
    reported = payload.get('coverage', {})
    recomputed = coverage(reported.get('expected', []), payload.get('tasks', {}))
    if recomputed != reported or not recomputed['complete']:
        raise Blocked('Missing, failed, unexpected, or manipulated task coverage')


def payload(status: str, expected: list[str], tasks: dict, **extra) -> dict:
    c = coverage(expected, tasks)
    return {'status': status if c['complete'] else 'INCOMPLETE_OR_FAILED',
            'coverage': c, 'tasks': tasks, **extra}


class Store:
    """Single-writer evidence directory, outside source; exact-context resume only.

    A task hash includes source, dependencies/environment, configuration, actual
    inputs/seeds, and requested algorithm. Scientific failure is cached as failure.
    """
    def __init__(self, root: Path, output: Path, config: dict, plan: dict, *,
                 development=False, resume=False, retry_failed=False):
        validate_plan(plan)
        self.root, self.output = root.resolve(), output.resolve()
        self.development, self.retry_failed = development, retry_failed
        if self.output == self.root or self.root in self.output.parents:
            raise Blocked('Put generated evidence OUTSIDE the source checkout')
        self.output.mkdir(parents=True, exist_ok=True)
        self.lock = self.output / '.writer.lock'
        self.fd = None
        try:
            self.fd = os.open(self.lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        except FileExistsError:
            raise Blocked('Another writer or an interrupted run owns .writer.lock; inspect before removing it') from None
        os.write(self.fd, canonical_json({'pid': os.getpid(), 'time': time.time()}))
        try:
            source_state = snapshot(self.root, not development)
            if not development:
                lock = self.root / 'requirements.lock'
                if not lock.is_file():
                    raise Blocked('Authoritative environment requires the committed requirements.lock')
                mismatches = []
                pinned_names = set()
                for line in lock.read_text().splitlines():
                    match = re.match(r'^([A-Za-z0-9_.-]+)==([^;\s]+)\s*$', line.strip())
                    if match:
                        pinned_names.add(match[1].lower().replace("_", "-"))
                        try:
                            actual = importlib.metadata.version(match[1])
                        except importlib.metadata.PackageNotFoundError:
                            actual = 'NOT_INSTALLED'
                        if actual != match[2]:
                            mismatches.append(f'{match[1]}: expected {match[2]}, found {actual}')
                if not {'numpy', 'scipy', 'mpmath'} <= pinned_names:
                    raise Blocked('The numerical numpy/scipy/mpmath versions must be explicitly pinned in the lock file')
                if mismatches:
                    raise Blocked('Pinned environment mismatch: ' + '; '.join(mismatches))
            identity = {**source_state, 'environment': environment(),
                        'authority': 'DEVELOPMENT_ONLY' if development else 'COMMITTED_SOURCE_RUN',
                        'config_sha256': digest(config), 'research_plan_sha256': digest(plan)}
            self.identity = identity
            self.context_hash = digest(identity)
            manifest = self.output / 'run.json'
            if manifest.exists():
                if not resume:
                    raise Blocked('Run directory already used: choose a fresh directory or --resume')
                self.run = read_json(manifest)
                if self.run['identity'] != identity or self.run['context_sha256'] != self.context_hash:
                    raise Blocked('Cannot resume across changed source/configuration/environment/authority')
            else:
                if resume:
                    raise Blocked('No run.json to resume')
                self.run = {'identity': identity, 'context_sha256': self.context_hash,
                            'run_id': str(uuid.uuid4()), 'started_unix': time.time()}
                atomic_json(manifest, self.run)
            (self.output / 'tasks').mkdir(exist_ok=True)
        except BaseException:
            self.close()
            raise

    def close(self):
        if self.fd is not None:
            os.close(self.fd)
            self.fd = None
            self.lock.unlink(missing_ok=True)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def unchanged(self):
        observed = snapshot(self.root, not self.development)
        if any(self.identity[k] != v for k, v in observed.items()) or environment() != self.identity['environment']:
            raise Blocked('Source or environment changed during execution; refuse to publish receipts')

    def task(self, request: dict, compute) -> dict:
        key = digest({'context': self.context_hash, 'request': request})
        path = self.output / 'tasks' / (key + '.json')
        record = None
        if path.exists():
            old = read_json(path)
            check = old.pop('receipt_sha256', None)
            if check != digest(old) or old.get('request') != request or old.get('context_sha256') != self.context_hash:
                raise Blocked('Corrupt/mismatched task receipt: ' + key)
            if not self.retry_failed or old['result'].get('status') == 'PASS':
                record = old
        if record is None:
            started = time.time()
            try:
                result = compute()
                canonical_json(result)
                if not isinstance(result, dict) or 'status' not in result:
                    raise ValueError('Task lacks an explicit status')
            except Exception as exc:
                result = {'status': 'FAILED', 'error_type': type(exc).__name__, 'reason': str(exc),
                          'traceback': traceback.format_exc(limit=8)}
            record = {'request': request, 'result': result, 'context_sha256': self.context_hash,
                      'started_unix': started, 'finished_unix': time.time()}
            atomic_json(path, {**record, 'receipt_sha256': digest(record)})
        return {**record['result'], 'receipt_key': key}

    def write(self, name: str, kind: str, data: dict, dependencies=()) -> dict:
        if Path(name).name != name:
            raise Blocked('Artifact name must be a basename')
        self.unchanged()
        record = {'schema_version': 2, 'kind': kind, 'provenance': self.run,
                  'payload': data, 'dependencies': {n: file_hash(self.output / n) for n in dependencies},
                  'task_receipts': {p.name: file_hash(p) for p in sorted((self.output / 'tasks').glob('*.json'))},
                  'assets': {p.relative_to(self.output).as_posix(): file_hash(p)
                             for directory in ('logs', 'lean-worktree/verified-objects')
                             for p in sorted((self.output / directory).rglob('*')) if p.is_file()},
                  'finished_unix': time.time()}
        record['integrity_sha256'] = digest(record)
        atomic_json(self.output / name, record)
        return record

    def verify(self, record: dict):
        if record.get('schema_version') != 2:
            raise Blocked('Only campaign evidence schema version 2 is accepted')
        if record.get('integrity_sha256') != digest({k: v for k, v in record.items() if k != 'integrity_sha256'}):
            raise Blocked('Artifact integrity mismatch')
        if record.get('provenance') != self.run:
            raise Blocked('Mixed source, environment, or run provenance')
        if record.get('finished_unix', 0) < self.run['started_unix']:
            raise Blocked('Artifact predates its source run')
        for name, checksum in record.get('assets', {}).items():
            path = (self.output / name).resolve()
            if self.output not in path.parents or not path.is_file() or file_hash(path) != checksum:
                raise Blocked('Build log or compiled proof object changed: ' + name)
        for name, checksum in record.get('dependencies', {}).items():
            if Path(name).name != name or file_hash(self.output / name) != checksum:
                raise Blocked('Dependency changed: ' + name)
        for name, checksum in record.get('task_receipts', {}).items():
            if Path(name).name != name or file_hash(self.output / 'tasks' / name) != checksum:
                raise Blocked('Raw evidence changed: ' + name)

    def read(self, name: str, accepted: str | None = None) -> dict:
        try:
            record = read_json(self.output / name)
        except OSError:
            raise Blocked('Missing prerequisite: ' + name) from None
        self.verify(record)
        data = record['payload']
        if accepted:
            if data.get('status') != accepted:
                raise Blocked(f'{name}: {data.get("status")} is not {accepted}')
            require_coverage(data)
        return data

    def blocked(self, name: str, reason: str, expected=None):
        kind = ARTIFACTS.get(name, ('campaign-diagnostic', ''))[0]
        ids = expected or ['prerequisite']
        return self.write(name, kind, payload('UNREACHABLE', ids,
            {k: {'status': 'BLOCKED', 'reason': reason} for k in ids},
            reason=reason, continuum_qnm_theorem=False))

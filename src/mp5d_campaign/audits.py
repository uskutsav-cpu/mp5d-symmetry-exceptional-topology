"""Benchmark, bibliography, optional ball arithmetic, and actual Lean builds."""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path

import mpmath as mp

from mp5d_science.convergence import Resolution

from .core import ARTIFACTS, Blocked, file_hash, payload, read_json
from .numerics import z


def command(command, cwd, log, timeout):
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    try:
        with log.open('w', encoding='utf-8') as handle:
            result = subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT,
                                    timeout=timeout, check=False)
        code, reason = result.returncode, None
    except (OSError, subprocess.TimeoutExpired) as exc:
        code, reason = None, str(exc)
        if not log.exists():
            log.write_text(reason + '\n')
    return {'status': 'PASS' if code == 0 else 'FAILED', 'command': list(command),
            'returncode': code, 'reason': reason, 'log_name': log.name,
            'log_sha256': file_hash(log), 'started_unix': started, 'finished_unix': time.time()}


def benchmarks(campaign):
    store, cfg = campaign.store, campaign.config
    static = read_json(store.root / cfg['static_references'])['records']
    rotating = read_json(store.root / cfg['rotating_references'])['records']
    refs = static + rotating
    expected, tasks = [], {}
    base = Resolution(**cfg['benchmark_resolution'])
    for ref in refs:
        with mp.workdps(100):
            scale = 2 * mp.pi if ref['normalization'].startswith('omega_tilde=') else 1
            if ref['normalization'] not in ('omega_tilde=2*pi*omega; rh=1', 'omega; rh=1', 'omega; M=1'):
                raise Blocked('Unrecognized benchmark frequency normalization')
            seed = [str(mp.mpf(ref['frequency'][0]) / scale), str(mp.mpf(ref['frequency'][1]) / scale)]
        for kind, field, settings in [('A', 'cf_depth', cfg['benchmark_depths']),
                                      ('C', 'radial_n', cfg['benchmark_radial_n'])]:
            for setting in settings:
                key = ref['id'] + f'/{kind}/{field}={setting}'
                expected.append(key)
                result = campaign.physics.solve(kind, ref['point'], ref['sector'], ref['N'], seed,
                                                 replace(base, **{field: setting}))
                error = abs(z(result['omega']) - z(seed)) if 'omega' in result else None
                tolerance = ref.get('benchmark_tolerance', cfg['benchmark_tolerance'])
                tasks[key] = {'status': 'PASS' if result['status'] == 'PASS' and error is not None and error <= tolerance else 'FAILED',
                    'reference': ref, 'solver': result, 'absolute_error': error, 'tolerance': tolerance,
                    'published_rounding_is_not_the_entire_method_error': True}
    # Exact symbolic/finite-algebra checks remain separate from numerical root agreement.
    key = 'foundation-submission-regressions'
    expected.append(key)
    paths = store.root / 'tests/submission'
    tasks[key] = command([sys.executable, '-m', 'pytest', '-q', str(paths), '-m', 'not slow', '-o', 'addopts=', '-o', 'cache_dir=' + str(store.output / 'pytest-cache')],
                         store.root, store.output / 'logs/benchmark-algebra-tests.log', cfg['pytest_timeout_seconds'])
    data = payload('SUBMISSION_BENCHMARK_PASS', expected, tasks,
                   checked_sources=[r['source_id'] for r in refs],
                   scope='CONFIGURED_NUMERICAL_REFERENCES_AND_FINITE_ALGEBRA_ONLY',
                   complete_literature_audit_established=False,
                   configured_reference_count=len(refs),
                   continuum_certificate=False)
    store.write('benchmarks.json', ARTIFACTS['benchmarks.json'][0], data)
    return data


def safe_source_file(root, name):
    if not isinstance(name, str) or not name:
        raise Blocked('Missing audited evidence path')
    path = (root / name).resolve()
    if root.resolve() not in path.parents or not path.is_file() or (root / name).is_symlink():
        raise Blocked('Primary evidence must be a regular, hashed file within committed source')
    return path


def bibliography(campaign):
    store, cfg = campaign.store, campaign.config
    review_manifest = read_json(store.root / cfg['source_reviews'])
    reviews = review_manifest.get('reviews', [])
    by_source = {}
    for review in reviews:
        source = review.get('source_id')
        if source in by_source:
            raise Blocked('Duplicate bibliography review: ' + str(source))
        by_source[source] = review
    required = list(cfg['required_primary_sources'])
    input_ids = []
    for name in (cfg['static_references'], cfg['rotating_references']):
        for row in read_json(store.root / name)['records']:
            input_ids.append(row['id'])
            if row['source_id'] not in required:
                required.append(row['source_id'])
    csv_file = store.root / cfg['bibliography_file']
    missing_csv = not csv_file.is_file()
    if not missing_csv:
        with csv_file.open(newline='', encoding='utf-8-sig') as handle:
            bibliography_rows = list(csv.DictReader(handle))
        # Every indexed citation needs either a review or an explicit, justified
        # not-used-by-submission disposition. An old metadata_verified flag alone
        # never counts as a fresh convention or numeric transcription review.
        for index, row in enumerate(bibliography_rows):
            key = row.get('key') or row.get('citation_key') or row.get('id') or f'csv-row-{index + 2}'
            if key not in required:
                required.append(key)
    else:
        bibliography_rows = []
    tasks = {}
    all_reviewed_inputs = []
    for source in required:
        review = by_source.get(source)
        try:
            if not review:
                raise Blocked('No fresh primary-source audit')
            notes = safe_source_file(store.root, review.get('verification_notes_file'))
            if file_hash(notes) != review.get('verification_notes_sha256'):
                raise Blocked('Verification notes hash mismatch')
            if review.get('disposition') == 'NOT_USED_BY_SUBMISSION':
                if source in cfg['required_primary_sources'] or not review.get('exclusion_reason'):
                    raise Blocked('A required source cannot be silently excluded')
            else:
                for field in ('verified_metadata', 'verified_conventions', 'verified_numeric_transcription',
                              'verified_domain_of_applicability'):
                    if review.get(field) is not True:
                        raise Blocked('Unverified primary-source dimension: ' + field)
                if not str(review.get('primary_source_url', '')).startswith('https://'):
                    raise Blocked('Missing primary-source attribution')
                evidence = safe_source_file(store.root, review.get('primary_evidence_file'))
                if file_hash(evidence) != review.get('primary_evidence_sha256'):
                    raise Blocked('Primary evidence hash mismatch')
                all_reviewed_inputs.extend(review.get('submission_inputs', []))
            tasks[source] = {'status': 'PASS', 'review': review}
        except (Blocked, OSError) as exc:
            tasks[source] = {'status': 'FAILED', 'reason': str(exc)}
    inventory_ok = (review_manifest.get('complete_submission_input_inventory') is True and
                    set(input_ids) <= set(all_reviewed_inputs) and not missing_csv)
    tasks['complete-input-inventory'] = {'status': 'PASS' if inventory_ok else 'FAILED',
        'configured_inputs': input_ids, 'unreviewed_inputs': sorted(set(input_ids) - set(all_reviewed_inputs)),
        'bibliography_csv_missing': missing_csv,
        'complete_submission_inventory_attested': review_manifest.get('complete_submission_input_inventory') is True}
    data = payload('ALL_SUBMISSION_INPUTS_VERIFIED', required + ['complete-input-inventory'], tasks,
                   audit_dimensions=['metadata', 'units/sign/normalization', 'numeric transcription', 'applicability', 'input inventory'])
    store.write('bibliography.json', ARTIFACTS['bibliography.json'][0], data)
    return data


def lean_build(store, config):
    source = store.root / 'lean'
    expected = ['toolchain', 'source-holes', 'dependency-resolution', 'lake-build', 'all-lean-files']
    tasks = {}
    files = sorted(source.rglob('*.lean')) if source.exists() else []
    files = [p for p in files if '.lake' not in p.parts]
    holes = []
    for path in files:
        text = re.sub(r'/\-.*?\-/', '', path.read_text(), flags=re.S)
        text = re.sub(r'--[^\n]*', '', text)
        if re.search(r'\b(sorry|sorryAx|admit|axiom|unsafeCast)\b', text):
            holes.append(path.relative_to(source).as_posix())
    tasks['source-holes'] = {'status': 'PASS' if files and not holes else 'FAILED',
                             'files': [p.relative_to(source).as_posix() for p in files], 'holes': holes,
                             'limitations': 'Lexical preflight is additional to, not a replacement for, the kernel build'}
    executable = shutil.which('lake')
    if not executable:
        for key in expected:
            tasks.setdefault(key, {'status': 'BLOCKED', 'reason': 'lake executable is unavailable'})
        data = payload('BUILD_PASS', expected, tasks, compilation_checked=False,
                       proof_scope='ONLY_THE_COMPILED_FINITE_STATEMENTS; NO_CONTINUUM_OPERATOR_THEOREM')
        data['status'] = 'NOT_BUILT'
        store.write('lean.json', ARTIFACTS['lean.json'][0], data)
        return data
    build = store.output / 'lean-worktree'
    if build.exists():
        shutil.rmtree(build)
    shutil.copytree(source, build, ignore=shutil.ignore_patterns('.lake'))
    timeout = config['lean_timeout_seconds']
    logdir = store.output / 'logs'
    tasks['toolchain'] = command([executable, 'env', 'lean', '--version'], build, logdir / 'lean-version.log', timeout)
    tasks['toolchain']['requested_toolchain'] = (source / 'lean-toolchain').read_text().strip()
    expected_version = tasks['toolchain']['requested_toolchain'].split(':')[-1].lstrip('v')
    if expected_version not in (logdir / 'lean-version.log').read_text():
        tasks['toolchain']['status'] = 'FAILED'
        tasks['toolchain']['reason'] = 'Actual Lean version does not match the pinned toolchain'
    tasks['dependency-resolution'] = command([executable, 'update'], build, logdir / 'lean-update.log', timeout)
    tasks['lake-build'] = command([executable, 'build'], build, logdir / 'lean-build.log', timeout)
    individual = {}
    if tasks['lake-build']['status'] == 'PASS' and not holes:
        for original in files:
            relative = original.relative_to(source)
            output = build / 'verified-objects' / relative.with_suffix('.olean')
            output.parent.mkdir(parents=True, exist_ok=True)
            key = relative.as_posix()
            report = command([executable, 'env', 'lean', '-o', str(output), str(relative)],
                             build, logdir / ('lean-' + key.replace('/', '_') + '.log'), timeout)
            if not output.is_file():
                report['status'] = 'FAILED'
            else:
                report['olean_sha256'] = file_hash(output)
            individual[key] = report
    tasks['all-lean-files'] = {'status': 'PASS' if len(individual) == len(files) and files and
                              all(r['status'] == 'PASS' for r in individual.values()) else 'FAILED',
                              'individual_builds': individual}
    data = payload('BUILD_PASS', expected, tasks,
                   compilation_checked=True, proof_scope='COMPILED_FINITE_STATEMENTS_ONLY',
                   full_bijection_operator_theorem_established=False)
    store.write('lean.json', ARTIFACTS['lean.json'][0], data)
    return data


def arb_coefficients(store):
    """Actually execute generic rotating ball arithmetic when python-flint exists.

    This checks coefficient evaluation only. It does not certify a root, a
    continuum error bound, a generic rotating Krawczyk inclusion, or a no-EP claim.
    """
    def compute():
        try:
            from flint import acb, ctx
        except ImportError:
            return {'status': 'NOT_EXECUTED', 'reason': 'python-flint is not installed'}
        from mp5d_science.radial_polynomial import coefficient_formula, horizons
        old = ctx.prec
        try:
            ctx.prec = 192
            a, b, mass = acb('0.31'), acb('0.13'), acb('0.4')
            w, lam = acb('2.0', '-0.5'), acb('24.0', '0.1')
            p, m = horizons(a, b, acb(1), lambda value: value.sqrt())
            k = (w * w - mass * mass).sqrt()
            A, B, C, remainder = coefficient_formula(a=a, b=b, mu=mass, m1=2, m2=0, omega=w, lam=lam, p=p, m=m, outgoing_k=k, imaginary_unit=acb(0, 1))
            finite = all(value.is_finite() for group in (A, B, C, remainder) for value in group)
            return {'status': 'PASS' if finite else 'FAILED', 'precision_bits': 192,
                    'parameters': {'a': '0.31', 'b': '0.13', 'mu': '0.4', 'M': '1', 'm1': 2, 'm2': 0},
                    'A': list(map(str, A)), 'B': list(map(str, B)), 'C': list(map(str, C)),
                    'division_remainder': list(map(str, remainder)),
                    'scope': 'FINITE_GENERIC_ROTATING_COEFFICIENT_EVALUATION_ONLY',
                    'root_inclusion_certified': False, 'continuum_error_bound': False}
        finally:
            ctx.prec = old
    result = store.task({'algorithm': 'generic-rotating-acb-coefficient-evaluation/v1'}, compute)
    store.write('arb.json', 'optional-finite-coefficient-ball-audit', result)
    return result

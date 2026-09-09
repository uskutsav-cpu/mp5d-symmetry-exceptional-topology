#!/usr/bin/env python3
"""Re-run standalone development checks, preserving failures and exit codes.

This is not the upstream full-repository submission workflow. No dependencies
are installed, no network calls are made and no Git operations are performed.
The overall exit is 2 while scientific or authority gates remain blocked.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
OVERLAY = HERE / 'overlay'
BASE = '72abf0bf68d41d07de0fbd5056c40b7a713a6cb5'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'evidence')
    parser.add_argument('--resume-candidates', action='store_true',
                        help='Reuse identical-source completed candidate records, with strict validation')
    args = parser.parse_args()
    out = args.output.resolve()
    if out.is_relative_to(OVERLAY):
        parser.error('evidence output must be outside the overlay source tree')
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ | {'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1',
                        'PYTHONPATH': str(OVERLAY / 'src')}
    py = sys.executable

    def script(name: str, *arguments: str) -> list[str]:
        return [py, str(OVERLAY / 'scripts' / name), *map(str, arguments)]

    tasks = [
        ('new-tests', [py, '-m', 'pytest', '-c', str(OVERLAY / 'pytest-science.ini'),
                        str(OVERLAY / 'tests/submission'), '-q', '--strict-markers',
                        '--junitxml='+str(out / 'new-tests.xml')]),
        ('bundle-installer-tests', [py, '-m', 'pytest', str(HERE / 'test_bundle_installer.py'),
                                    '-q', '--junitxml='+str(out / 'bundle-installer-tests.xml')]),
        ('candidates', script('run_candidate_audit.py', '--development', '--dps', '45',
                              '--output', out / 'candidates.json')),
        ('benchmarks', script('run_benchmark_audit.py', '--development',
                              '--output', out / 'benchmarks.json')),
        ('boundary', script('run_boundary_or_ladders.py', 'boundary', '--development',
                            '--output', out / 'boundary.json')),
        ('convergence', script('run_boundary_or_ladders.py', 'ladders', '--development',
                               '--output', out / 'convergence.json')),
        ('lean', script('check_lean.py', '--development', '--output', out / 'lean.json')),
        ('claims', script('check_claim_consistency.py')),
    ]

    if args.resume_candidates:
        tasks = [(name, [py, str(HERE / 'RESUME_CANDIDATES.py'), '--output', str(out / 'candidates.json')]
                  if name == 'candidates' else command) for name, command in tasks]

    def run(task: tuple[str, list[str]]) -> dict:
        name, command = task
        started = time.monotonic()
        with (out / (name+'.log')).open('w') as log:
            try:
                result = subprocess.run(command, cwd=HERE, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=3600, check=False)
                code = result.returncode
                failure = None
            except (OSError, subprocess.TimeoutExpired) as error:
                code = None
                failure = str(error)
                log.write('\nExecution error: '+failure+'\n')
        record = {'name': name, 'command': command, 'returncode': code,
                  'runtime_seconds': time.monotonic()-started, 'execution_error': failure}
        print(f'{name}: returncode={code}', flush=True)
        return record

    with ThreadPoolExecutor(max_workers=2) as executor:
        records = list(executor.map(run, tasks))
    receipts = [out / (name+'.json') for name in ('candidates','benchmarks','boundary','convergence','lean')]
    records.append(run(('schema', script('check_result_schema.py', *receipts))))
    records.append(run(('provenance-authority', script('check_provenance.py', *receipts,
                                                      '--expected-commit', BASE))))
    records.append(run(('freeze', script('check_science_freeze.py', '--artifact-dir', out,
                                         '--expected-commit', BASE, '--output', out / 'freeze.json'))))
    software_pass = all(r['returncode'] == 0 for r in records if r['name'] in
                        {'new-tests','bundle-installer-tests','claims','schema'})
    all_pass = all(r['returncode'] == 0 for r in records)
    report = {'created_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'STANDALONE_DEVELOPMENT_CHECKS_NOT_FULL_UPSTREAM_INTEGRATION',
              'software_and_schema_checks_pass': software_pass,
              'all_scientific_and_authority_gates_pass': all_pass,
              'overall_status': 'ALL_CHECKS_PASS' if all_pass else 'SCIENCE_FREEZE_BLOCKED',
              'checks': records}
    (out / 'execution_summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return 0 if all_pass else 2


if __name__ == '__main__':
    raise SystemExit(main())

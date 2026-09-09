#!/usr/bin/env python3
"""Resume a development candidate audit only under identical source and settings.

Completed records are reused only after integrity, source, environment, fixture,
ordering and per-solver settings checks. This cannot promote results to a freeze.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent / 'overlay'
sys.path.insert(0, str(ROOT / 'src'))
from mp5d_science.provenance import (
    read_json, check_result, source_snapshot, environment, envelope, atomic_json,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    old = read_json(args.output)
    errors = check_result(old, expected_source_digest=source_snapshot(ROOT)['sha256'])
    if old['provenance']['authority'] != 'DEVELOPMENT_ONLY':
        errors.append('resume helper accepts only development receipts')
    if old['provenance']['environment']['sha256'] != environment()['sha256']:
        errors.append('environment changed')
    fixture = read_json(ROOT / 'data/regressions/ten_interactions.json')
    payload = old['payload']
    results = payload.get('results', [])
    if payload.get('fixture_source') != fixture['source']:
        errors.append('fixture provenance changed')
    if payload.get('completed') != len(results) or payload.get('expected') != len(fixture['records']):
        errors.append('invalid completed/expected counts')
    if len(results) > len(fixture['records']):
        errors.append('too many records')
    for result, item in zip(results, fixture['records']):
        if result['id'] != item['id']:
            errors.append('record order differs from fixture')
        for branch in item['branches']:
            rows = [r for r in result['records'] if r['branch'] == branch]
            expected = [{'solver': 'A', 'depth': n, 'dps': 45} for n in (160, 240, 320)]
            expected += [{'solver': 'C', 'radial_n': n} for n in (160, 220, 280)]
            if [r['requested'] for r in rows] != expected:
                errors.append('completed record settings differ')
    if errors:
        raise ValueError('; '.join(errors))
    spec = importlib.util.spec_from_file_location('candidate_runner', ROOT / 'scripts/run_candidate_audit.py')
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    for item in fixture['records'][len(results):]:
        result = runner.audit_candidate(item, [160, 240, 320], 45, [160, 220, 280])
        results.append(result)
        payload.update(status='RUNNING', completed=len(results), results=results)
        atomic_json(args.output, envelope('candidate-regression-audit', payload, ROOT, allow_development=True))
        print(result['id'], result['status'], flush=True)
    valid = all(r['status'] == 'POINTWISE_ROOT_CHECKS_PASSED' for r in results)
    payload['status'] = 'POINTWISE_CHECKS_PASS' if valid else 'SCIENTIFIC_REVALIDATION_REQUIRED'
    atomic_json(args.output, envelope('candidate-regression-audit', payload, ROOT, allow_development=True))
    return 0 if valid else 2


if __name__ == '__main__':
    raise SystemExit(main())

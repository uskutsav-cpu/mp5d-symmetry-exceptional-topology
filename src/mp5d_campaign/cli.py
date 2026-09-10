"""Run with PYTHONPATH=src python -m mp5d_campaign.cli --help."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from mp5d_science.convergence import Resolution

from . import audits
from .campaign import Campaign
from .core import ARTIFACTS, Blocked, Store, atomic_json, read_json, validate_plan
from .numerics import Physics, distinct, z
from .release import prepare_claims, release


def diagnose(store, config):
    physics = Physics(store, config)
    res = Resolution(cf_depth=240, angular_n=40, precision_dps=35, radial_n=220,
                     contour_length=60, scaling_angle_deg=70, continuation_step=.01)
    sector = {'m1': 2, 'm2': 0, 'ell': 4}
    roots = []
    for n in (2, 3):
        seed = 2.5 - (n + .5) * .70710678j
        a = physics.solve('A', [0., 0., 0.], sector, n, seed, res)
        c = physics.solve('C', [0., 0., 0.], sector, n, seed, res)
        roots.append({'algorithm_inversion': n, 'A': a, 'C': c})
    identity = 'NOT_CHECKED'
    if all(r['A']['status'] == r['C']['status'] == 'PASS' for r in roots):
        try:
            distinct([r['A']['omega'] for r in roots], [abs(z(r['A']['omega']) - z(r['C']['omega'])) for r in roots])
            identity = 'DISTINCT_ROOTS_BUT_OVERTONE_LABELS_NOT_CERTIFIED'
        except Blocked:
            identity = 'DUPLICATE_OR_UNRESOLVED_ROOTS'
    result = {'status': 'DEVELOPMENT_DIAGNOSTIC', 'roots': roots, 'identity_result': identity,
              'validates_atlas': False, 'science_ready': False}
    store.write('diagnostics.json', 'inversion-is-not-identity-regression', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['plan', 'run', 'diagnose', 'prepare-claims', 'release'])
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--config', default='config/campaign/campaign.json')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--development', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--stages', default='all', help='Comma-separated stages, or all')
    parser.add_argument('--point-limit', type=int)
    parser.add_argument('--tag', help='Explicit local numerical-v... release tag; never pushed automatically')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        config = read_json(root / args.config)
        plan = read_json(root / 'config/science/research_plan.json')
        validate_plan(plan)
        if args.command == 'plan':
            print(json.dumps({'axes': plan['axes'], 'sectors': config['sectors'], 'scope': config['scope_note'],
                              'formal_overtones': [0, 1, 2, 3], 'higher_modes': plan['higher_modes'],
                              'quasiresonant_scope': 'EXCLUDED', 'science_ready': False}, indent=2))
            return 0
        if not args.output:
            parser.error('--output is required and must be outside the source checkout')
        if args.point_limit is not None and args.point_limit < 1:
            parser.error('--point-limit must be positive')
        if args.tag and args.command != 'release':
            parser.error('--tag is only valid for release')
        with Store(root, args.output, config, plan, development=args.development or args.command == 'diagnose',
                   resume=args.resume, retry_failed=args.retry_failed) as store:
            if args.command == 'diagnose':
                print(json.dumps(diagnose(store, config), indent=2))
                return 0
            if args.command == 'prepare-claims':
                print(json.dumps(prepare_claims(store, config, plan), indent=2))
                return 0
            if args.command == 'release':
                print(json.dumps(release(store, config, plan, args.tag), indent=2))
                return 0
            campaign = Campaign(store, config, plan)
            stages = {
                'benchmarks': lambda: audits.benchmarks(campaign),
                'atlas': lambda: campaign.atlas(args.point_limit),
                'candidates': campaign.candidates,
                'convergence': campaign.convergence,
                'boundary': campaign.boundary,
                'robustness': campaign.robustness,
                'bibliography': lambda: audits.bibliography(campaign),
                'lean': lambda: audits.lean_build(store, config),
                'arb': lambda: audits.arb_coefficients(store),
            }
            selected = list(stages) if args.stages == 'all' else args.stages.split(',')
            if not selected or any(stage not in stages for stage in selected) or len(selected) != len(set(selected)):
                raise Blocked('Unknown or duplicate stage names')
            statuses = {}
            for name in selected:
                print('RUNNING ' + name, flush=True)
                try:
                    result = stages[name]()
                except Exception as exc:
                    record = store.blocked(name + '.json', type(exc).__name__ + ': ' + str(exc))
                    result = record['payload']
                statuses[name] = result['status']
                print(name + ': ' + str(result['status']), flush=True)
                atomic_json(store.output / 'progress.json', {'stages': statuses, 'science_ready': False,
                            'authority': store.identity['authority'], 'run_id': store.run['run_id']})
            all_required = all(statuses.get(Path(name).stem) == accepted for name, (_, accepted) in ARTIFACTS.items())
            print(json.dumps({'stages': statuses, 'all_required_stage_statuses_pass': all_required,
                              'science_ready': False, 'claim_promotion_is_a_separate_verified_step': True}, indent=2))
            return 0 if all_required else 2
    except (Blocked, OSError, ValueError, KeyError) as exc:
        print('BLOCKED: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())

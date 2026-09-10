"""Two-phase claim promotion; publication outputs only after strict acceptance."""
from __future__ import annotations

import copy
import math
import re

from .core import ARTIFACTS, Blocked, atomic_json, file_hash, git, read_json, require_coverage
from .numerics import z, sid, pid
from .campaign import pairs_at


def _global_claims_safe(value):
    if isinstance(value, dict):
        identity = str(value.get('id', value.get('claim_id', ''))).lower()
        if ('global' in identity and ('ep' in identity or 'exception' in identity)):
            if value.get('active') is True or value.get('status') not in (None, 'NOT_ESTABLISHED', 'UNSUPPORTED'):
                return False
        return all(_global_claims_safe(v) for v in value.values())
    if isinstance(value, list):
        return all(_global_claims_safe(v) for v in value)
    return True


def check_physical_row(row, tolerance):
    if row.get('status') != 'PASS' or not row.get('solver_reports'):
        raise Blocked('Missing physical A/C evidence')
    if len(row['overtones']) != len(row['frequencies']) or len(row['solver_reports']) != len(row['frequencies']):
        raise Blocked('Malformed branch/solver evidence cardinality')
    for report, frequency in zip(row['solver_reports'], row['frequencies']):
        a, c = report['A'], report['C']
        if a.get('status') != 'PASS' or c.get('status') != 'PASS' or not a.get('converged') or not c.get('converged'):
            raise Blocked('An unconverged root was mislabeled PASS')
        if not a.get('receipt_key') or not c.get('receipt_key'):
            raise Blocked('Missing raw solver receipts')
        if abs(z(a['omega']) - z(c['omega'])) > tolerance or abs(z(a['omega']) - z(frequency)) > tolerance:
            raise Blocked('Reported root disagrees with its solver evidence')
        if z(frequency).imag >= 0 or not math.isfinite(abs(z(frequency))):
            raise Blocked('Invalid damped QNM frequency')
    approaches = row.get('approaches', [])
    if len(approaches) != 2 or {r.get('name') for r in approaches} != {'forward', 'reverse'}:
        raise Blocked('Two independently continued approaches are missing')
    hashes = {r['transport']['path_sha256'] for r in approaches}
    if len(hashes) != 2 or not all(r['transport'].get('adaptive_complete') for r in approaches):
        raise Blocked('Identical or incomplete approaches cannot certify identity')


def accepted_evidence(store, config, plan):
    store.unchanged()
    if store.identity['authority'] != 'COMMITTED_SOURCE_RUN' or store.identity['dirty_checkout']:
        raise Blocked('Development evidence can never freeze claims, set readiness, or create a release')
    if config.get('quasiresonant_scope') != 'EXCLUDED' or config.get('rigorous_claim_scope') != 'FINITE_DISCRETIZATION_ONLY':
        raise Blocked('This producer does not implement a physical hyperboloidal Solver F or a continuum certificate')
    artifacts = {name: store.read(name, accepted) for name, (_, accepted) in ARTIFACTS.items()}
    atlas = artifacts['atlas.json']
    if atlas.get('limited_execution') or atlas.get('inherited_legacy_labels') is not False:
        raise Blocked('Partial or inherited-label atlas')
    if atlas.get('scope', {}).get('kind') != 'DECLARED_FINITE_PARAMETER_ATLAS':
        raise Blocked('Unsupported atlas scope')
    expected_points = {sid(s) + '@' + pid(p) for s in config['sectors'] for p in atlas['scope']['points']}
    if expected_points != set(atlas['tasks']):
        raise Blocked('Atlas does not cover its declared sector/point product')
    for row in atlas['tasks'].values():
        check_physical_row(row, config['frequency_tolerance'])
    actual_minimum = min((p for row in atlas['tasks'].values() for p in pairs_at(row)), key=lambda p: p['gap'])
    recorded = atlas['minimum_candidate']
    if (recorded['point'], recorded['sector'], recorded['overtones']) != (actual_minimum['point'], actual_minimum['sector'], actual_minimum['overtones']):
        raise Blocked('Minimum candidate is not the fresh finite-atlas minimum')
    if abs(recorded['gap'] - actual_minimum['gap']) > 1e-12:
        raise Blocked('Minimum gap does not match the root data')
    convergence = artifacts['convergence.json']
    if convergence.get('axes') != plan['axes']:
        raise Blocked('Convergence ladders differ from the merged plan')
    for study in convergence['tasks'].values():
        require_coverage(study)
        for axis, values in plan['axes'].items():
            result = study['tasks'][axis]
            if result['requested_values'] != values or set(result['rungs']) != set(map(str, values)):
                raise Blocked('A convergence rung was skipped/substituted')
            if not result.get('all_rungs_executed_and_passed') or any(r.get('status') != 'PASS' or r.get('fixed_label_identity_check', {}).get('status') != 'PASS' for r in result['rungs'].values()):
                raise Blocked('A skipped or failed rung was reported converged')
    boundary = artifacts['boundary.json']
    if boundary.get('global_minimum_proved') is not False:
        raise Blocked('Local refinement does not establish a global theorem')
    polls = boundary.get('history', [])
    if len(polls) < 2 or any(not p.get('poll_complete') or p.get('significant_improvement') for p in polls[-2:]):
        raise Blocked('Missing two complete stable terminal boundary polls')
    if any(max(p['radius']) > plan['adaptive_boundary']['parameter_tolerance'] for p in polls[-2:]):
        raise Blocked('Boundary radii did not reach the configured parameter tolerance')
    minimum = boundary['minimum_candidate']
    robust = artifacts['robustness.json']
    expected_modes = {sid({**minimum['sector'], 'ell': minimum['sector']['ell'] + offset}) + '/N' + str(n)
                      for offset in plan['higher_modes']['ell_offsets'] for n in plan['higher_modes']['overtones']}
    if set(robust['tasks']) != expected_modes:
        raise Blocked('Higher-mode coverage is incomplete')
    for row in robust['tasks'].values():
        check_physical_row(row, config['frequency_tolerance'])
        require_coverage(row['resolution_study'])
        if set(row['resolution_study']['tasks']) != set(config['robustness_axes']):
            raise Blocked('High-mode resolution study is incomplete')
        for axis, values in config['robustness_axes'].items():
            result = row['resolution_study']['tasks'][axis]
            if result['requested_values'] != values or set(result['rungs']) != set(map(str,values)):
                raise Blocked('High-mode resolution rung skipped or substituted')
    lean = artifacts['lean.json']
    if not lean.get('compilation_checked') or lean['tasks']['lake-build'].get('returncode') != 0:
        raise Blocked('Lean BUILD_PASS requires a real successful compiler process')
    if not artifacts['bibliography.json']['tasks']['complete-input-inventory'].get('complete_submission_inventory_attested'):
        raise Blocked('The bibliography audit covers only a subset of submission inputs')
    return artifacts


def prepare_claims(store, config, plan):
    evidence = accepted_evidence(store, config, plan)
    source = read_json(store.root / 'config/science/claims.json')
    if not _global_claims_safe(source):
        raise Blocked('Global no-EP claim must remain inactive and NOT_ESTABLISHED')
    proposed = copy.deepcopy(source)
    proposed['science_ready'] = True
    proposed['quasiresonant_scope'] = 'EXCLUDED'
    proposed['campaign_v2'] = {
        'scope': evidence['atlas.json']['scope'], 'claim_level': 'NUMERICAL_NOT_CONTINUUM_RIGOROUS',
        'local_minimum_only': True, 'global_no_ep': 'NOT_ESTABLISHED',
        'preflight_source_commit': store.identity['source_commit'],
        'requires_fresh_reproduction_after_registry_commit': True,
    }
    atomic_json(store.output / 'claims.proposed.json', proposed)
    return {'status': 'CLAIM_PROPOSAL_READY', 'file': 'claims.proposed.json',
            'source_registry_modified': False,
            'next_step': 'Review and commit this registry change, then rerun ALL science from that new clean commit. Old receipts cannot validate the new commit.'}


def publication_figures(store, evidence):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    directory = store.output / 'publication-figures'
    directory.mkdir(exist_ok=True)
    rows = sorted(evidence['atlas.json']['tasks'].values(), key=lambda r: (sid(r['sector']), pid(r['point'])))
    values = [min(pairs_at(row), key=lambda p: p['gap']) for row in rows]
    figure, axis = plt.subplots(figsize=(10, 5.5))
    axis.errorbar(range(len(values)), [r['gap'] for r in values],
                  yerr=[r['empirical_error'] for r in values], fmt='o', capsize=2)
    axis.set_xlabel('Declared finite-atlas sample index')
    axis.set_ylabel('Smallest within-sector frequency separation')
    axis.set_title('Validated finite-atlas samples — not a global no-EP theorem')
    figure.tight_layout()
    outputs = []
    for suffix in ('png', 'svg'):
        output = directory / ('finite-atlas-gaps.' + suffix)
        figure.savefig(output, dpi=180)
        outputs.append(output)
    plt.close(figure)
    atomic_json(directory / 'plotted-data.json', {'source_commit': store.identity['source_commit'],
                                                 'samples': values, 'uncertainty_type': 'EMPIRICAL_NOT_INTERVAL'})
    outputs.append(directory / 'plotted-data.json')
    return {p.relative_to(store.output).as_posix(): file_hash(p) for p in outputs}


def release(store, config, plan, tag=None):
    evidence = accepted_evidence(store, config, plan)
    registry = read_json(store.root / 'config/science/claims.json')
    if registry.get('science_ready') is not True or not registry.get('campaign_v2') or not _global_claims_safe(registry):
        raise Blocked('Commit the reviewed campaign_v2 claim registry, then reproduce again before releasing')
    if tag and not re.fullmatch(r'numerical-v[0-9A-Za-z][0-9A-Za-z._-]*', tag):
        raise Blocked('Use an explicit numerical-v... tag; no forced or ambiguous tag updates')
    if tag and git(store.root, 'tag', '--list', tag):
        raise Blocked('Tag already exists; never overwrite a scientific release')
    figures = publication_figures(store, evidence)
    record = {'science_ready': True, 'source_commit': store.identity['source_commit'],
              'source_tree_sha256': store.identity['source_tree_sha256'],
              'run_id': store.run['run_id'], 'quasiresonant_scope': 'EXCLUDED',
              'claim_level': 'NUMERICAL', 'global_no_ep': 'NOT_ESTABLISHED',
              'continuum_certificate': False, 'local_minimum_only': True,
              'artifact_hashes': {name: file_hash(store.output / name) for name in ARTIFACTS},
              'figure_hashes': figures, 'registry_hash': file_hash(store.root / 'config/science/claims.json'),
              'local_tag_requested': tag}
    atomic_json(store.output / 'release.json', record)
    atomic_json(store.output / 'claims.frozen.json', registry)
    if tag:
        store.unchanged()
        git(store.root, 'tag', '-a', tag, store.identity['source_commit'], '-m',
            'Numerical-only MP5D release. Evidence SHA256: ' + file_hash(store.output / 'release.json'))
        atomic_json(store.output / 'tag-receipt.json', {'local_tag_created': tag,
                    'remote_tag_pushed': False, 'release_sha256': file_hash(store.output / 'release.json'),
                    'source_commit': store.identity['source_commit']})
    return record

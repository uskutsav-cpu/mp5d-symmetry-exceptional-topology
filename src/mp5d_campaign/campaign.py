"""Fresh atlas, exact seven-axis study, local refinement, and higher modes."""
from __future__ import annotations

import itertools
import math
from dataclasses import asdict, replace

import mpmath as mp

from mp5d_science.convergence import Resolution

from .core import ARTIFACTS, Blocked, file_hash, payload, read_json
from .numerics import Physics, diameter, distinct, pid, sid, z


def build_points(config, historical):
    grid = config['atlas_grid']
    result = {pid(p): list(map(float, p)) for p in itertools.product(*(grid[k] for k in ('s', 'delta', 'mu')))}
    # Historical frequencies, branch names, and gap estimates are never read.
    for row in historical:
        point = row['point']
        result[pid(point)] = list(map(float, point))
    for point in result.values():
        if any(not lo <= v <= hi for v, (lo, hi) in zip(point, config['bounds'])):
            raise Blocked('Atlas point outside the declared formal domain')
    return [result[k] for k in sorted(result)]


def pairs_at(record):
    frequencies = record['frequencies']
    errors = record['empirical_errors']
    for i, j in itertools.combinations(range(len(frequencies)), 2):
        yield {'gap': abs(z(frequencies[i]) - z(frequencies[j])),
               'empirical_error': errors[i] + errors[j], 'indices': [i, j],
               'point': record['point'], 'sector': record['sector'],
               'overtones': [record['overtones'][i], record['overtones'][j]],
               'frequencies': [frequencies[i], frequencies[j]],
               'branch_ids': [sid(record['sector']) + '/N' + str(record['overtones'][k]) for k in (i, j)]}


class Campaign:
    def __init__(self, store, config, plan):
        self.store, self.config, self.plan = store, config, plan
        self.physics = Physics(store, config)
        self.resolution = Resolution(**config['resolution'])
        self.atlas_resolution = Resolution(**config.get('atlas_resolution', config['resolution']))
        self.historical = read_json(store.root / config['dangerous_points_file'])['points']
        self.points = build_points(config, self.historical)
        self.sectors = config['sectors']
        self.formal_n = plan['formal_domain_proposal']['overtones']
        expected_bounds = [plan['formal_domain_proposal'][k] for k in ('s', 'delta', 'mu')]
        if config['bounds'] != expected_bounds:
            raise Blocked('Campaign bounds do not match the merged formal-domain proposal')
        if not self.sectors or len({sid(s) for s in self.sectors}) != len(self.sectors):
            raise Blocked('Sector manifest must be nonempty and unique')
        for s in self.sectors:
            from .numerics import params
            params([0, 0, 0], s)  # validates admissible ell and azimuthal integers

    def reference(self, overtone):
        entries = read_json(self.store.root / self.config['static_references'])['records']
        rows = [r for r in entries if r['ell'] == 2 and r['N'] == overtone]
        if not rows:
            raise Blocked(f'No trusted, source-audited static N={overtone} anchor; inversion alone cannot label it')
        row = rows[0]
        if row.get('label_provenance') != 'PUBLISHED_TABLE':
            audit = row.get('label_audit', {})
            evidence = self.store.root / audit.get('evidence_path', '__missing__')
            if audit.get('status') != 'VERIFIED' or not evidence.is_file() or file_hash(evidence) != audit.get('sha256'):
                raise Blocked('Additional overtone anchor lacks a verifiable mode-label audit')
        with mp.workdps(100):
            scale = 2 * mp.pi if row['normalization'] == 'omega_tilde=2*pi*omega; rh=1' else mp.mpf(1)
            if row['normalization'] not in ('omega_tilde=2*pi*omega; rh=1', 'omega; rh=1'):
                raise Blocked('Unknown frequency normalization')
            root = [str(mp.mpf(row['frequency'][0]) / scale), str(mp.mpf(row['frequency'][1]) / scale)]
        return root, row

    def anchor(self, sector, overtones):
        request = {'algorithm': 'fresh-published-anchor-with-lambda-transport/v1',
                   'sector': sector, 'overtones': overtones, 'anchor_settings': self.config['anchor_resolution']}

        def compute():
            references = [self.reference(n) for n in overtones]
            seeds = [r[0] for r in references]
            source_sector = {'m1': 2, 'm2': 0, 'ell': 2}
            base = Resolution(**self.config['anchor_resolution'])
            rows = []
            for depth, radial in zip(self.config['anchor_depths'], self.config['anchor_radial_n']):
                res = replace(base, cf_depth=depth, radial_n=radial)
                row = self.physics.group([0., 0., 0.], source_sector, overtones, seeds, res,
                                         self.config['anchor_tolerance'])
                if any(abs(z(a) - z(b)) > self.config['anchor_tolerance'] for a, b in zip(row['frequencies'], seeds)):
                    raise Blocked('Static root moved outside its published anchor neighborhood')
                rows.append(row)
            drift = diameter([r['frequency_decimals'] for r in rows])
            if max(drift) > self.config['anchor_tolerance']:
                raise Blocked('Trusted static anchor failed the anchor convergence preflight')
            transport = self.physics.lambda_transport(overtones, rows[-1]['frequencies'], sector['ell'], base)
            if transport['status'] != 'PASS':
                raise Blocked(transport.get('reason', 'Static angular homotopy failed'))
            endpoint = self.physics.group([0., 0., 0.], sector, overtones,
                                          transport['frequencies'], base, self.config['anchor_tolerance'])
            if any(abs(z(a) - z(b)) > self.config['anchor_tolerance'] for a, b in
                   zip(endpoint['frequencies'], transport['frequencies'])):
                raise Blocked('Physical integer-ell endpoint lost its transported identity')
            return {**endpoint, 'sources': [r[1] for r in references], 'anchor_preflight': rows,
                    'anchor_empirical_drift': drift, 'angular_transport': transport,
                    'inherited_legacy_labels': False, 'inversion_index_is_identity_proof': False}
        return self.store.task(request, compute)

    def atlas(self, point_limit=None):
        if len(self.historical) != self.config['required_historical_points']:
            raise Blocked('Atlas cannot omit any of the ten historical dangerous coordinates')
        if point_limit is not None and not self.store.development:
            raise Blocked('A point budget is development-only; it can never validate an atlas')
        expected, tasks, anchors = [], {}, {}
        executed = 0
        for sector in self.sectors:
            key = sid(sector)
            anchor = self.anchor(sector, self.formal_n)
            anchors[key] = anchor
            for point in self.points:
                task_id = key + '@' + pid(point)
                expected.append(task_id)
                if anchor['status'] != 'PASS':
                    tasks[task_id] = {'status': 'BLOCKED', 'reason': anchor.get('reason', 'Unvalidated anchor')}
                elif point_limit is not None and executed >= point_limit:
                    tasks[task_id] = {'status': 'NOT_RUN', 'reason': 'Explicit development point budget'}
                else:
                    executed += 1
                    tasks[task_id] = self.store.task(
                        {'algorithm': 'fresh-atlas-point/v1', 'anchor': anchor, 'target': point,
                         'resolution': asdict(self.atlas_resolution)},
                        lambda a=anchor, p=point: self.physics.approaches(a, p, self.atlas_resolution))
        candidates = [p for row in tasks.values() if row['status'] == 'PASS' for p in pairs_at(row)]
        minimum = min(candidates, key=lambda p: p['gap']) if candidates else None
        result = payload('VALIDATED', expected, tasks, anchors=anchors, minimum_candidate=minimum,
            scope={'kind': 'DECLARED_FINITE_PARAMETER_ATLAS', 'bounds': self.config['bounds'],
                   'sectors': self.sectors, 'overtones': self.formal_n, 'points': self.points, 'quasiresonant_exclusion': self.config.get('quasiresonant_exclusion')},
            inherited_legacy_labels=False, global_domain_coverage_proved=False,
            global_no_ep='NOT_ESTABLISHED', quasiresonant_scope='EXCLUDED',
            limited_execution=point_limit is not None)
        if point_limit is not None:
            result['status'] = 'DEVELOPMENT_PARTIAL'
        self.store.write('atlas.json', ARTIFACTS['atlas.json'][0], result)
        return result

    def candidates(self):
        if len(self.historical) != self.config['required_historical_points']:
            raise Blocked('Historical dangerous-point inventory is incomplete; expected all ten coordinates')
        atlas = self.store.read('atlas.json', 'VALIDATED')
        tasks = {}
        for row in self.historical:
            candidates = [value for value in atlas['tasks'].values()
                          if value['status'] == 'PASS' and pid(value['point']) == pid(row['point'])]
            tasks[row['id']] = {'status': 'PASS' if candidates else 'FAILED',
                               'point': row['point'], 'fresh_spectra': candidates,
                               'legacy_frequencies_and_branch_labels_used': False,
                               'legacy_close_pair_reaffirmed': False}
        data = payload('POINTWISE_CHECKS_PASS', [r['id'] for r in self.historical], tasks)
        self.store.write('candidates.json', ARTIFACTS['candidates.json'][0], data, ['atlas.json'])
        return data

    def dangerous(self, atlas):
        rows = list(atlas['tasks'].values())
        chosen = {sid(r['sector']) + '@' + pid(r['point']): r for r in rows
                  if any(pid(r['point']) == pid(h['point']) for h in self.historical)}
        # Also audit the smallest gaps, highest empirical errors, and every sampled
        # boundary point. No stale candidate is silently treated as the minimum.
        for row in sorted(rows, key=lambda r: min(p['gap'] for p in pairs_at(r)))[:self.config['small_gap_points']]:
            chosen[sid(row['sector']) + '@' + pid(row['point'])] = row
        for row in sorted(rows, key=lambda r: max(r['empirical_errors']), reverse=True)[:self.config['uncertain_points']]:
            chosen[sid(row['sector']) + '@' + pid(row['point'])] = row
        for row in rows:
            if any(abs(v - lo) < 1e-12 or abs(v - hi) < 1e-12 for v, (lo, hi) in zip(row['point'], self.config['bounds'])):
                chosen[sid(row['sector']) + '@' + pid(row['point'])] = row
        return [chosen[k] for k in sorted(chosen)]

    def ladders(self, row, anchor, axes=None, base=None):
        axes = self.plan["axes"] if axes is None else axes
        base = self.resolution if base is None else base
        # Restrict a full atlas anchor to the exact requested pair after boundary
        # refinement; never compare different branch sets across axes.
        indices = [anchor['overtones'].index(n) for n in row['overtones']]
        anchor = {**anchor, 'overtones': list(row['overtones']),
                  'frequencies': [anchor['frequencies'][i] for i in indices]}
        local_anchor = {'point': row['point'], 'sector': row['sector'],
                        'overtones': row['overtones'], 'frequencies': row['frequencies'],
                        'status': 'PASS'}
        tasks = {}
        for axis, values in axes.items():
            rungs = {}
            for value in values:
                resolution = replace(base, **{axis: value})
                request = {'algorithm': 'seven-axis-rung/v1', 'axis': axis, 'value': value,
                           'point': row['point'], 'sector': row['sector'], 'overtones': row['overtones'],
                           'seeds': row['frequencies'], 'resolution': asdict(resolution),
                           'anchor': local_anchor if axis == 'continuation_step' else None}
                def evaluate_rung(res=resolution, ax=axis):
                    result = (self.physics.approaches(local_anchor, row['point'], res)
                              if ax == 'continuation_step' else
                              self.physics.group(row['point'], row['sector'], row['overtones'], row['frequencies'], res))
                    baseline = list(map(z, row['frequencies']))
                    displacements = []
                    for index, frequency in enumerate(map(z, result['frequencies'])):
                        movement = abs(frequency - baseline[index])
                        own_gap = min(abs(baseline[index] - other) for j, other in enumerate(baseline) if j != index)
                        competing = min(abs(frequency - other) for j, other in enumerate(baseline) if j != index)
                        if movement >= .45 * own_gap or competing <= movement:
                            raise Blocked('Resolution rung left its fixed-label identity neighborhood')
                        displacements.append(movement)
                    result['fixed_label_identity_check'] = {'status': 'PASS', 'displacements': displacements}
                    if ax == 'continuation_step':
                        result['continuation_scope'] = 'LOCAL_NONTRIVIAL_CLOSED_PATHS_FROM_VALIDATED_POINT'
                    return result
                rungs[str(value)] = self.store.task(request, evaluate_rung)
            complete = all(r['status'] == 'PASS' for r in rungs.values())
            spread = None
            stable = False
            if complete:
                spread = diameter([r['frequency_decimals'] for r in rungs.values()])
                # Both discretizations matter: track C's own tail as well, not
                # only the A root that does not change during a radial sweep.
                c_spread = diameter([[list(map(str, report['C']['omega'])) for report in r['solver_reports']]
                                     for r in rungs.values()])
                spread = [max(a, c) for a, c in zip(spread, c_spread)]
                stable = max(spread) <= self.config['frequency_tolerance']
                try:
                    last = list(rungs.values())[-1]
                    distinct(last['frequencies'], [max(e, d) for e, d in zip(last['empirical_errors'], spread)])
                except Blocked:
                    stable = False
            tasks[axis] = {'status': 'PASS' if complete and stable else 'FAILED',
                           'requested_values': values, 'rungs': rungs, 'tail_pairwise_diameter': spread,
                           'all_rungs_executed_and_passed': complete,
                           'precision_applies_to': 'A-direct; C remains float64' if axis == 'precision_dps' else 'see effective_resolution',
                           'empirical_stability_only': True}
        return payload('PASS', list(axes), tasks,
                       point=row['point'], sector=row['sector'], branch_ids=row.get('branch_ids'),
                       continuum_error_bound=False)

    def convergence(self):
        atlas = self.store.read('atlas.json', 'VALIDATED')
        tasks = {}
        rows = self.dangerous(atlas)
        for row in rows:
            key = sid(row['sector']) + '@' + pid(row['point'])
            tasks[key] = self.ladders(row, atlas['anchors'][sid(row['sector'])])
        data = payload('PASS', [sid(r['sector']) + '@' + pid(r['point']) for r in rows], tasks,
                       axes=self.plan['axes'], selection='historical coordinates + smallest gaps + highest errors + all boundary samples')
        self.store.write('convergence.json', ARTIFACTS['convergence.json'][0], data, ['atlas.json'])
        return data

    def boundary(self):
        atlas = self.store.read('atlas.json', 'VALIDATED')
        self.store.read('convergence.json', 'PASS')
        candidate = atlas['minimum_candidate']
        if not candidate:
            raise Blocked('No validated finite-atlas minimum')
        spec = self.plan['adaptive_boundary']
        anchor = {**candidate, 'status': 'PASS'}
        observations = {}
        radius = list(map(float, spec['initial_radius']))
        max_evaluations = int(spec['max_evaluations'])
        stable_polls = 0
        best = None
        history = []
        exhausted = False

        def evaluate(point):
            key = pid(point)
            if key not in observations:
                if len(observations) >= max_evaluations:
                    raise Blocked('Boundary evaluation budget exhausted')
                record = self.store.task({'algorithm': 'boundary-two-approach-observation/v1',
                    'anchor': anchor, 'target': list(point), 'resolution': asdict(self.resolution)},
                    lambda: self.physics.approaches(anchor, point, self.resolution))
                if record['status'] != 'PASS':
                    observations[key] = record
                    raise Blocked('Boundary point failed branch/solver checks: ' + record.get('reason', key))
                pair_record = list(pairs_at(record))[0]
                observations[key] = {**record, **pair_record}
            result = observations[key]
            if result['status'] != 'PASS':
                raise Blocked('Cached boundary point failed')
            if result['empirical_error'] > spec['gap_tolerance'] / 10:
                raise Blocked('Numerical discrepancy is too large for the requested boundary gap tolerance')
            return result

        try:
            best = evaluate(candidate['point'])
            while len(history) < self.config['boundary_max_polls']:
                center = best['point']
                stencil = {pid(p): p for offsets in itertools.product((-1, 0, 1), repeat=3)
                           for p in [[min(hi, max(lo, x + o * h)) for x, o, h, (lo, hi) in
                                      zip(center, offsets, radius, self.config['bounds'])]]}
                evaluated = [evaluate(p) for p in stencil.values()]
                trial = min(evaluated, key=lambda x: x['gap'])
                improved = (trial['gap'] + trial['empirical_error'] + spec['gap_tolerance'] <
                            best['gap'] - best['empirical_error'])
                history.append({'center': center, 'radius': radius.copy(), 'poll_complete': True,
                                'poll_points': sorted(stencil), 'significant_improvement': improved,
                                'best_gap': trial['gap']})
                if improved:
                    # Do not shrink while descending; otherwise shrinking alone
                    # can manufacture a spurious numerical stationary point.
                    best = trial
                    stable_polls = 0
                else:
                    if max(radius) <= spec['parameter_tolerance']:
                        stable_polls += 1
                    if stable_polls >= 2:
                        break
                    radius = [h / 2 for h in radius]
            if stable_polls < 2:
                raise Blocked('No two complete stable terminal polls')
        except Blocked as exc:
            exhausted = True
            failure = str(exc)
        # A moved minimum gets a NEW full seven-axis study; initial evidence is
        # not automatically inherited by a refined point.
        final_study = self.ladders(best, atlas['anchors'][sid(best['sector'])]) if best and not exhausted else None
        tasks = {'complete_local_polls': {'status': 'PASS' if not exhausted else 'FAILED',
                    'reason': failure if exhausted else None, 'stable_terminal_polls': stable_polls},
                 'refined_minimum_seven_axes': {'status': 'PASS' if final_study and final_study['status'] == 'PASS' else 'FAILED',
                                               'study': final_study}}
        data = payload('LOCALLY_STABLE_NUMERICAL_MINIMUM', list(tasks), tasks,
                       minimum_candidate=best, observations=observations, history=history,
                       parameter_tolerance=spec['parameter_tolerance'], gap_tolerance=spec['gap_tolerance'],
                       independent_directions=['forward', 'reverse'], global_minimum_proved=False,
                       scope='LOCAL_NUMERICAL_STENCIL_STABILITY_ONLY')
        self.store.write('boundary.json', ARTIFACTS['boundary.json'][0], data, ['atlas.json', 'convergence.json'])
        return data

    def robustness(self):
        boundary = self.store.read('boundary.json', 'LOCALLY_STABLE_NUMERICAL_MINIMUM')
        minimum = boundary['minimum_candidate']
        expected, tasks = [], {}
        resolution = Resolution(**self.config['robustness_resolution'])
        for offset in self.plan['higher_modes']['ell_offsets']:
            sector = {**minimum['sector'], 'ell': minimum['sector']['ell'] + offset}
            # N=4 and N=5 MUST be transported jointly, not in separate runs that
            # could converge onto the same root and each incorrectly pass.
            overtones = [0] + list(self.plan['higher_modes']['overtones'])
            keys = [sid(sector) + '/N' + str(n) for n in self.plan['higher_modes']['overtones']]
            expected.extend(keys)
            anchor = self.anchor(sector, overtones)
            if anchor['status'] != 'PASS':
                for key in keys:
                    tasks[key] = {'status': 'BLOCKED', 'reason': anchor.get('reason', 'Joint higher-mode anchor failed'),
                                  'anchor': anchor, 'required_joint_labels': overtones}
                continue
            result = self.store.task({'algorithm': 'joint-higher-mode-two-approach-probe/v1', 'anchor': anchor,
                'point': minimum['point'], 'resolution': asdict(resolution)},
                lambda a=anchor: self.physics.approaches(a, minimum['point'], resolution))
            if result['status'] == 'PASS':
                study = self.ladders(result, anchor, self.config['robustness_axes'], resolution)
                result = {**result, 'status': 'PASS' if study['status'] == 'PASS' else 'FAILED', 'resolution_study': study}
            for key in keys:
                tasks[key] = {**result, 'required_joint_labels': overtones}
        data = payload('VALIDATED_OUT_OF_DOMAIN_PROBE', expected, tasks,
                       point=minimum['point'], scope='OUT_OF_DOMAIN_ROBUSTNESS_PROBE',
                       probe_resolution=asdict(resolution), probe_axes=self.config['robustness_axes'],
                       substitutes_for_core_seven_axis_study=False,
                       expands_formal_overtone_domain=False, global_no_ep='NOT_ESTABLISHED')
        self.store.write('robustness.json', ARTIFACTS['robustness.json'][0], data, ['boundary.json'])
        return data

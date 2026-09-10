"""Actual MP5D solvers, joint branch transport, and independent approaches.

An inversion index is an algorithm setting, NEVER evidence of mode identity.
Errors below are empirical numerical discrepancies, not interval enclosures.
"""
from __future__ import annotations

import itertools
import cmath
from decimal import Decimal
import math
from dataclasses import asdict, replace, is_dataclass
from pathlib import Path

import mpmath as mp
import numpy as np

from mp5d_science.convergence import Resolution
from mp5d_science.physical import solve_a, solve_c
from mp5d_science.radial_polynomial import RadialParameters, coefficient_formula, mul, power
from mp5d_science.tracking import Mode, TrackingPolicy, continue_branches

from .core import Blocked, digest


def plain(value):
    """Lossless finite diagnostic JSON; nonfinite diagnostics stay explicit strings."""
    if is_dataclass(value):
        return plain(asdict(value))
    if isinstance(value, dict):
        return {str(k): plain(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [plain(v) for v in value]
    if isinstance(value, complex):
        return pair(value)
    if hasattr(value, 'tolist'):
        return plain(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if value is None or isinstance(value, (str,int,float,bool)):
        return value
    return repr(value)


def z(value):
    return complex(*map(float, value)) if isinstance(value, (tuple, list)) else complex(value)


def pair(value):
    value = complex(value)
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise Blocked('Nonfinite frequency')
    return [value.real, value.imag]


def sid(sector):
    return 'm%d_%d_ell%d' % (sector['m1'], sector['m2'], sector['ell'])


def pid(point):
    return ','.join(format(float(x), '.12g') for x in point)


def params(point, sector):
    s, delta, mass = (Decimal(str(value)) for value in point)
    if not all(value.is_finite() for value in (s, delta, mass)):
        raise Blocked('Nonfinite physical parameter')
    return RadialParameters(a=str(s + delta), b=str(s - delta), mu=str(mass),
                            m1=sector['m1'], m2=sector['m2'], ell=sector['ell'], M='1')


def distinct(frequencies, errors, floor=1e-8):
    if len(frequencies) != len(errors) or not frequencies:
        raise Blocked('Malformed root/error vector')
    if any(not math.isfinite(abs(z(w))) for w in frequencies) or any(not math.isfinite(e) or e < 0 for e in errors):
        raise Blocked('Nonfinite root or invalid numerical discrepancy')
    for i, j in itertools.combinations(range(len(frequencies)), 2):
        gap = abs(z(frequencies[i]) - z(frequencies[j]))
        if gap <= max(floor, 10 * (errors[i] + errors[j])):
            raise Blocked(f'Unresolved or duplicate branch roots {i}/{j}: gap={gap:.9g}')


def route(start, end, order, step):
    if step <= 0 or not math.isfinite(step):
        raise ValueError('Continuation step must be finite and positive')
    here = list(map(float, start))
    points = [here.copy()]
    for axis in order:
        target = float(end[axis])
        n = max(1, math.ceil(abs(target - here[axis]) / step))
        origin = here[axis]
        for k in range(1, n + 1):
            candidate = here.copy()
            candidate[axis] = target if k == n else origin + (target - origin) * k / n
            if candidate != points[-1]:
                points.append(candidate)
        here[axis] = target
    return points


def join_routes(*routes):
    result = []
    for piece in routes:
        for point in piece:
            if not result or point != result[-1]:
                result.append(point)
    return result


def approach_routes(start, end, bounds, step, distance=0.01):
    """Different terminal approach directions, not merely different early detours.

    Both paths start at the same validated anchor and approach the target along
    different coordinate axes. The penultimate waypoint is near the TARGET, so
    paths cannot merge into the same final continuation segment unnoticed.
    """
    if not math.isfinite(distance) or distance <= 0:
        raise Blocked('Approach detour must be finite and positive')
    paths = []
    for axis, order in [(0, (0, 1, 2)), (1, (2, 1, 0))]:
        waypoint = list(map(float, end))
        lo, hi = bounds[axis]
        room_up, room_down = hi - waypoint[axis], waypoint[axis] - lo
        if max(room_up, room_down) <= 1e-12:
            raise Blocked('No independent terminal approach available in the configured domain')
        sign = 1 if room_up >= room_down else -1
        waypoint[axis] += sign * min(distance, max(room_up, room_down))
        if all(abs(a-b) < 1e-12 for a,b in zip(waypoint,start)):
            waypoint[axis] = float(end[axis]) + sign * min(2*distance, max(room_up,room_down))
        path = join_routes(route(start, waypoint, order, step), route(waypoint, end, order, step))
        if len(path) < 3:
            raise Blocked('Trivial continuation route cannot test a step-size ladder')
        paths.append(path)
    if digest(paths[0]) == digest(paths[1]):
        raise Blocked('Independent approach paths are identical')
    final_axes = [{i for i,(a,b) in enumerate(zip(path[-2],path[-1])) if abs(a-b)>1e-12} for path in paths]
    if final_axes != [{0},{1}]:
        raise Blocked('Terminal approach directions are not independently resolved')
    return paths


def diameter(rows):
    """Full pairwise diameter of final three resolutions, using decimal roots.

    Adjacent-rung differences alone undercount monotone drift. Decimal strings
    avoid erasing arbitrary precision improvements through float64 conversion.
    """
    if len(rows) < 3:
        raise Blocked('Need three terminal rungs')
    count = len(rows[-1])
    if any(len(row) != count for row in rows):
        raise Blocked('Branch count changed between rungs')
    with mp.workdps(110):
        values = [[mp.mpc(str(v[0]), str(v[1])) for v in row] for row in rows[-3:]]
        if any(not mp.isfinite(w) for row in values for w in row):
            raise Blocked('Nonfinite arbitrary-precision root')
        return [float(max(abs(a[k] - b[k]) for a, b in itertools.combinations(values, 2)))
                for k in range(count)]


class Physics:
    def __init__(self, store, config):
        self.store, self.config = store, config
        if not store.development:
            import mp5d_science.physical as imported
            expected = (store.root / 'src/mp5d_science/physical.py').resolve()
            if Path(imported.__file__).resolve() != expected:
                raise Blocked('Physics was imported from a different source checkout')

    def solve(self, kind, point, sector, overtone, seed, resolution):
        r = asdict(resolution)
        effective = ({k: r[k] for k in ('cf_depth', 'angular_n', 'precision_dps')}
                     if kind == 'A' else
                     {**{k: r[k] for k in ('radial_n', 'angular_n', 'contour_length', 'scaling_angle_deg')},
                      'precision_dps': 16})
        request = {'algorithm': 'MP5D-A-direct' if kind == 'A' else 'MP5D-C-bordered',
                   'point': list(point), 'sector': sector, 'inversion': overtone,
                   'seed': list(seed) if isinstance(seed, (tuple, list)) else pair(seed),
                   'effective_resolution': effective}

        def compute():
            p = params(point, sector)
            if kind == 'A':
                root = solve_a(p, tuple(request['seed']), overtone=overtone, depth=r['cf_depth'],
                               angular_n=r['angular_n'], dps=r['precision_dps'])
            else:
                root = solve_c(p, z(seed), radial_n=r['radial_n'], angular_n=r['angular_n'],
                               length=r['contour_length'], angle_deg=r['scaling_angle_deg'])
            data = root.to_dict()
            data['status'] = 'PASS' if root.converged and root.residual <= self.config['residual_tolerance'] else 'FAILED'
            data['effective_resolution'] = effective
            data['mode_label_certified_by_inversion'] = False
            return data
        return self.store.task(request, compute)

    def group(self, point, sector, overtones, seeds, resolution, tolerance=None):
        if len(overtones) != len(seeds) or not seeds:
            raise Blocked('Nonempty one-to-one overtone/seed list required')
        tolerance = tolerance or self.config['frequency_tolerance']
        modes, errors, decimals, reports = [], [], [], []
        for overtone, seed in zip(overtones, seeds):
            # C starts independently from the supplied predictor, not A's result.
            a = self.solve('A', point, sector, overtone, seed, resolution)
            c = self.solve('C', point, sector, overtone, seed, resolution)
            reports.append({'overtone_transport_label': overtone, 'A': a, 'C': c})
            if a['status'] != 'PASS' or c['status'] != 'PASS':
                raise Blocked(f'Solver check failed for transported N={overtone}: A={a.get("reason", a["status"])}, C={c.get("reason", c["status"])}')
            error = abs(z(a['omega']) - z(c['omega']))
            if error > tolerance:
                raise Blocked(f'A/C disagreement for transported N={overtone}: {error:.9g} > {tolerance:.9g}')
            cut = self.config.get('quasiresonant_exclusion')
            if cut:
                frequency = z(a['omega'])
                outgoing_k = cmath.sqrt(frequency ** 2 - float(point[2]) ** 2)
                if -frequency.imag < cut['minimum_damping'] or abs(outgoing_k) < cut['minimum_abs_outgoing_wave_number']:
                    raise Blocked('EXCLUDED_QUASIRESONANT_REGIME: outside the registered numerical claim domain')
            modes.append(a['omega'])
            decimals.append(a.get('omega_decimal', list(map(str, a['omega']))))
            errors.append(error)
        distinct(modes, errors)
        return {'status': 'PASS', 'point': list(point), 'sector': sector,
                'overtones': list(overtones), 'frequencies': modes, 'frequency_decimals': decimals,
                'empirical_errors': errors, 'solver_reports': reports,
                'independent_radial_formulations': 2, 'independent_angular_formulations': 1,
                'uncertainty_is_rigorous_enclosure': False}

    def track(self, sector, overtones, seeds, path, resolution):
        request = {'algorithm': 'joint-adaptive-one-to-one/v1', 'sector': sector,
                   'overtones': list(overtones), 'initial_frequencies': [pair(x) for x in seeds],
                   'requested_path': path, 'resolution': asdict(resolution),
                   'policy': self.config['tracking_policy']}

        def compute():
            def solver(point, predictors):
                values = []
                for overtone, predictor in zip(overtones, predictors):
                    root = self.solve('A', list(map(float, point)), sector, overtone, predictor, resolution)
                    if root['status'] != 'PASS':
                        raise Blocked(root.get('reason', 'A root failed during joint continuation'))
                    values.append(Mode(z(root['omega']), root['residual'], converged=True))
                return values
            policy = TrackingPolicy(**self.config['tracking_policy'])
            result = continue_branches(solver, path,
                                       [Mode(z(x), 0.0, converged=True) for x in seeds], policy=policy)
            history = [{'point': list(map(float, p)), 'frequencies': [pair(m.omega) for m in row]}
                       for p, row in zip(result.parameters, result.modes)]
            if not result.complete or not history:
                return {'status': 'FAILED', 'reason': 'Adaptive continuation did not complete: ' + str(result.failure),
                        'history': history, 'retries': plain(result.retries),
                        'evaluations': int(result.evaluations), 'adaptive_complete': False}
            if not np.allclose(history[-1]['point'], path[-1], rtol=0, atol=1e-12):
                raise Blocked('Continuation endpoint mismatch')
            return {'status': 'PASS', 'history': history,
                    'accepted_times': list(map(float, result.accepted_times)),
                    'evaluations': int(result.evaluations), 'retries': plain(result.retries),
                    'requested_vertices_completed': result.requested_vertices_completed,
                    'frequencies': history[-1]['frequencies'],
                    'path_sha256': digest(path), 'adaptive_complete': result.complete}
        return self.store.task(request, compute)

    def approaches(self, anchor, target, resolution):
        paths = approach_routes(anchor['point'], target, self.config['bounds'],
                                resolution.continuation_step, self.config['approach_detour'])
        groups = []
        for name, path in zip(('forward', 'reverse'), paths):
            transport = self.track(anchor['sector'], anchor['overtones'], anchor['frequencies'], path, resolution)
            if transport['status'] != 'PASS':
                raise Blocked(f'{name} approach: {transport.get("reason", "failed")}')
            endpoint = self.group(target, anchor['sector'], anchor['overtones'],
                                  transport['frequencies'], resolution)
            # Terminal A/C roots must still be within the tracked root's identity radius.
            for previous, new in zip(transport['frequencies'], endpoint['frequencies']):
                if abs(z(previous) - z(new)) > self.config['frequency_tolerance']:
                    raise Blocked('Endpoint polishing left the tracked branch')
            groups.append({'name': name, 'transport': transport, 'endpoint': endpoint})
        drift = [abs(z(a) - z(b)) for a, b in zip(groups[0]['endpoint']['frequencies'],
                                                 groups[1]['endpoint']['frequencies'])]
        if max(drift, default=math.inf) > self.config['frequency_tolerance']:
            raise Blocked('Independent paths disagree on branch identity')
        endpoint = groups[0]['endpoint']
        endpoint['empirical_errors'] = [max(e, d) for e, d in zip(endpoint['empirical_errors'], drift)]
        distinct(endpoint['frequencies'], endpoint['empirical_errors'])
        return {**endpoint, 'approaches': groups, 'approach_discrepancies': drift,
                'independent_approach_directions': ['forward', 'reverse'],
                'branch_ids': [sid(anchor['sector']) + f'/N{n}' for n in anchor['overtones']]}

    def static_lambda_root(self, lam, seed, overtone, resolution):
        """Transport a trusted static overtone in continuous separation constant.

        Fractional Lambda is only a transport homotopy, not an admissible physical
        ell mode. Every integer-ell endpoint is checked with the physical A and C
        solvers. Common polynomial factors match A-direct exactly.
        """
        request = {'algorithm': 'static-separation-constant-homotopy/v1', 'lambda': str(lam),
                   'seed': pair(seed), 'inversion': overtone,
                   'depth': resolution.cf_depth, 'dps': resolution.precision_dps}

        def compute():
            from mp5d_science.physical import recurrence_cf
            with mp.workdps(resolution.precision_dps):
                L = mp.mpf(str(lam))
                w0 = mp.mpc(str(z(seed).real), str(z(seed).imag))
                lift = power([mp.mpf(1), mp.mpf(-1)], 4)

                def residual(w, inversion=overtone):
                    A, B, C, _ = coefficient_formula(a=mp.mpf(0), b=mp.mpf(0), mu=mp.mpf(0), m1=0, m2=0,
                        omega=w, lam=L, p=mp.mpf(1), m=mp.mpf(0), outgoing_k=w, imaginary_unit=mp.j)
                    return recurrence_cf(mul(A, lift), mul(B, lift), mul(C, lift),
                                         depth=resolution.cf_depth, inversion=inversion)
                root = mp.findroot(residual, (w0 * (1 - mp.mpf('.00013')), w0 * (1 + mp.mpf('.00017'))),
                                   tol=mp.power(10, -resolution.precision_dps + 8), maxsteps=50)
                cross = max(abs(residual(root, i)) for i in range(max(4, overtone + 2)))
                valid = root.real > 0 and root.imag < 0 and cross < self.config['residual_tolerance']
                return {'status': 'PASS' if valid else 'FAILED', 'omega': pair(root),
                        'omega_decimal': [str(root.real), str(root.imag)], 'residual': float(cross)}
        return self.store.task(request, compute)

    def lambda_transport(self, overtones, seeds, target_ell, resolution):
        endpoint = target_ell * (target_ell + 2)
        steps = max(1, math.ceil(abs(endpoint - 8) / self.config['lambda_step']))
        path = [[8 + (endpoint - 8) * k / steps] for k in range(steps + 1)]
        if endpoint == 8:
            path = [[8.0], [8.5], [8.0]]

        def compute():
            def solver(point, predictors):
                result = []
                for n, seed in zip(overtones, predictors):
                    root = self.static_lambda_root(point[0], seed, n, resolution)
                    if root['status'] != 'PASS':
                        raise Blocked('Static Lambda transport solver failed')
                    result.append(Mode(z(root['omega']), root['residual'], converged=True))
                return result
            tracked = continue_branches(solver, path, [Mode(z(x), 0.0, converged=True) for x in seeds],
                                         policy=TrackingPolicy(**self.config['tracking_policy']))
            if not tracked.complete:
                return {'status': 'FAILED', 'reason': 'Static Lambda transport incomplete: ' + str(tracked.failure),
                        'retries': plain(tracked.retries), 'accepted_parameters': plain(tracked.parameters),
                        'accepted_modes': plain(tracked.modes)}
            return {'status': 'PASS', 'frequencies': [pair(m.omega) for m in tracked.modes[-1]],
                    'history': [{'lambda': float(p[0]), 'frequencies': [pair(m.omega) for m in row]}
                                for p, row in zip(tracked.parameters, tracked.modes)],
                    'retries': plain(tracked.retries), 'not_a_physical_fractional_ell_claim': True}
        return self.store.task({'algorithm': 'joint-lambda-transport/v1', 'path': path,
                               'overtones': overtones, 'seeds': seeds, 'resolution': asdict(resolution)}, compute)


def rung_resolution(base, axis, value):
    return replace(base, **{axis: value})

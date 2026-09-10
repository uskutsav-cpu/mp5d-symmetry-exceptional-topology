import copy

import pytest

from mp5d_campaign.core import (AXES, Blocked, Store, atomic_json, coverage, digest,
                                payload, read_json, require_coverage, validate_plan)


def plan():
    return {'axes': copy.deepcopy(AXES),
            'higher_modes': {'overtones': [4, 5], 'ell_offsets': [0, 2, 4]},
            'formal_domain_proposal': {'overtones': [0, 1, 2, 3]},
            'adaptive_boundary': {'independent_continuation_directions': ['forward', 'reverse']}}


def source(tmp_path):
    root = tmp_path / 'src'
    root.mkdir()
    (root / 'hello.py').write_text('x = 1\n')
    return root


def test_exact_plan():
    validate_plan(plan())


@pytest.mark.parametrize('axis', AXES)
def test_no_dropped_rungs(axis):
    changed = plan()
    changed['axes'][axis] = changed['axes'][axis][:-1]
    with pytest.raises(Blocked):
        validate_plan(changed)


def test_missing_higher_mode():
    changed = plan()
    changed['higher_modes']['overtones'] = [4]
    with pytest.raises(Blocked):
        validate_plan(changed)


@pytest.mark.parametrize('expected', [[], ['a', 'a']])
def test_invalid_coverage_inventory(expected):
    with pytest.raises(Blocked):
        coverage(expected, {})


@pytest.mark.parametrize('tasks', [{}, {'a': {'status': 'NOT_RUN'}},
                                  {'a': {'status': 'FAILED'}},
                                  {'a': {'status': 'PASS'}, 'extra': {'status': 'PASS'}}])
def test_incomplete_never_passes(tasks):
    report = payload('VALIDATED', ['a'], tasks)
    assert report['status'] != 'VALIDATED'
    with pytest.raises(Blocked):
        require_coverage(report)


def test_coverage_is_recomputed():
    report = payload('PASS', ['a'], {'a': {'status': 'PASS'}})
    report['tasks']['a']['status'] = 'NOT_RUN'
    with pytest.raises(Blocked):
        require_coverage(report)


def test_no_nan_receipts(tmp_path):
    with pytest.raises(ValueError):
        atomic_json(tmp_path / 'bad.json', {'x': float('nan')})
    assert not (tmp_path / 'bad.json').exists()


def test_hash_is_order_independent():
    assert digest({'a': 1, 'b': 2}) == digest({'b': 2, 'a': 1})


def test_task_cache_and_negative_cache(tmp_path):
    root = source(tmp_path)
    count = []
    with Store(root, tmp_path / 'run', {}, plan(), development=True) as store:
        def compute():
            count.append(1)
            return {'status': 'FAILED', 'reason': 'numerical failure'}
        assert store.task({'x': 1}, compute)['status'] == 'FAILED'
        assert store.task({'x': 1}, compute)['status'] == 'FAILED'
        assert len(count) == 1
        assert store.identity['authority'] == 'DEVELOPMENT_ONLY'


def test_retry_failed_is_explicit(tmp_path):
    root, out = source(tmp_path), tmp_path / 'run'
    with Store(root, out, {}, plan(), development=True) as store:
        store.task({'x': 1}, lambda: {'status': 'FAILED'})
    with Store(root, out, {}, plan(), development=True, resume=True, retry_failed=True) as store:
        assert store.task({'x': 1}, lambda: {'status': 'PASS'})['status'] == 'PASS'


def test_cache_corruption_rejected(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {}, plan(), development=True) as store:
        result = store.task({'x': 1}, lambda: {'status': 'PASS', 'value': 1})
        path = store.output / 'tasks' / (result['receipt_key'] + '.json')
        data = read_json(path)
        data['result']['value'] = 999
        atomic_json(path, data)
        with pytest.raises(Blocked):
            store.task({'x': 1}, lambda: {'status': 'PASS'})


def test_output_inside_source_rejected(tmp_path):
    root = source(tmp_path)
    with pytest.raises(Blocked):
        Store(root, root / 'results', {}, plan(), development=True)


def test_strict_rejects_non_git_source(tmp_path):
    with pytest.raises(Blocked):
        Store(source(tmp_path), tmp_path / 'run', {}, plan())


def test_single_writer_lock(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {}, plan(), development=True):
        with pytest.raises(Blocked):
            Store(root, tmp_path / 'run', {}, plan(), development=True, resume=True)
    assert not (tmp_path / 'run/.writer.lock').exists()


def test_source_change_invalidates_run(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {}, plan(), development=True) as store:
        (root / 'hello.py').write_text('x = 2\n')
        with pytest.raises(Blocked):
            store.write('x.json', 'test', {'status': 'PASS'})


def test_changed_configuration_cannot_resume(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {'x': 1}, plan(), development=True):
        pass
    with pytest.raises(Blocked):
        Store(root, tmp_path / 'run', {'x': 2}, plan(), development=True, resume=True)
    assert not (tmp_path / 'run/.writer.lock').exists()


def test_old_artifact_cannot_be_read_into_new_run(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'one', {}, plan(), development=True) as first:
        first.write('x.json', 'test', payload('PASS', ['x'], {'x': {'status': 'PASS'}}))
    with Store(root, tmp_path / 'two', {}, plan(), development=True) as second:
        (second.output / 'x.json').write_bytes((tmp_path / 'one/x.json').read_bytes())
        with pytest.raises(Blocked):
            second.read('x.json', 'PASS')


def test_dependency_mutation_invalidates_child(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {}, plan(), development=True) as store:
        store.write('parent.json', 'test', {'status': 'PASS'})
        store.write('child.json', 'test', {'status': 'PASS'}, ['parent.json'])
        (store.output / 'parent.json').write_text('{}')
        with pytest.raises(Blocked):
            store.read('child.json')


def test_artifact_raw_evidence_mutation_rejected(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {}, plan(), development=True) as store:
        task = store.task({'x': 1}, lambda: {'status': 'PASS'})
        store.write('x.json', 'test', {'status': 'PASS'})
        (store.output / 'tasks' / (task['receipt_key'] + '.json')).write_text('{}')
        with pytest.raises(Blocked):
            store.read('x.json')


def test_schema_one_not_accepted(tmp_path):
    root = source(tmp_path)
    with Store(root, tmp_path / 'run', {}, plan(), development=True) as store:
        record = store.write('x.json', 'test', {'status': 'PASS'})
        record['schema_version'] = 1
        with pytest.raises(Blocked):
            store.verify(record)

import math

import pytest

from mp5d_campaign.core import Blocked, digest
from mp5d_campaign.numerics import approach_routes, diameter, distinct, route


def test_duplicate_root_cannot_be_two_overtones():
    with pytest.raises(Blocked):
        distinct([2.291490104347315-1.8328077227013557j]*2, [1e-10, 1e-10])


def test_error_larger_than_separation_is_unresolved():
    with pytest.raises(Blocked):
        distinct([1-1j, 1.0001-1j], [1e-5, 1e-5])


def test_well_separated_roots():
    distinct([1-1j, 2-1j], [1e-7, 1e-7])


@pytest.mark.parametrize('frequencies, errors', [([complex('nan')], [0]), ([1-1j], [float('nan')]),
                                              ([1-1j], [-1]), ([1-1j, 2-1j], [0])])
def test_nonfinite_and_malformed_identity_vectors(frequencies, errors):
    with pytest.raises(Blocked):
        distinct(frequencies, errors)


def test_tail_uses_full_diameter_not_adjacent_difference():
    assert diameter([[['0', '-1']], [['.75', '-1']], [['1.5', '-1']]]) == [1.5]


def test_decimal_precision_is_preserved():
    rows = [[['1', '-1']], [['1.00000000000000000000000000000000000000000000000001', '-1']],
            [['1.00000000000000000000000000000000000000000000000002', '-1']]]
    assert math.isclose(diameter(rows)[0], 2e-50, rel_tol=1e-12)


def test_nonfinite_decimal_ladder_rejected():
    with pytest.raises(Blocked):
        diameter([[['1', '-1']], [['nan', '-1']], [['1', '-1']]])


def test_two_rungs_are_not_a_convergence_tail():
    with pytest.raises(Blocked):
        diameter([[['1', '-1']], [['1', '-1']]])


def test_routes_are_distinct_even_for_identical_endpoints():
    start = [0, 0, 0]
    paths = approach_routes(start, start, [[0, .42], [-.2, .2], [0, 1.8]], .0025)
    assert len({digest(path) for path in paths}) == 2
    assert all(len(path) >= 3 and path[-1] == start for path in paths)


def test_route_respects_exact_step_and_endpoint():
    path = route([0, 0, 0], [.42, -.2, 1.8], (2, 1, 0), .005)
    assert path[-1] == [.42, -.2, 1.8]
    assert max(max(abs(a-b) for a,b in zip(x,y)) for x,y in zip(path,path[1:])) <= .005000000001


def test_no_fake_independent_directions_on_zero_width_domain():
    with pytest.raises(Blocked):
        approach_routes([0,0,0], [0,0,0], [[0,0],[0,0],[0,0]], .01)


@pytest.mark.parametrize('step', [0, -1, float('nan')])
def test_bad_continuation_step(step):
    with pytest.raises(ValueError):
        route([0,0,0], [1,1,1], (0,1,2), step)


def test_independent_terminal_axes_not_only_early_detours():
    paths=approach_routes([0,0,0],[.42,-.2,1.8],[[0,.42],[-.2,.2],[0,1.8]],.01)
    axes=[{i for i,(a,b) in enumerate(zip(p[-2],p[-1])) if abs(a-b)>1e-12} for p in paths]
    assert axes==[{0},{1}]

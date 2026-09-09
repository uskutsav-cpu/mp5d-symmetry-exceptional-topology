import cmath

import numpy as np
import pytest

from mp5d_science.tracking import (
    Mode,
    TrackingFailure,
    TrackingPolicy,
    assign,
    continue_branches,
    loop_permutation,
    predict,
)


def modes(values):
    return [Mode(complex(x), 1e-14) for x in values]


def test_predictor_rejects_old_off_by_one_coordinates():
    with pytest.raises(ValueError, match="same positive length"):
        predict([[0j], [1 + 0j]], [0.0, 1.0, 2.0], 2.0)


def test_linear_predictor_is_not_shifted():
    assert predict([[2 + 3j], [5 + 3j]], [0.0, 1.0], 2.0)[0] == 8 + 3j


@pytest.mark.parametrize("offset", [0.0, 1e6, 1e10])
def test_centered_quadratic_predictor(offset):
    assert (
        abs(predict([[0], [1], [4]], [offset, offset + 1, offset + 2], offset + 3)[0] - 9) < 1e-12
    )


def test_predictor_duplicate_coordinates():
    with pytest.raises(ValueError):
        predict([[1], [2]], [1, 1], 2)


def test_global_assignment_beats_greedy_with_explicit_nonlocal_use():
    # Greedy maps 0 -> 0.9, then 1 -> -2 (total 3.9).
    # Hungarian maps 0 -> -2, 1 -> 0.9 (total 2.1).
    policy = TrackingPolicy(predictor_rtol=10)
    ordered, report = assign(modes([0, 1]), modes([0.9, -2]), [0, 1], policy, enforce_motion=False)
    assert [m.omega for m in ordered] == [-2, 0.9]
    assert report["permutation"] == [1, 0]


def test_duplicate_root_is_not_an_ep():
    with pytest.raises(TrackingFailure, match="duplicate"):
        assign(modes([1, 2]), modes([1.5, 1.5]), [1.5, 1.5])


def test_ambiguous_assignment_fails():
    with pytest.raises(TrackingFailure, match="ambiguous"):
        assign(modes([-1, 1]), modes([-1j, 1j]), [0, 0])


def test_branch_jump_is_rejected():
    with pytest.raises(TrackingFailure, match="branch jump"):
        assign(modes([1, 3]), modes([2, 4]), [1, 3])


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_nonfinite_residual_never_passes(bad):
    with pytest.raises(TrackingFailure):
        assign(modes([1, 2]), [Mode(1, bad), Mode(2, 1e-15)], [1, 2])


def line_solver(point, seeds):
    t = point[0]
    return modes([2 + 0.2 * t + 1j * t, -2 + 0.1 * t - 0.5j * t])[::-1]


@pytest.mark.parametrize("steps", [4, 8, 16])
def test_forward_reverse_and_step_halving(steps):
    path = np.linspace(0, 1, steps + 1)[:, None]
    start = modes([-2, 2])
    forward = continue_branches(line_solver, path, start)
    assert forward.complete
    assert abs(forward.modes[-1][0].omega - (-1.9 - 0.5j)) < 1e-12
    reverse = continue_branches(line_solver, path[::-1], forward.modes[-1])
    assert reverse.complete
    assert max(abs(a.omega - b.omega) for a, b in zip(reverse.modes[-1], start)) < 1e-12


def test_adaptive_retry_is_real_not_just_documented():
    def solver(point, seeds):
        return modes([-1 + point[0], 1 + point[0]])

    result = continue_branches(solver, [[0], [1]], modes([-1, 1]))
    assert result.complete and result.retries
    assert len(result.parameters) > 2
    assert result.requested_vertices_completed == 2


def test_label_collapse_budget_exhaustion_is_incomplete():
    def collapsed(point, seeds):
        return modes([1, 1])

    result = continue_branches(collapsed, [[0], [1]], modes([0, 2]), TrackingPolicy(max_halvings=3))
    assert not result.complete
    assert len(result.modes) == 1
    assert result.requested_vertices_completed == 1
    assert result.failure


def test_ep_loop_swap_and_two_loop_restoration():
    def solver(point, seeds):
        w = cmath.sqrt(complex(*point))
        return modes([w, -w])

    for turns, expected in [(1, [1, 0]), (2, [0, 1])]:
        theta = np.linspace(0, 2 * np.pi * turns, 48 * turns + 1)
        path = np.column_stack([np.cos(theta), np.sin(theta)])
        result = continue_branches(solver, path, modes([-1, 1]))
        assert result.complete
        assert loop_permutation(result) == expected


def test_eigenvectors_preserve_identity_at_real_part_crossing():
    policy = TrackingPolicy(predictor_rtol=1)
    before = [Mode(-0.01 + 1j, 0, np.array([1, 0])), Mode(0.01 - 1j, 0, np.array([0, 1]))]
    after = [Mode(-0.01 - 1j, 0, np.array([0, 1])), Mode(0.01 + 1j, 0, np.array([1, 0]))]
    ordered, _ = assign(before, after, [0.01 + 1j, -0.01 - 1j], policy)
    assert ordered[0].omega == 0.01 + 1j


def test_nonclosed_loop_rejected():
    r = continue_branches(line_solver, [[0], [0.1]], modes([-2, 2]))
    with pytest.raises(TrackingFailure, match="not closed"):
        loop_permutation(r)

"""Validate the EP machinery on systems whose exceptional structure is exact.

A detector that has never been shown to fire on a real EP cannot support a
negative result.  Every routine in ``mp5d.exceptional.diagnostics`` and
``mp5d.exceptional.verification`` is therefore exercised here against closed-form
problems:

* ``lambda^2 - t``            -- EP2 at ``t = 0`` (square-root branching);
* ``lambda^2 - (t^2 + g^2)``  -- avoided crossing, *no* EP for ``g != 0``;
* ``lambda^2 - t^2``          -- ordinary transversal crossing (diabolic point);
* ``lambda^3 - t``            -- EP3 (cube-root branching, 3-cycle monodromy);
* Jordan blocks               -- defective vs. semisimple degeneracy.

The negative controls matter as much as the positive ones: the tests assert
that the machinery *does not* report an EP for the avoided crossing or for the
semisimple degeneracy.
"""

from __future__ import annotations

import numpy as np
import pytest

from mp5d.exceptional.diagnostics import (
    cauchy_taylor,
    root_separation,
    winding_zero_count,
    zero_count,
)
from mp5d.exceptional.verification import (
    fit_puiseux,
    jordan_chain_matrix,
    monodromy,
    roots_in_disc,
)

# --------------------------------------------------------------------------
# Cauchy machinery
# --------------------------------------------------------------------------


def test_cauchy_taylor_recovers_known_coefficients():
    # f(z) = 3 + 2(z-1) - 5(z-1)^2 + 0.5(z-1)^3
    def f(z):
        w = z - 1.0
        return 3.0 + 2.0 * w - 5.0 * w**2 + 0.5 * w**3

    td = cauchy_taylor(f, 1.0, radius=0.3, n=64, n_coeffs=5)
    assert abs(td.coeffs[0] - 3.0) < 1e-12
    assert abs(td.coeffs[1] - 2.0) < 1e-12
    assert abs(td.coeffs[2] + 5.0) < 1e-12
    assert abs(td.coeffs[3] - 0.5) < 1e-12
    assert abs(td.coeffs[4]) < 1e-12


def test_winding_counts_zeros_and_poles():
    assert abs(winding_zero_count(lambda z: z, 0.0, 1.0) - 1.0) < 1e-9
    assert abs(winding_zero_count(lambda z: z**3, 0.0, 1.0) - 3.0) < 1e-9
    # a pole contributes -1
    assert abs(winding_zero_count(lambda z: 1.0 / z, 0.0, 1.0) + 1.0) < 1e-9
    # zero outside the contour is not counted
    assert abs(winding_zero_count(lambda z: z - 5.0, 0.0, 1.0)) < 1e-9


def test_zero_count_raises_when_unresolved():
    # a zero sitting essentially on the contour is not resolvable
    with pytest.raises((ValueError, FloatingPointError)):
        zero_count(lambda z: z - 1.0, 0.0, 1.0, n=64)


# --------------------------------------------------------------------------
# Root separation: the scale-free EP2 diagnostic
# --------------------------------------------------------------------------


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3, 1e-4])
def test_root_separation_measures_true_distance(d):
    """``|a1/a2|`` must recover the actual distance between two roots."""

    def f(z):
        return (z - 0.0) * (z - d)

    out = root_separation(f, 0.0, radius=min(0.4, 10 * d))
    assert out["separation"] == pytest.approx(d, rel=1e-6)


def test_root_separation_is_invariant_under_rescaling():
    """The whole point: rescaling F must not change the diagnostic.

    ``|dF/domega|`` fails this test, which is why it cannot support a bound.
    """
    d = 1e-3

    def f(z):
        return (z - 0.0) * (z - d)

    def g(z):  # same zeros, wildly different scale and an analytic factor
        return 1e7 * np.exp(3.0 * z) * f(z)

    sep_f = root_separation(f, 0.0, radius=1e-2)["separation"]
    sep_g = root_separation(g, 0.0, radius=1e-2)["separation"]
    assert sep_g == pytest.approx(sep_f, rel=2e-2)

    # ... whereas the naive derivative diagnostic differs by ~1e7
    df = abs(cauchy_taylor(f, 0.0, 1e-2).a1)
    dg = abs(cauchy_taylor(g, 0.0, 1e-2).a1)
    assert dg / df > 1e5


def test_root_separation_large_for_avoided_crossing():
    """Negative control: an avoided crossing has separation bounded below."""
    g = 0.05

    def f(z):  # lambda^2 - g^2 at the closest approach t = 0
        return z**2 - g**2

    out = root_separation(f, -g, radius=0.02)
    assert out["separation"] == pytest.approx(2 * g, rel=1e-6)


# --------------------------------------------------------------------------
# Root finding by contour moments
# --------------------------------------------------------------------------


def test_roots_in_disc_finds_both_roots_of_a_pair():
    d = 1e-3
    rts = roots_in_disc(lambda z: (z - 0.0) * (z - d), 0.0, radius=0.05)
    assert len(rts) == 2
    rts = sorted(rts, key=lambda z: z.real)
    assert abs(rts[0] - 0.0) < 1e-9
    assert abs(rts[1] - d) < 1e-9


def test_roots_in_disc_finds_three_roots():
    t = 1e-3
    rts = roots_in_disc(lambda z: z**3 - t, 0.0, radius=0.5, max_roots=4)
    assert len(rts) == 3
    for z in rts:
        assert abs(z**3 - t) < 1e-9


# --------------------------------------------------------------------------
# Monodromy
# --------------------------------------------------------------------------


def _root_fn_sqrt(radius_disc: float):
    """Roots of lambda^2 - t for t = p[0] + i p[1]."""

    def root_fn(p):
        t = complex(p[0], p[1])
        return roots_in_disc(lambda z: z**2 - t, 0.0, radius=radius_disc)

    return root_fn


def test_monodromy_swaps_roots_around_an_ep2():
    res = monodromy(_root_fn_sqrt(0.5), center_params=[0.0, 0.0], radius=1e-2,
                    n_points=64, n_roots=2)
    assert res.is_swap, f"expected a transposition, got {res.permutation}"
    assert res.closed_error < 1e-6


def test_two_loops_restore_the_roots():
    """One loop swaps, two loops must return the original assignment."""
    root_fn = _root_fn_sqrt(0.5)
    loop = []
    for k in range(129):
        th = 2 * 2.0 * np.pi * k / 128  # two full turns
        loop.append([1e-2 * np.cos(th), 1e-2 * np.sin(th)])
    from mp5d.exceptional.verification import track_roots_on_loop

    hist = track_roots_on_loop(root_fn, loop, n_roots=2)
    assert abs(hist[-1][0] - hist[0][0]) < 1e-6
    assert abs(hist[-1][1] - hist[0][1]) < 1e-6


def test_monodromy_does_not_swap_for_an_avoided_crossing():
    """Negative control: no EP enclosed, so no permutation."""
    g = 0.1

    def root_fn(p):
        t = complex(p[0], p[1])
        return roots_in_disc(lambda z: z**2 - (t**2 + g**2), 0.0, radius=0.6)

    res = monodromy(root_fn, center_params=[0.0, 0.0], radius=1e-2,
                    n_points=64, n_roots=2)
    assert not res.is_swap
    assert res.permutation == [0, 1]


def test_monodromy_three_cycle_for_ep3():
    def root_fn(p):
        t = complex(p[0], p[1])
        return roots_in_disc(lambda z: z**3 - t, 0.0, radius=0.5, max_roots=4)

    res = monodromy(root_fn, center_params=[0.0, 0.0], radius=1e-3,
                    n_points=96, n_roots=3)
    perm = res.permutation
    assert sorted(perm) == [0, 1, 2]
    assert perm != [0, 1, 2], "an EP3 loop must permute the three branches"
    # a 3-cycle has no fixed point
    assert all(perm[k] != k for k in range(3))


# --------------------------------------------------------------------------
# Puiseux
# --------------------------------------------------------------------------


def test_puiseux_detects_square_root_branching():
    ts = np.geomspace(1e-6, 1e-3, 12)
    split = [2.0 * np.sqrt(t) for t in ts]  # |+sqrt(t) - (-sqrt(t))|
    fit = fit_puiseux(ts, split)
    assert fit.exponent == pytest.approx(0.5, abs=1e-6)
    assert fit.prefers_sqrt


def test_puiseux_rejects_square_root_for_linear_splitting():
    ts = np.geomspace(1e-6, 1e-3, 12)
    split = [3.0 * t for t in ts]
    fit = fit_puiseux(ts, split)
    assert fit.exponent == pytest.approx(1.0, abs=1e-6)
    assert not fit.prefers_sqrt


def test_puiseux_on_actual_roots_of_the_ep2_normal_form():
    """End-to-end: extract the splitting from the roots themselves."""
    ts, split = [], []
    for t in np.geomspace(1e-8, 1e-4, 10):
        rts = roots_in_disc(lambda z, tt=t: z**2 - tt, 0.0, radius=0.05)
        ts.append(t)
        split.append(rts[0] - rts[1])
    fit = fit_puiseux(ts, split)
    assert fit.exponent == pytest.approx(0.5, abs=1e-3)
    assert fit.prefers_sqrt


# --------------------------------------------------------------------------
# Jordan chains
# --------------------------------------------------------------------------


def test_jordan_chain_detects_defective_degeneracy():
    """A 2x2 Jordan block is defective: geometric multiplicity 1."""
    A = np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    T = -A                      # T(omega) = omega I - A at omega = 0
    dT = np.eye(2, dtype=complex)
    out = jordan_chain_matrix(T, dT)
    assert out["sigma_min"] < 1e-12
    assert out["defective"], out
    assert out["chain_residual"] < 1e-10


def test_jordan_chain_rejects_semisimple_degeneracy():
    """Negative control: a repeated eigenvalue that is *not* defective.

    This is the case the equal-spin U(2) multiplet falls into, so getting it
    right is directly load-bearing for claim C8.
    """
    A = np.zeros((2, 2), dtype=complex)   # omega = 0 twice, diagonalizable
    T = -A
    dT = np.eye(2, dtype=complex)
    out = jordan_chain_matrix(T, dT)
    assert out["sigma_min"] < 1e-12
    assert not out["defective"], out


def test_jordan_chain_on_a_block_diagonal_pair_is_not_defective():
    """Two *decoupled* blocks sharing an eigenvalue stay semisimple.

    This is exactly the structure of the ``(m1,m2)`` sector decomposition: a
    coincidence of eigenvalues across uncoupled blocks cannot manufacture a
    Jordan chain.  Claim C8 in matrix form.
    """
    A = np.diag([0.0, 0.0]).astype(complex)  # two 1x1 blocks, same eigenvalue
    out = jordan_chain_matrix(-A, np.eye(2, dtype=complex))
    assert not out["defective"]

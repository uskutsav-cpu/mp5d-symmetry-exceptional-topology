"""The full discrete symmetry group of the separated MP5D scalar problem.

Two generators act on ``(parameters) x (sector labels)``:

    E : (a, b, m1, m2) -> (b, a, m2, m1)          exchange
    P : (a, b, m1, m2) -> (-a, -b, -m1, -m2)      simultaneous parity

In adapted coordinates ``s = (a+b)/2``, ``delta = (a-b)/2``:

    E : (s, delta) -> (s, -delta),  (m1, m2) -> (m2, m1)
    P : (s, delta) -> (-s, -delta), (m1, m2) -> (-m1, -m2)

``{1, E, P, EP}`` is a Klein four-group.  Neither generator acts on the
parameters alone -- each also moves the sector -- so a symmetry of the
*parameter space at fixed sector* exists only in sectors the corresponding
group element stabilizes:

    m1 = m2   (diagonal)      : E stabilizes  -> spectrum even in delta
    m1 = -m2  (anti-diagonal)  : EP stabilizes -> spectrum even in s
    (0, 0)                     : all of them  -> even in delta AND in s

This settles claim C13, which conjectured that exceptional sets come in
quadruples under ``s -> -s`` and ``delta -> -delta`` independently.  As stated
that is **false**: it holds exactly in the sector ``(0,0)``, and nowhere else.
"""

from __future__ import annotations

import pytest

from mp5d.radial.qnm import solve_qnm_cf

SOLVE = dict(depth=200, depth_schedule=(100, 200), angular_N=40)
TOL = 1e-12


def w(a, b, mu, m1, m2, ell, seed=1.6 - 0.36j):
    return solve_qnm_cf(a, b, mu, m1, m2, ell, 0, initial_frequency=seed, **SOLVE).omega


# --------------------------------------------------------------------------
# Generators
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m1,m2,ell", [(1, 1, 2), (2, 1, 3), (1, 0, 1), (2, 2, 4)])
def test_parity_generator_P(m1, m2, ell):
    """``(a,b,m1,m2) -> (-a,-b,-m1,-m2)`` leaves the frequency invariant.

    This is the reflection ``phi -> -phi``, ``psi -> -psi`` composed with
    reversing both rotation senses; it is an isometry of the metric that maps
    the sector to its negative.
    """
    lhs = w(0.25, 0.15, 0.4, m1, m2, ell)
    rhs = w(-0.25, -0.15, 0.4, -m1, -m2, ell)
    assert abs(lhs - rhs) < TOL, (lhs, rhs)


@pytest.mark.parametrize("m1,m2,ell", [(2, 1, 3), (1, 0, 1), (2, 0, 2)])
def test_exchange_generator_E(m1, m2, ell):
    lhs = w(0.25, 0.15, 0.4, m1, m2, ell)
    rhs = w(0.15, 0.25, 0.4, m2, m1, ell)
    assert abs(lhs - rhs) < TOL, (lhs, rhs)


# --------------------------------------------------------------------------
# Stabilizers: which sectors inherit a parameter-space symmetry
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m,ell", [(1, 2), (2, 4), (0, 2)])
def test_diagonal_sectors_are_even_in_delta(m, ell):
    """C10: ``m1 = m2`` is stabilized by E, so ``omega`` is even in ``delta``."""
    plus = w(0.25, 0.15, 0.4, m, m, ell)      # s=0.20, delta=+0.05
    minus = w(0.15, 0.25, 0.4, m, m, ell)     # s=0.20, delta=-0.05
    assert abs(plus - minus) < TOL


@pytest.mark.parametrize("m,ell", [(1, 2), (2, 4)])
def test_antidiagonal_sectors_are_even_in_s(m, ell):
    """``m1 = -m2`` is stabilized by ``EP``, so ``omega`` is even in ``s``.

    Previously unrecorded.  ``EP`` maps ``(s, delta) -> (-s, delta)`` and
    ``(m1, m2) -> (-m2, -m1)``, which fixes ``(m, -m)``.
    """
    plus = w(0.25, 0.15, 0.4, m, -m, ell)     # s=+0.20, delta=0.05
    minus = w(-0.15, -0.25, 0.4, m, -m, ell)  # s=-0.20, delta=0.05
    assert abs(plus - minus) < TOL


def test_sector_00_has_the_full_klein_group():
    """C13, corrected: the quadruple structure holds in ``(0,0)`` only."""
    base = w(0.25, 0.15, 0.4, 0, 0, 2)
    assert abs(w(0.15, 0.25, 0.4, 0, 0, 2) - base) < TOL      # delta -> -delta
    assert abs(w(-0.15, -0.25, 0.4, 0, 0, 2) - base) < TOL    # s -> -s
    assert abs(w(-0.25, -0.15, 0.4, 0, 0, 2) - base) < TOL    # both


def test_delta_evenness_fails_in_a_non_diagonal_sector():
    """Negative control: the symmetry is a property of the sector, not of delta.

    Without this, the diagonal-sector result could be an artifact of a solver
    that ignores ``delta`` altogether.
    """
    plus = w(0.25, 0.15, 0.4, 2, 0, 2)
    minus = w(0.15, 0.25, 0.4, 2, 0, 2)
    assert abs(plus - minus) > 1e-3, "off-diagonal sector must NOT be delta-even"


# --------------------------------------------------------------------------
# O7: the equal-spin multiplet count, for all l
# --------------------------------------------------------------------------


def multiplet_at(ell: int, m: int) -> list[tuple[int, int]]:
    """All ``(m1, m2)`` degenerate with each other at ``delta = 0``, fixed ``(l, m)``.

    Degeneracy condition (C5, C6): ``l = 2n + |m1| + |m2|`` with ``n >= 0`` and
    ``m = m1 + m2``.
    """
    out = []
    for m1 in range(-ell - 1, ell + 2):
        m2 = m - m1
        t = abs(m1) + abs(m2)
        if t <= ell and (ell - t) % 2 == 0:
            out.append((m1, m2))
    return out


@pytest.mark.parametrize("ell", list(range(0, 41)))
def test_multiplet_dimension_is_l_plus_1_for_every_allowed_m(ell):
    """O7 discharged for all ``l`` (checked to ``l = 40``; proof in the docs).

    Counting argument: with ``m >= 0``, ``|m1| + |m2|`` takes the value ``m``
    for the ``m+1`` pairs with ``m2 in [0, m]``, and the value ``m + 2j`` for
    exactly two pairs (``m2 = m+j`` and ``m2 = -j``) for each ``j >= 1``.  The
    constraint ``|m1|+|m2| <= l`` with matching parity allows ``j <= (l-m)/2``,
    giving ``(m+1) + 2 (l-m)/2 = l + 1`` -- independent of ``m``.
    """
    ms = [m for m in range(-ell, ell + 1) if (ell - abs(m)) % 2 == 0]
    assert len(ms) == ell + 1
    for m in ms:
        assert len(multiplet_at(ell, m)) == ell + 1, (ell, m)


@pytest.mark.parametrize("ell", list(range(0, 21)))
def test_total_multiplicity_is_l_plus_1_squared(ell):
    """The SO(4) scalar harmonic level has dimension ``(l+1)^2``."""
    ms = [m for m in range(-ell, ell + 1) if (ell - abs(m)) % 2 == 0]
    assert sum(len(multiplet_at(ell, m)) for m in ms) == (ell + 1) ** 2

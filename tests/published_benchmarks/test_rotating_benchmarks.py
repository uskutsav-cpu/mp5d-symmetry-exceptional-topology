"""Regression tests against Huang-Huang, arXiv:2502.11764.

Tolerance policy: the paper prints 5-6 significant digits, and its own two
methods disagree by up to 1.9e-4 at an identical parameter point (Table III
row 1: MM 1.6812-0.347026i vs CFM 1.68112-0.3472i).  Agreement is therefore
judged against its continued-fraction column where available -- the closest
methodological match to Solver A -- at 2e-4, and against its matrix-method
grids at 3e-4, which is the paper's own method spread rather than a relaxed
goalpost.
"""

import pytest

from mp5d.radial.qnm import solve_qnm_cf

SCH = (100, 200, 400)


# (a, b, mu, k, m1, m2, omega_CFM)
TABLE_II = [
    (0.0, 0.0, 0.0, 0, 0, 0, 0.533838 - 0.383387j),
    (0.0, 0.0, 0.0, 0, 1, 0, 1.01602 - 0.362328j),
    (0.0, 0.0, 0.0, 0, 1, 1, 1.51057 - 0.357537j),
    (0.2, 0.0, 0.0, 1, 1, 1, 2.56514 - 0.353072j),
    (0.5, 0.0, 0.0, 1, 1, 0, 2.17828 - 0.340943j),
    # Table II's key is {a, mu, k, m1, m2} with b = 0 fixed by the caption:
    # the printed {0.3, 0.3, 1, 1, 1} is a = 0.3, mu = 0.3, NOT a = b = 0.3.
    (0.3, 0.0, 0.3, 1, 1, 1, 2.60913 - 0.348797j),
]

TABLE_III = [
    (0.2, 0.3, 0.1, 0, 1.68112 - 0.3472j),
    (0.4, 0.2, 0.9, 0, 1.82222 - 0.311152j),
    (0.3, 0.1, 0.1, 1, 2.6344 - 0.349097j),
]


@pytest.mark.parametrize("a,b,mu,k,m1,m2,cfm", TABLE_II)
def test_table_ii(a, b, mu, k, m1, m2, cfm):
    ell = 2 * k + abs(m1) + abs(m2)
    s = solve_qnm_cf(a, b, mu, m1, m2, ell, 0,
                     initial_frequency=cfm * 1.02, depth_schedule=SCH)
    assert s.converged
    assert abs(s.omega - cfm) < 2e-4, f"got {s.omega}, paper CFM {cfm}"
    assert s.omega.imag < 0


@pytest.mark.parametrize("a,b,mu,k,cfm", TABLE_III)
def test_table_iii_two_spin(a, b, mu, k, cfm):
    """The three arbitrary two-spin points -- the core rotating benchmark."""
    ell = 2 * k + 2
    s = solve_qnm_cf(a, b, mu, 1, 1, ell, 0,
                     initial_frequency=cfm * 1.02, depth_schedule=SCH)
    assert s.converged
    assert abs(s.omega - cfm) < 2e-4, f"got {s.omega}, paper CFM {cfm}"
    assert s.ode_residual < 1e-9


@pytest.mark.parametrize("a,b", [(0.2, 0.1), (0.3, 0.1), (0.3, 0.2), (0.4, 0.2)])
def test_table_v_vi_exchange_pairs(a, b):
    """omega(a,b; k,m1,m2) = omega(b,a; k,m2,m1) -- Tables V and VI are transposes."""
    P = solve_qnm_cf(a, b, 0.1, 1, 0, 3, 0,
                     initial_frequency=2.1 - 0.35j, depth_schedule=SCH)
    Q = solve_qnm_cf(b, a, 0.1, 0, 1, 3, 0,
                     initial_frequency=2.1 - 0.35j, depth_schedule=SCH)
    assert abs(P.omega - Q.omega) < 1e-10


@pytest.mark.parametrize("a,b,ell,pub", [
    (0.2, 0.3, 2, 1.6812 - 0.347026j),     # Table IV  (k=0)
    (0.3, 0.4, 2, 1.78795 - 0.329626j),    # Table IV  (k=0)
    (0.2, 0.3, 4, 2.67298 - 0.346169j),    # Table VII (k=1)
    (0.3, 0.4, 4, 2.77349 - 0.332904j),    # Table VII (k=1), r2 = 0.140 > their 0.1 limit
])
def test_equal_spin_family_off_diagonal(a, b, ell, pub):
    s = solve_qnm_cf(a, b, 0.1, 1, 1, ell, 0,
                     initial_frequency=pub * 1.02, depth_schedule=SCH)
    assert abs(s.omega - pub) < 3e-4, f"got {s.omega}, paper MM {pub}"


@pytest.mark.parametrize("x,published_diagonal", [
    (0.2, 1.63886 - 0.351547j),
    (0.3, 1.72958 - 0.340538j),
])
def test_table_vii_diagonal_is_a_source_artifact(x, published_diagonal):
    """The k=1 equal-spin diagonal of Table VII does not match its own trend.

    Our value must (a) disagree with the printed diagonal, and (b) lie between
    the paper's own off-diagonal neighbours -- which is the evidence that the
    printed diagonal carries k=0 values.  Documented in the manifest's
    ``suspicious_entries``.
    """
    neighbours = {
        0.2: (2.59614 - 0.352072j, 2.67298 - 0.346169j),
        0.3: (2.67298 - 0.346169j, 2.77349 - 0.332904j),
    }[x]
    s = solve_qnm_cf(x, x, 0.1, 1, 1, 4, 0,
                     initial_frequency=2.6 - 0.35j, depth_schedule=SCH)
    assert abs(s.omega - published_diagonal) > 0.5, "diagonal unexpectedly matches"
    lo, hi = sorted(c.real for c in neighbours)
    assert lo < s.omega.real < hi, (
        f"our {s.omega.real} does not interpolate the paper's neighbours ({lo}, {hi})"
    )

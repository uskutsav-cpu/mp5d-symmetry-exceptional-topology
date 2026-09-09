import json
from pathlib import Path

import numpy as np
import pytest

from mp5d_science.hyperboloidal import chebyshev_lobatto
from mp5d_science.hyperboloidal import testbed_benchmark as run_testbed_benchmark
from mp5d_science.physical import ContourProblem, solve_a, solve_c
from mp5d_science.radial_polynomial import RadialParameters
from mp5d_science.symbolic_audit import audit_radial_division

ROOT = Path(__file__).resolve().parents[2]
STATIC = RadialParameters("0", "0", "0", 0, 0, 0)


def test_exact_radial_symbolic_remainder():
    result = audit_radial_division()
    assert result["status"] == "PASS"


@pytest.mark.parametrize("degree", range(1, 6))
def test_lobatto_differentiation(degree):
    y, D = chebyshev_lobatto(12)
    assert np.max(abs(D @ (y**degree) - degree * y ** (degree - 1))) < 1e-11


def test_hyperboloidal_testbed_is_not_claimed_as_mp5d():
    result = run_testbed_benchmark()
    assert result["status"] == "TESTBED_PASS"
    assert not result["mp5d_solver_F_validated"]


def test_corrected_c_static_frequency():
    sol = solve_c(STATIC, 0.534 - 0.384j, radial_n=180, angular_n=16, length=80)
    assert sol.converged
    assert abs(sol.omega - (0.533835574268 - 0.383375368512j)) < 1e-9


def test_fixed_bordered_residual_is_locally_holomorphic():
    prob = ContourProblem(STATIC, 0.54 - 0.39j, radial_n=80, angular_n=16, length=80)
    z = 0.541 - 0.391j
    h = 1e-5
    dx = (prob.bordered_residual(z + h) - prob.bordered_residual(z - h)) / (2 * h)
    dy = (prob.bordered_residual(z + 1j * h) - prob.bordered_residual(z - 1j * h)) / (2j * h)
    assert abs(dx - dy) < 1e-5 * max(1, abs(dx))


def test_old_svd_residual_is_a_nonnegative_singular_value():
    prob = ContourProblem(STATIC, 0.54 - 0.39j, radial_n=40, angular_n=16, length=80)
    mat = prob.matrix(0.541 - 0.391j)
    U, S, Vh = np.linalg.svd(mat)
    old = U[:, -1].conj() @ mat @ Vh[-1].conj()
    assert abs(old - S[-1]) < 1e-14
    assert old.real > 0 and abs(old.imag) < 1e-14


def test_quasiresonant_contour_is_rejected():
    params = RadialParameters(".62", ".22", "1.8", 1, 1, 2)
    with pytest.raises(ValueError, match="radiation|suppression"):
        solve_c(params, 1.70277948 - 0.00026860j, angle_deg=80)


def test_ten_historical_candidate_gaps_and_labels_not_fabricated():
    fixture = json.loads((ROOT / "data/regressions/ten_interactions.json").read_text())
    assert len(fixture["records"]) == 10
    for row in fixture["records"]:
        a, b = map(lambda x: complex(*x), row["saved_frequencies"])
        assert abs(abs(a - b) - row["saved_gap"]) < 1e-12
        assert len(set(row["branches"])) == 2
        assert row["expected_classification"] == "REVALIDATE_NOT_ASSUME_PHYSICAL"


@pytest.mark.slow
def test_static_cf_cross_inversions_agree():
    sol = solve_a(STATIC, 0.534 - 0.384j, depth=240, angular_n=16, dps=40)
    assert sol.converged
    assert sol.diagnostics["cross_inversion_consistent"]


@pytest.mark.slow
def test_historical_inversion_specific_root_is_not_accepted():
    p = RadialParameters(".22", ".62", "1.8", 2, 0, 4)
    sol = solve_a(
        p, 2.5094126155849015 - 1.5480379310838372j, overtone=2, depth=200, angular_n=24, dps=45
    )
    assert sol.residual < 1e-25
    assert not sol.converged
    assert not sol.diagnostics["cross_inversion_consistent"]
    assert max(float(v) for v in sol.diagnostics["cross_inversion_residuals"].values()) > 1e-4


def test_legacy_cf_guard_accepts_the_valid_static_root():
    from mp5d_science.physical import check_existing_cf_root

    result = check_existing_cf_root(
        STATIC, 0.533835574268 - 0.383375368512j, depth=240, angular_n=16
    )
    assert result["valid"]


def test_legacy_cf_guard_rejects_the_historical_inversion_specific_root():
    from mp5d_science.physical import check_existing_cf_root

    p = RadialParameters(".22", ".62", "1.8", 2, 0, 4)
    result = check_existing_cf_root(
        p, 2.5094126155849015 - 1.5480379310838372j, depth=200, angular_n=24, overtone=2
    )
    assert not result["valid"]

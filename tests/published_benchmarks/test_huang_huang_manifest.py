"""Second-check the transcribed Huang-Huang manifest against itself.

No solver is involved.  These tests verify that the transcription obeys the
symmetries the *paper's own equations* require, so a typing error in the
manifest is caught before it can be mistaken for a solver failure.
"""

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
MAN = ROOT / "data" / "published_benchmarks" / "huang_huang_2025"


def load(name):
    if not name.endswith(".json"):
        name += ".json"
    return json.loads((MAN / name).read_text())


def index(records):
    return {(r["a"], r["b"]): complex(*r["omega"]) for r in records}


def test_all_manifest_files_exist():
    for n in [
        "table_ii",
        "table_iii",
        "table_iv",
        "table_v",
        "table_vi",
        "table_vii",
        "provenance",
    ]:
        assert (MAN / f"{n}.json").exists(), n


def test_table_v_and_vi_are_exchange_partners():
    """omega(a,b; k,m1,m2) = omega(b,a; k,m2,m1).

    Table V is (k=1, m1=1, m2=0) and Table VI is (k=1, m1=0, m2=1), so
    Table VI must be the transpose of Table V.  This is the sharpest available
    check on the transcription.
    """
    V, W = index(load("table_v")["records"]), index(load("table_vi")["records"])
    checked = 0
    for (a, b), w in V.items():
        if (b, a) in W:
            assert abs(w - W[(b, a)]) < 1e-9, f"V({a},{b}) != VI({b},{a})"
            checked += 1
    assert checked >= 20


@pytest.mark.parametrize("name", ["table_iv", "table_vii"])
def test_m1_equals_m2_tables_are_symmetric_under_a_swap_b(name):
    """With m1 = m2 the paper notes an exact a <-> b symmetry."""
    G = index(load(name)["records"])
    checked = 0
    for (a, b), w in G.items():
        if (b, a) in G:
            assert abs(w - G[(b, a)]) < 1e-9, f"{name}: ({a},{b}) != ({b},{a})"
            checked += 1
    assert checked >= 20


def test_table_iii_agrees_with_table_iv_where_they_overlap():
    """Table III row 1 (a=0.2, b=0.3, mu=0.1, k=0) also appears in Table IV."""
    iii = load("table_iii")["records"][0]
    iv = index(load("table_iv")["records"])
    assert abs(complex(*iii["omega_MM"]) - iv[(0.2, 0.3)]) < 1e-9


def test_inner_horizon_matches_the_papers_own_table_i():
    """Our r2 computation must reproduce the paper's Table I."""
    prov = load("provenance")
    for key, r2 in prov["table_i_inner_horizons"].items():
        a = float(key.split(",")[0].split("=")[1])
        b = float(key.split(",")[1].split("=")[1])
        disc = (1.0 - a * a - b * b) ** 2 - 4 * a * a * b * b
        z2 = 0.5 * ((1.0 - a * a - b * b) - disc**0.5)
        assert abs(z2**0.5 - r2) < 2e-5, f"r2({a},{b}): got {z2**0.5}, paper {r2}"


def test_suspicious_equal_spin_diagonals_are_flagged():
    """The k=1 diagonals are recorded as suspect, not as pass criteria."""
    prov = load("provenance")
    assert prov["suspicious_entries"], "the diagonal anomaly must be documented"
    entry = prov["suspicious_entries"][0]
    assert "Tables V, VI, VII" in entry["where"]
    for r in load("table_vii")["records"]:
        if r["a"] == r["b"]:
            assert "equal_spin_diagonal" in r["flags"]


def test_table_vii_diagonal_duplicates_table_iv_diagonal():
    """Pin the concrete anomaly: k=1 diagonal equals the k=0 diagonal.

    If a future revision of the paper fixes this, the test fails and the
    manifest must be re-derived.
    """
    iv = index(load("table_iv")["records"])
    vii = index(load("table_vii")["records"])
    for x in (0.1, 0.2, 0.3, 0.4):
        assert abs(iv[(x, x)] - vii[(x, x)]) < 1e-9, (
            f"diagonal at a=b={x} no longer duplicates Table IV"
        )

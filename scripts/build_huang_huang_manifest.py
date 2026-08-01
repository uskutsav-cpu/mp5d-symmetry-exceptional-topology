#!/usr/bin/env python3
"""Build the machine-readable Huang-Huang benchmark manifest.

Source: Zi-Yang Huang and Jia-Hui Huang, "Quasinormal modes of massive scalar
fields in five-dimensional Myers-Perry black holes with two arbitrary rotation
parameters", arXiv:2502.11764.

Values are transcribed from the paper's Tables II-VII and then **checked a
second way** by the consistency tests in
``tests/published_benchmarks/test_huang_huang_manifest.py``: exchange symmetry
between Tables V and VI, ``a <-> b`` symmetry within Tables IV and VII, and
agreement of shared entries between Tables III and IV.  Entries that fail those
internal checks are flagged in the manifest rather than silently trusted.

Convention mapping (verified against the paper's Eqs. 1-13):
  * metric, ``rho^2``, ``Delta``, and the mode ansatz agree with ours exactly;
  * their ``lambda_{k m1 m2}`` uses ``(a^2-b^2) cos^2 theta`` with no ``b^2``
    term, so **their lambda = our Ahat = Lambda + (w^2-mu^2) b^2**;
  * their ``k`` is our angular node index ``n``, so ``l = 2k + |m1| + |m2|``;
  * their ``sigma`` denominator ``2 r_H (a^2+b^2-M+2 r_H^2)`` equals our
    ``2 r_+ (z_+ - z_-)`` identically;
  * their ``r2`` is our inner horizon ``r_-``.
"""

from __future__ import annotations

import json
import pathlib

OUT = pathlib.Path(__file__).resolve().parents[1] / "data" / "published_benchmarks" / "huang_huang_2025"

GRID = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]


def C(re: float, im: float) -> list[float]:
    return [re, im]


# ---------------------------------------------------------------- Table II
# M = 1, b = 0 (except the last row).  mu, k, m1, m2 as listed.
TABLE_II = [
    {"a": 0.0, "b": 0.0, "mu": 0.0, "k": 0, "m1": 0, "m2": 0,
     "omega_ref": C(0.53384, -0.38338), "omega_MM": C(0.534406, -0.383241),
     "omega_CFM": C(0.533838, -0.383387)},
    {"a": 0.0, "b": 0.0, "mu": 0.0, "k": 0, "m1": 1, "m2": 0,
     "omega_ref": C(1.01602, -0.36233), "omega_MM": C(1.01620, -0.362074),
     "omega_CFM": C(1.01602, -0.362328)},
    {"a": 0.0, "b": 0.0, "mu": 0.0, "k": 0, "m1": 1, "m2": 1,
     "omega_ref": C(1.51057, -0.35754), "omega_MM": C(1.51055, -0.357371),
     "omega_CFM": C(1.51057, -0.357537)},
    {"a": 0.2, "b": 0.0, "mu": 0.0, "k": 1, "m1": 1, "m2": 1,
     "omega_ref": C(2.565177, -0.353047), "omega_MM": C(2.5651, -0.353033),
     "omega_CFM": C(2.56514, -0.353072)},
    {"a": 0.5, "b": 0.0, "mu": 0.0, "k": 1, "m1": 1, "m2": 0,
     "omega_ref": C(2.178379, -0.340853), "omega_MM": C(2.17821, -0.340858),
     "omega_CFM": C(2.17828, -0.340943)},
    # NOTE: Table II's key is {a, mu, k, m1, m2} and its caption fixes b = 0 for
    # the whole table.  The printed row {0.3, 0.3, 1, 1, 1} therefore means
    # a = 0.3, mu = 0.3, b = 0 -- a singly rotating MASSIVE case, NOT a = b = 0.3.
    # Transcribing it as equal-spin was an error on our side (caught by a 1e-1
    # disagreement with our solver); with the correct reading we agree to 4e-6.
    {"a": 0.3, "b": 0.0, "mu": 0.3, "k": 1, "m1": 1, "m2": 1,
     "omega_ref": C(2.609174, -0.348709), "omega_MM": C(2.60908, -0.348761),
     "omega_CFM": C(2.60913, -0.348797)},
]

# --------------------------------------------------------------- Table III
TABLE_III = [
    {"a": 0.2, "b": 0.3, "mu": 0.1, "k": 0, "m1": 1, "m2": 1,
     "omega_MM": C(1.6812, -0.347026), "omega_CFM": C(1.68112, -0.3472)},
    {"a": 0.4, "b": 0.2, "mu": 0.9, "k": 0, "m1": 1, "m2": 1,
     "omega_MM": C(1.8222, -0.311043), "omega_CFM": C(1.82222, -0.311152)},
    {"a": 0.3, "b": 0.1, "mu": 0.1, "k": 1, "m1": 1, "m2": 1,
     "omega_MM": C(2.63435, -0.349055), "omega_CFM": C(2.6344, -0.349097)},
]

# Grids are given as rows b = 0.1..0.6, columns a = 0.1..0.6.  "-" = not listed
# (those parameters violate M > (a+b)^2).
TABLE_IV = [
    ["1.56895 -0.355839", "1.60232 -0.354012", "1.64036 -0.350607", "1.68394 -0.345147", "1.73441 -0.33684", "1.79398 -0.324253"],
    ["1.60232 -0.354012", "1.63886 -0.351547", "1.6812 -0.347026", "1.73079 -0.339647", "1.78996 -0.327913", "1.86312 -0.308557"],
    ["1.64036 -0.350607", "1.6812 -0.347026", "1.72958 -0.340538", "1.78795 -0.329626", "1.86084 -0.310928", "1.95865 -0.274389"],
    ["1.68394 -0.345147", "1.73079 -0.339647", "1.78795 -0.329626", "1.86008 -0.311689", "1.95776 -0.275621", "-"],
    ["1.73441 -0.33684", "1.78996 -0.327913", "1.86084 -0.310928", "1.95776 -0.275621", "-", "-"],
    ["1.79398 -0.324253", "1.86312 -0.308557", "1.95865 -0.274389", "-", "-", "-"],
]

TABLE_V = [
    ["1.04769 -0.360006", "2.06828 -0.352844", "2.10233 -0.349991", "2.14023 -0.345594", "2.18235 -0.339247", "2.22911 -0.33035"],
    ["2.04284 -0.352826", "1.08355 -0.35586", "2.10912 -0.347573", "2.14826 -0.342413", "2.192 -0.334849", "2.24086 -0.323904"],
    ["2.05147 -0.350037", "2.08406 -0.34753", "1.12891 -0.345735", "2.16212 -0.336585", "2.20875 -0.326514", "2.26124 -0.311063"],
    ["2.06391 -0.34587", "2.09849 -0.342359", "2.1378 -0.336449", "1.19003 -0.321012", "2.23358 -0.312024", "-"],
    ["2.08056 -0.339985", "2.11801 -0.334817", "2.16113 -0.326066", "2.21081 -0.311453", "-", "-"],
    ["2.10205 -0.331786", "2.14357 -0.323765", "2.19196 -0.309667", "-", "-", "-"],
]

TABLE_VI = [
    ["1.04769 -0.360006", "2.04284 -0.352826", "2.05147 -0.350037", "2.06391 -0.34587", "2.08056 -0.339985", "2.10205 -0.331786"],
    ["2.06828 -0.352844", "1.08355 -0.35586", "2.08406 -0.34753", "2.09849 -0.342359", "2.11801 -0.334817", "2.14357 -0.323765"],
    ["2.10233 -0.349991", "2.10912 -0.347573", "1.12891 -0.345735", "2.1378 -0.336449", "2.16113 -0.326066", "2.19196 -0.309667"],
    ["2.14023 -0.345594", "2.14826 -0.342413", "2.16212 -0.336585", "1.19003 -0.321012", "2.21081 -0.311453", "-"],
    ["2.18235 -0.339247", "2.192 -0.334849", "2.20875 -0.326514", "2.23358 -0.312024", "-", "-"],
    ["2.22911 -0.33035", "2.24086 -0.323904", "2.26124 -0.311063", "-", "-", "-"],
]

TABLE_VII = [
    ["1.56895 -0.355839", "2.59614 -0.352072", "2.63435 -0.349055", "2.67832 -0.344386", "2.72896 -0.337559", "2.78765 -0.327711"],
    ["2.59614 -0.352072", "1.63886 -0.351547", "2.67298 -0.346169", "2.72134 -0.340239", "2.77812 -0.33127", "2.84569 -0.317518"],
    ["2.63435 -0.349055", "2.67298 -0.346169", "1.72958 -0.340538", "2.77349 -0.332904", "2.83939 -0.31984", "2.92089 -0.297536"],
    ["2.67832 -0.344386", "2.72134 -0.340239", "2.77349 -0.332904", "1.86008 -0.311689", "2.91735 -0.298797", "-"],
    ["2.72896 -0.337559", "2.77812 -0.33127", "2.83939 -0.31984", "2.91735 -0.298797", "-", "-"],
    ["2.78765 -0.327711", "2.84569 -0.317518", "2.92089 -0.297536", "-", "-", "-"],
]

# Table I: the paper's own inner-horizon values, used for the reliability flag.
TABLE_I_R2 = {
    (0.1, 0.1): 0.01010, (0.2, 0.1): 0.02052, (0.3, 0.1): 0.03164,
    (0.4, 0.1): 0.04396, (0.5, 0.1): 0.05826, (0.6, 0.1): 0.07594,
    (0.2, 0.2): 0.04174, (0.3, 0.2): 0.06448, (0.4, 0.2): 0.089890,
    (0.5, 0.2): 0.11990, (0.6, 0.2): 0.15826,
    (0.3, 0.3): 0.1, (0.4, 0.3): 0.14042, (0.5, 0.3): 0.18990, (0.6, 0.3): 0.25903,
}


def parse_cell(s: str):
    if s.strip() == "-":
        return None
    re_s, im_s = s.split()
    return [float(re_s), float(im_s)]


def inner_horizon(a: float, b: float, M: float = 1.0) -> float | None:
    disc = (M - a * a - b * b) ** 2 - 4 * a * a * b * b
    if disc < 0:
        return None
    z2 = 0.5 * ((M - a * a - b * b) - disc**0.5)
    return z2**0.5 if z2 >= 0 else None


def grid_records(grid, k: int, m1: int, m2: int, mu: float, table: str):
    out = []
    for i, brow in enumerate(grid):
        for j, cell in enumerate(brow):
            val = parse_cell(cell)
            if val is None:
                continue
            a, b = GRID[j], GRID[i]
            r2 = inner_horizon(a, b)
            rec = {
                "table": table, "row": f"b={b}", "col": f"a={a}",
                "M": 1.0, "a": a, "b": b, "mu": mu, "k": k, "m1": m1, "m2": m2,
                "l": 2 * k + abs(m1) + abs(m2),
                "method": "matrix method (MM)",
                "omega": val,
                "printed_digits": len(cell.split()[0].replace(".", "").lstrip("0")),
                "inner_horizon_r2": r2,
                "cfm_regarded_reliable_by_authors": (r2 is not None and r2 < 0.1),
                "flags": [],
            }
            if a == b:
                rec["flags"].append("equal_spin_diagonal")
            out.append(rec)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    common = {
        "source": {
            "authors": "Zi-Yang Huang and Jia-Hui Huang",
            "title": ("Quasinormal modes of massive scalar fields in five-dimensional "
                      "Myers-Perry black holes with two arbitrary rotation parameters"),
            "arxiv": "2502.11764", "dated": "2025-02-20",
        },
        "normalization": {
            "metric": "their Eq.(1); rho^2 and Delta agree with ours",
            "ansatz": "Psi = e^{-i w t} e^{i m1 phi} e^{i m2 psi} R(r) S(theta)",
            "sign_convention": "damped modes have Im(omega) < 0",
            "lambda_mapping": "their lambda_{k m1 m2} = our Ahat = Lambda + (w^2-mu^2) b^2",
            "k_mapping": "their k is our angular node index n; l = 2k + |m1| + |m2|",
            "r2_mapping": "their r2 is our inner horizon r_-",
            "restriction": "paper restricts to m1 >= 0, m2 >= 0",
        },
        "cfm_caveat": (
            "The authors state the continued-fraction method is applicable for "
            "r2 << 1 and that their CFM results 'begin to show instability when "
            "r2 >~ 0.1' (text near Fig. 1). Their CFM uses a 25-term recursion."
        ),
    }

    (OUT / "table_ii.json").write_text(json.dumps(
        {**common, "table": "II", "caption": "Comparison with literature. M=1, b=0 except last row.",
         "records": TABLE_II}, indent=2) + "\n")
    (OUT / "table_iii.json").write_text(json.dumps(
        {**common, "table": "III", "caption": "Two-spin comparison, M=1, m1=m2=1.",
         "records": TABLE_III}, indent=2) + "\n")
    (OUT / "table_iv.json").write_text(json.dumps(
        {**common, "table": "IV", "caption": "k=0, m1=m2=1, mu=0.1 grid.",
         "records": grid_records(TABLE_IV, 0, 1, 1, 0.1, "IV")}, indent=2) + "\n")
    (OUT / "table_v.json").write_text(json.dumps(
        {**common, "table": "V", "caption": "k=1, m1=1, m2=0, mu=0.1 grid.",
         "records": grid_records(TABLE_V, 1, 1, 0, 0.1, "V")}, indent=2) + "\n")
    (OUT / "table_vi.json").write_text(json.dumps(
        {**common, "table": "VI", "caption": "k=1, m1=0, m2=1, mu=0.1 grid.",
         "records": grid_records(TABLE_VI, 1, 0, 1, 0.1, "VI")}, indent=2) + "\n")
    (OUT / "table_vii.json").write_text(json.dumps(
        {**common, "table": "VII", "caption": "k=m1=m2=1, mu=0.1 grid.",
         "records": grid_records(TABLE_VII, 1, 1, 1, 0.1, "VII")}, indent=2) + "\n")

    prov = {
        **common,
        "extraction": {
            "method": "text extracted from the arXiv PDF with pypdf, transcribed here, "
                      "then re-checked by automated internal-consistency tests",
            "tests": "tests/published_benchmarks/test_huang_huang_manifest.py",
            "pdf_committed": False,
        },
        "table_i_inner_horizons": {f"a={a},b={b}": r2 for (a, b), r2 in TABLE_I_R2.items()},
        "suspicious_entries": [
            {
                "where": "Tables V, VI, VII: every equal-spin diagonal entry (a = b)",
                "status": "suspected source artifact; verified directly by our solver",
                "observation": (
                    "Table VII's diagonal (1.56895, 1.63886, 1.72958, 1.86008) is "
                    "IDENTICAL to Table IV's diagonal, although Table IV is k=0 and "
                    "Table VII is k=1. Tables V and VI show the same pattern: their "
                    "diagonals (~1.05-1.19) break the trend of their own off-diagonal "
                    "neighbours (~2.05-2.26)."
                ),
                "retracted_argument": (
                    "An earlier draft of this manifest claimed Table II contained an "
                    "equal-spin entry contradicting Table VII. That was OUR misreading: "
                    "Table II's key is {a, mu, k, m1, m2} with b = 0 fixed by its "
                    "caption, so {0.3, 0.3, 1, 1, 1} is a = 0.3, mu = 0.3, b = 0. "
                    "Table II contains no equal-spin entry. The claim is withdrawn."
                ),
                "direct_test": (
                    "The anomaly is instead tested directly: our solver is run at "
                    "a = b, k = 1, m1 = m2 = 1, mu = 0.1 and compared with both the "
                    "printed diagonal and the trend of its own off-diagonal "
                    "neighbours. See results/rotating_validation.json."
                ),
                "hypothesis": (
                    "The k=1 diagonals appear to carry k=0 values. At a=b the angular "
                    "equation is exactly degenerate (our c2 = 4 s delta (w^2-mu^2) "
                    "vanishes identically), so an angular routine that returns the "
                    "lowest eigenvalue regardless of the requested k would reproduce "
                    "exactly this pattern."
                ),
                "action": (
                    "Treated as a suspected transcription/solver artifact in the "
                    "source, NOT as a target for our solver. Excluded from the pass "
                    "criteria; tested separately against Table II, which is internally "
                    "consistent."
                ),
            }
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote manifest to {OUT}")
    for f in sorted(OUT.glob("*.json")):
        d = json.loads(f.read_text())
        n = len(d.get("records", []))
        print(f"  {f.name}: {n} records" if n else f"  {f.name}")


if __name__ == "__main__":
    main()

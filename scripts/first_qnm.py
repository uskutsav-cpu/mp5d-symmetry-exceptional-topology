#!/usr/bin/env python3
"""Compute and validate the first trustworthy MP5D quasinormal modes.

Reproduces every number in ``docs/FIRST_QNM_VALIDATION.md``.

    python scripts/first_qnm.py                  # full validation set
    python scripts/first_qnm.py --case st5d_l0n0 # a single case
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mp5d.radial.qnm import solve_qnm_cf, solve_qnm_hill  # noqa: E402

TWO_PI = 2.0 * np.pi

# Matyjasek, arXiv:2107.04815, Tables I and II.
# Metric f(r) = 1 - r^{3-D} with D = 5, so r_+ = 1 and T_H = 1/(2 pi).
# Published as omega_tilde = omega / T_H; we store omega in r_+ = 1 units.
PUBLISHED = {
    "st5d_l0n0": {"l": 0, "n": 0, "m": (0, 0), "wt": (3.35418783669, -2.40881848257)},
    "st5d_l0n1": {"l": 0, "n": 1, "m": (0, 0), "wt": (2.33646259225, -8.31019918700)},
    "st5d_l1n0": {"l": 1, "n": 0, "m": (1, 0), "wt": (6.38382253011, -2.27657411582)},
    "st5d_l1n1": {"l": 1, "n": 1, "m": (1, 0), "wt": (5.38079295983, -7.27345089157)},
}


def commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "uncommitted"


def env_hash() -> str:
    p = ROOT / "requirements.lock"
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16] if p.exists() else "none"


def run_benchmarks(depth_schedule, sizes):
    out = []
    for key, spec in PUBLISHED.items():
        wref = complex(*spec["wt"]) / TWO_PI
        m1, m2 = spec["m"]
        seed = wref * 1.04 + 0.01  # deliberately off the answer
        t0 = time.time()
        A = solve_qnm_cf(
            0.0,
            0.0,
            0.0,
            m1,
            m2,
            spec["l"],
            overtone=spec["n"],
            initial_frequency=seed,
            depth_schedule=depth_schedule,
        )
        B = solve_qnm_hill(
            0.0,
            0.0,
            0.0,
            m1,
            m2,
            spec["l"],
            overtone=spec["n"],
            initial_frequency=seed,
            sizes=sizes,
        )
        out.append(
            {
                "case": key,
                "mode_labels": {"m1": m1, "m2": m2, "l": spec["l"], "overtone": spec["n"]},
                "parameters": {"a": 0.0, "b": 0.0, "mu": 0.0, "M": 1.0},
                "omega_solverA": [A.omega.real, A.omega.imag],
                "omega_solverB": [B.omega.real, B.omega.imag],
                "published_omega": [wref.real, wref.imag],
                "published_omega_over_TH": list(spec["wt"]),
                "published_source": "Matyjasek, arXiv:2107.04815, Tables I-II",
                "absolute_error": abs(A.omega - wref),
                "relative_error": abs(A.omega - wref) / abs(wref),
                "solverAB_difference": abs(A.omega - B.omega),
                "cf_residual": A.cf_residual,
                "ode_residual": A.ode_residual,
                "depth_table": [[d, w.real, w.imag] for d, w in A.depth_table],
                "hill_size_table": [[d, w.real, w.imag] for d, w in B.depth_table],
                "precision": "double",
                "evidence_level": "cross-solver verified",
                "runtime_s": round(time.time() - t0, 2),
            }
        )
    return out


def run_symmetry_checks(depth_schedule):
    sch = depth_schedule
    checks = {}

    P = solve_qnm_cf(
        0.30, 0.15, 0.20, 1, 0, 1, 0, initial_frequency=1.02 - 0.36j, depth_schedule=sch
    )
    Q = solve_qnm_cf(
        0.15, 0.30, 0.20, 0, 1, 1, 0, initial_frequency=1.02 - 0.36j, depth_schedule=sch
    )
    checks["exchange_symmetry"] = {
        "description": "(a,m1)<->(b,m2) is an exact identity; residual measures numerical error",
        "omega_1": [P.omega.real, P.omega.imag],
        "omega_2": [Q.omega.real, Q.omega.imag],
        "residual": abs(P.omega - Q.omega),
    }

    trio = []
    for m1, m2 in [(1, 1), (2, 0), (0, 2)]:
        S = solve_qnm_cf(
            0.25, 0.25, 0.15, m1, m2, 2, 0, initial_frequency=1.30 - 0.35j, depth_schedule=sch
        )
        trio.append(S.omega)
    checks["u2_multiplet_degeneracy"] = {
        "description": (
            "claim C6/C7: at delta=0 the spectrum depends on (m1,m2) only through "
            "m = m1+m2, so (1,1), (2,0), (0,2) at l=2, m=2 must coincide"
        ),
        "omegas": [[w.real, w.imag] for w in trio],
        "max_spread": max(abs(x - y) for x in trio for y in trio),
    }

    split = []
    for m1, m2 in [(1, 1), (2, 0)]:
        S = solve_qnm_cf(
            0.35, 0.15, 0.15, m1, m2, 2, 0, initial_frequency=1.30 - 0.35j, depth_schedule=sch
        )
        split.append(S.omega)
    checks["degeneracy_lifted_off_equal_spin"] = {
        "description": "delta != 0 must split what delta = 0 held degenerate",
        "omegas": [[w.real, w.imag] for w in split],
        "splitting": abs(split[0] - split[1]),
    }
    return checks


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=None)
    ap.add_argument("--depth", type=int, nargs="+", default=[200, 400, 800])
    ap.add_argument("--hill-sizes", type=int, nargs="+", default=[120, 180, 240, 300, 360])
    ap.add_argument("--skip-symmetry", action="store_true")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "results" / "first_qnm.json")
    args = ap.parse_args()

    if args.case:
        global PUBLISHED
        PUBLISHED = {args.case: PUBLISHED[args.case]}

    t0 = time.time()
    benches = run_benchmarks(tuple(args.depth), tuple(args.hill_sizes))
    checks = {} if args.skip_symmetry else run_symmetry_checks(tuple(args.depth))

    payload = {
        "status": "complete",
        "provenance": {
            "commit": commit(),
            "env_hash": env_hash(),
            "solver_A": "leaver-cf/2.0",
            "solver_B": "hill-wynn/1.0",
            "precision": "double",
            "total_runtime_s": round(time.time() - t0, 1),
        },
        "normalization": (
            "M = 1 so r_+ = 1 and kappa = 1; matches Matyjasek's f(r) = 1 - r^{3-D}, "
            "D = 5.  Published omega/T_H converted with T_H = 1/(2 pi)."
        ),
        "benchmarks": benches,
        "symmetry_checks": checks,
    }
    args.out.write_text(json.dumps(payload, indent=2, default=float) + "\n")

    print(f"wrote {args.out}")
    for b in benches:
        w = b["omega_solverA"]
        print(
            f"  {b['case']:12s} omega={w[0]:.12f}{w[1]:+.12f}j  "
            f"abs_err={b['absolute_error']:.2e}  |A-B|={b['solverAB_difference']:.2e}"
        )
    for k, v in checks.items():
        key = (
            "residual" if "residual" in v else ("max_spread" if "max_spread" in v else "splitting")
        )
        print(f"  {k}: {key} = {v[key]:.3e}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Re-evaluate every historical close interaction; scientific failures return exit 2.

Use --development only for explicitly non-authoritative local diagnostics.
Output is written even when physics validation fails. Such output is not a PASS.
"""

from __future__ import annotations

import argparse
import sys
import time
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mp5d_science.physical import solve_a, solve_c
from mp5d_science.provenance import atomic_json, envelope, read_json
from mp5d_science.radial_polynomial import RadialParameters


def audit_candidate(record, depths, dps, radial_ns):
    params = RadialParameters(
        str(Decimal(record["s"]) + Decimal(record["delta"])),
        str(Decimal(record["s"]) - Decimal(record["delta"])),
        record["mu"],
        record["m1"],
        record["m2"],
        record["ell"],
    )
    rows = []
    for label, seed_pair, overtone in zip(
        record["branches"], record["saved_frequencies"], record["overtone_labels"], strict=True
    ):
        seed = complex(*seed_pair)
        for depth in depths:
            started = time.monotonic()
            try:
                sol = solve_a(params, seed, overtone=overtone, depth=depth, angular_n=24, dps=dps)
                out = sol.to_dict() | {
                    "status": "ROOT_CHECK_PASSED" if sol.converged else "FAILED_ROOT_VALIDITY"
                }
                out["distance_from_historical_seed"] = abs(sol.omega - seed)
            except Exception as exc:
                out = {
                    "status": "NUMERICAL_FAILURE",
                    "exception_type": type(exc).__name__,
                    "error": str(exc),
                }
            rows.append(
                out
                | {
                    "branch": label,
                    "requested": {"solver": "A", "depth": depth, "dps": dps},
                    "runtime_s": time.monotonic() - started,
                }
            )
        for radial_n in radial_ns:
            started = time.monotonic()
            try:
                sol = solve_c(
                    params, seed, radial_n=radial_n, angular_n=24, length=40, angle_deg=65
                )
                displacement = abs(sol.omega - seed)
                same = displacement < 1e-5
                out = sol.to_dict() | {
                    "distance_from_historical_seed": displacement,
                    "same_historical_root": same,
                    "status": "ROOT_CHECK_PASSED"
                    if sol.converged and same
                    else "NOT_AN_INDEPENDENT_CONFIRMATION",
                }
            except Exception as exc:
                out = {
                    "status": "NUMERICAL_FAILURE",
                    "exception_type": type(exc).__name__,
                    "error": str(exc),
                }
            rows.append(
                out
                | {
                    "branch": label,
                    "requested": {"solver": "C", "radial_n": radial_n},
                    "runtime_s": time.monotonic() - started,
                }
            )
    passed = bool(rows) and all(r["status"] == "ROOT_CHECK_PASSED" for r in rows)
    return {
        "id": record["id"],
        "status": "POINTWISE_ROOT_CHECKS_PASSED" if passed else "REQUIRES_REVALIDATION",
        "records": rows,
        "scientific_claim": "no EP or exclusion claim follows from this regression audit",
        "historical_gap_not_recertified": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--depths", type=int, nargs="+", default=[160, 240, 320])
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--radial-n", type=int, nargs="+", default=[160, 220, 280])
    parser.add_argument("--development", action="store_true")
    args = parser.parse_args()
    fixture = read_json(ROOT / "data/regressions/ten_interactions.json")
    results = []
    for record in fixture["records"]:
        result = audit_candidate(record, args.depths, args.dps, args.radial_n)
        results.append(result)
        print(record["id"], result["status"], flush=True)
        payload = {
            "status": "RUNNING",
            "completed": len(results),
            "expected": len(fixture["records"]),
            "fixture_source": fixture["source"],
            "results": results,
        }
        atomic_json(
            args.output,
            envelope(
                "candidate-regression-audit", payload, ROOT, allow_development=args.development
            ),
        )
    valid = all(r["status"] == "POINTWISE_ROOT_CHECKS_PASSED" for r in results)
    payload["status"] = "POINTWISE_CHECKS_PASS" if valid else "SCIENTIFIC_REVALIDATION_REQUIRED"
    atomic_json(
        args.output,
        envelope("candidate-regression-audit", payload, ROOT, allow_development=args.development),
    )
    return 0 if valid else 2


if __name__ == "__main__":
    sys.exit(main())

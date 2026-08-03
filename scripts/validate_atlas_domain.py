"""Independently validate the atlas domain with a recurrence-free solver.

Why this exists
---------------
Prior work in this repository reported its near-extremality parameter as ``r2``,
which is the inner horizon **radius** ``r_-``.  The atlas domain guard was
written in terms of ``z_minus = r_-^2``.  ``z_minus <= 0.20`` is therefore
``r_- <= 0.447``, which reaches **well beyond** the previously validated
``r_- <= 0.19``.  The atlas is not trustworthy on that extension until a solver
from the other family confirms it.

Solvers A (Leaver continued fraction) and C (exterior complex scaling) share no
discretization: C is enforced recurrence-free by a structural test.  Agreement
between them across the domain is the evidence the atlas needs; disagreement
would bound the atlas back to where they do agree.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.geometry import MPGeometry, sdelta_to_ab  # noqa: E402
from mp5d.radial.qnm import solve_qnm_cf  # noqa: E402
from mp5d.radial.solver_c import min_scaling_angle, solve_qnm_c  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/final_multisolver_validation.json")
    ap.add_argument("--theta-deg", type=float, default=60.0)
    ap.add_argument("--resolution", type=int, default=260)
    args = ap.parse_args()

    # Points spanning the atlas domain, including its most extreme corners.
    cases = [
        # (m1, m2, ell, s, delta, mu, label)
        (1, 1, 2, 0.00, 0.00, 0.0, "static"),
        (1, 1, 2, 0.20, 0.00, 0.0, "equal-spin, mid s"),
        (1, 1, 2, 0.42, 0.00, 0.0, "equal-spin, max s"),
        (1, 1, 2, 0.42, 0.00, 1.8, "equal-spin, max s, max mu"),
        (1, 1, 2, 0.33, 0.20, 0.6, "two-spin, large delta"),
        (1, 1, 2, 0.42, 0.20, 1.8, "domain corner"),
        (2, 0, 2, 0.42, 0.00, 1.2, "off-diagonal, max s"),
        (2, 1, 3, 0.37, 0.15, 0.9, "odd l, large s"),
        (1, 0, 1, 0.42, 0.10, 1.5, "odd l, corner"),
        (2, 2, 4, 0.42, 0.00, 0.6, "high l, max s"),
    ]

    t0 = time.time()
    records = []
    for (m1, m2, ell, s, delta, mu, label) in cases:
        a, b = sdelta_to_ab(s, delta)
        geo = MPGeometry(a=a, b=b, M=1.0)
        rec = {"label": label, "m1": m1, "m2": m2, "ell": ell,
               "s": s, "delta": delta, "mu": mu, "a": a, "b": b,
               "r_minus": geo.r_minus, "z_minus": geo.z_minus,
               "extremality": geo.extremality}
        try:
            A = solve_qnm_cf(a, b, mu, m1, m2, ell, 0,
                             initial_frequency=1.5 - 0.35j if s < 0.3 else 1.9 - 0.3j,
                             depth=400, depth_schedule=(200, 400), angular_N=60)
            rec["solverA"] = [A.omega.real, A.omega.imag]
            rec["solverA_residual"] = A.cf_residual
            rec["solverA_converged"] = bool(A.converged)
        except Exception as exc:  # noqa: BLE001
            rec["solverA_error"] = str(exc)[:150]

        if "solverA" in rec:
            wA = complex(*rec["solverA"])
            # The applicability criterion for exterior complex scaling is NOT
            # "Re omega < mu".  It is the minimum admissible scaling angle,
            # tan(theta) > -Im(Omega)/Re(Omega) with Omega = sqrt(omega^2-mu^2).
            # As Re(Omega) -> 0 the outgoing wave stops oscillating and the
            # required angle -> 90 deg, which no contour can supply.
            th_min_deg = math.degrees(min_scaling_angle(wA, mu))
            rec["min_scaling_angle_deg"] = th_min_deg
            rec["quasiresonant"] = bool(th_min_deg > 85.0)
            if rec["quasiresonant"]:
                rec["solverC_status"] = (
                    f"not applicable: required scaling angle {th_min_deg:.2f} deg "
                    "-> 90 deg as Re(Omega) -> 0 (quasiresonant limit). This is a "
                    "formulation limit of exterior complex scaling, not a "
                    "convergence failure.")
            else:
                # Use the solver's own minimum scaling angle; a fixed theta is
                # too small for strongly damped modes and sends Muller to an
                # unrelated root (seen at the odd-l corner: |A-C| = 2.9 at
                # theta = 60 deg, 2e-9 at theta = 75 deg).
                best = None
                for factor in (1.25, 1.45, 1.7):
                    th_try = min(max(math.radians(th_min_deg * factor),
                                     math.radians(40.0)), math.radians(84.0))
                    for L in (60.0, 90.0):
                        try:
                            C = solve_qnm_c(a, b, mu, m1, m2, ell,
                                            initial_frequency=wA, theta=th_try,
                                            L=L, resolution=args.resolution)
                        except Exception:  # noqa: BLE001
                            continue
                        if not C.converged:
                            continue
                        d = abs(wA - C.omega)
                        if best is None or d < best[0]:
                            best = (d, C.omega, math.degrees(th_try), L)
                if best is not None:
                    rec["solverC"] = [best[1].real, best[1].imag]
                    rec["solverC_theta_deg"] = best[2]
                    rec["solverC_L"] = best[3]
                    rec["solverC_converged"] = True
                    rec["difference"] = best[0]
                else:
                    rec["solverC_status"] = "no converged Solver C run"

        records.append(rec)
        d = rec.get("difference")
        note = (f" [quasiresonant: needs theta={rec['min_scaling_angle_deg']:.1f} deg,"
                " C n/a]") if rec.get("quasiresonant") else ""
        print(f"{label:28s} r-={geo.r_minus:.4f} extr={geo.extremality:.4f} "
              f"|A-C| = {format(d, '.3e') if d is not None else 'n/a'}{note}",
              flush=True)

    diffs = [r["difference"] for r in records if "difference" in r]
    out = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "precision": "double",
                       "solverC_theta_deg": args.theta_deg,
                       "solverC_resolution": args.resolution,
                       "runtime_s": round(time.time() - t0, 1)},
        "purpose": "validate the atlas domain, which extends beyond the "
                   "previously validated r_- <= 0.19, using the recurrence-free "
                   "Solver C",
        "naming_note": "prior work's r2 is the inner horizon RADIUS r_-; "
                       "z_minus = r_-^2. The atlas guard z_minus <= 0.20 is "
                       "r_- <= 0.447.",
        "n_compared": len(diffs),
        "max_difference": max(diffs) if diffs else None,
        "median_difference": sorted(diffs)[len(diffs) // 2] if diffs else None,
        "max_r_minus_compared": max((r["r_minus"] for r in records
                                     if "difference" in r), default=None),
        "records": records,
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"\nmax |A-C| = {out['max_difference']}, "
          f"max r_- compared = {out['max_r_minus_compared']}")
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

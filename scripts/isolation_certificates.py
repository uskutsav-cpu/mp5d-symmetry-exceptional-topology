"""Exactly invariant isolation certificates from the argument principle.

Why this exists
---------------
The coefficient ratio ``|a_1/a_2|`` is only *asymptotically* insensitive to the
normalisation of the spectral condition.  Writing ``F -> gF`` with ``g``
holomorphic and non-vanishing, and expanding about a simple root,

    a~_1 = g_0 a_1 ,      a~_2 = g_0 a_2 + g_1 a_1 ,

so that

    a~_1/a~_2 = (a_1/a_2) / (1 + (g_1/g_0)(a_1/a_2)) ,

which returns ``a_1/a_2`` only in the limit ``a_1/a_2 -> 0``.  The relative
distortion is ``O(|g'/g| d)`` with ``d`` the root separation, and it is
measurable: for ``g = A e^{k w}`` the ratio is exactly
``|1 + k d|^{-1}``, which we verify numerically.

By contrast the **number and location of the zeros** of ``F`` inside a contour
are unchanged by ``F -> gF``, because a non-vanishing ``g`` contributes neither
zeros nor poles.  The argument principle therefore yields a diagnostic with no
normalisation dependence whatever, and it is the one used for the reported
bound.

The certificate
---------------
For a tracked root ``w`` and radius ``R``, evaluate

    n(R) = (1/2 pi i) oint_{|z-w|=R} F'/F dz .

If ``n(R) = 1`` then ``F`` has exactly one zero within distance ``R`` of ``w``,
so no partner root lies closer than ``R`` and no coalescence can be in progress
there.  The *isolation radius* is the largest tested ``R`` for which
``n(R) = 1``; the minimum isolation radius over the search domain is a lower
bound on the spectral separation that does not depend on any root finder, on
any local expansion, or on the normalisation of ``F``.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.exceptional.diagnostics import winding_zero_count  # noqa: E402
from mp5d.exceptional.ep_solver import spectral_condition  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def isolation_radius(rec: dict, radii: list[float], depth: int,
                     n_contour: int) -> dict:
    """Largest tested radius whose contour encloses exactly one zero."""
    w = complex(*rec["omega_a"])
    ell = rec["ell_a"]
    inv = rec["N_a"]          # the branch's own inversion; see collision_search

    def F(z: complex) -> complex:
        return spectral_condition(z, rec["s"], rec["delta"], rec["mu"],
                                  rec["m1"], rec["m2"], ell, inv, depth=depth)

    best, trace = 0.0, []
    for R in sorted(radii):
        try:
            wind = winding_zero_count(F, w, R, n=n_contour)
        except Exception as exc:  # noqa: BLE001
            trace.append({"radius": R, "error": str(exc)[:90]})
            break
        k = round(wind)
        resolved = abs(wind - k) < 0.08
        trace.append({"radius": R, "winding": wind, "count": k if resolved else None})
        if not resolved:
            break
        if k == 1:
            best = R
        else:
            break          # a second zero has entered; stop growing the contour
    return {"isolation_radius": best, "trace": trace, "inversion": inv,
            "omega": [w.real, w.imag]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="results/all_collision_candidates.json")
    ap.add_argument("--out", default="results/isolation_certificates.json")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--subsample", type=int, default=25)
    ap.add_argument("--seed", type=int, default=20260803)
    ap.add_argument("--depth", type=int, default=200)
    ap.add_argument("--n-contour", type=int, default=192)
    ap.add_argument("--radii", type=float, nargs="+",
                    default=[0.05, 0.10, 0.15, 0.20, 0.25])
    args = ap.parse_args()

    t0 = time.time()
    data = json.loads(pathlib.Path(args.candidates).read_text())
    cands = data["candidates"]

    idx = list(range(min(args.top, len(cands))))
    rng = random.Random(args.seed)
    rest = list(range(len(idx), len(cands)))
    rng.shuffle(rest)
    idx += rest[:args.subsample]

    records = []
    for i in idx:
        c = cands[i]
        out = isolation_radius(c, args.radii, args.depth, args.n_contour)
        out["parameters"] = {k: c[k] for k in ("m1", "m2", "s", "delta", "mu")}
        out["branches"] = [c["branch_a"], c["branch_b"]]
        out["tracked_gap"] = c["gap"]
        records.append(out)
        print(f"  [{i}] ({c['m1']},{c['m2']}) gap={c['gap']:.4f} "
              f"isolation R>={out['isolation_radius']:.3f}", flush=True)

    radii = [r["isolation_radius"] for r in records]
    ok = [r for r in radii if r > 0]
    result = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "solver": "leaver-cf",
                       "precision": "double", "depth": args.depth,
                       "n_contour": args.n_contour,
                       "runtime_s": round(time.time() - t0, 1)},
        "method": "argument principle; the winding of arg F on |z-w|=R counts "
                  "zeros minus poles. Invariant under F -> gF for non-vanishing "
                  "holomorphic g, hence free of the normalisation dependence "
                  "that afflicts |a1/a2| and |dF/domega|.",
        "certificate": "count = 1 at radius R implies no second root within R "
                       "of the tracked root.",
        "radii_tested": args.radii,
        "n_points": len(records),
        "n_certified": len(ok),
        "min_isolation_radius": min(ok) if ok else None,
        "records": records,
    }
    pathlib.Path(args.out).write_text(json.dumps(result, indent=1) + "\n")
    print(f"\ncertified {len(ok)}/{len(records)} points; "
          f"min isolation radius = {result['min_isolation_radius']}")
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

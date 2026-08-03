"""Build the fixed-sector branch atlas by nested continuation.

Structure of the walk (never a blind uniform cube):

    stage 1   continue in mu   at (s=0, delta=0), seeded from the static table
    stage 2   continue in s    at delta=0,        seeded from stage 1
    stage 3   continue in delta  outward from 0,  seeded from stage 2

Every point is therefore reached by a chain of small steps from a mode whose
identity is known at the static endpoint, which is what makes the branch label
meaningful.  Nothing in the atlas is labelled by sorting frequencies.

Sharding is by sector: each sector writes its own immutable shard, so the run is
resumable and shards can be produced in parallel without coordination.

Domain guard: points outside the validated region (inner horizon ``z_minus``
above ``--z-minus-max``, or extremality below ``--extremality-min``) are skipped
and counted, not silently computed.  The validated domain is a property of the
solvers (``docs/EXECUTION_GATES.md`` gate 4), and quietly stepping outside it is
how the previously refuted near-extremal candidate arose.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.continuation.branch import BranchLabel, continue_along  # noqa: E402
from mp5d.geometry import MPGeometry, sdelta_to_ab  # noqa: E402

# Independent sectors.  (0,1), (0,2), (1,2) are the exchange images of (1,0),
# (2,0), (2,1) and are spot-checked rather than recomputed.
SECTORS = [(0, 0), (1, 1), (-1, -1), (1, 0), (2, 0), (2, 1), (2, 2)]

MU_GRID = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8]
# Step size in s is the binding constraint on branch identity, at BOTH ends:
#   * the first step off s=0 has only one history point, so the predictor is
#     trivial (l=6 N=3 was captured by N=2 on a 0 -> 0.12 first step);
#   * at large s the branches move fast, and a 0.27 -> 0.34 step let l=4 N=3
#     jump (Im omega reversing from -2.530 to -2.762).
# Steps of <= 0.05 hold every branch tested.  Verified resolution-independent:
# an 11-point (<=0.05) and a 12-point (<=0.04) grid give identical endpoints.
S_GRID = [0.0, 0.04, 0.08, 0.12, 0.16, 0.20, 0.25, 0.29, 0.33, 0.37, 0.42]
DELTA_POS = [0.05, 0.10, 0.15, 0.20]


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def in_domain(s: float, delta: float, z_minus_max: float,
              extremality_min: float, M: float = 1.0) -> bool:
    a, b = sdelta_to_ab(s, delta)
    geo = MPGeometry(a=a, b=b, M=M)
    if not geo.has_horizon:
        return False
    return (geo.z_minus <= z_minus_max) and (geo.extremality >= extremality_min)


def build_sector(m1: int, m2: int, seed_table: dict, args, commit: str) -> dict:
    lmin = abs(m1) + abs(m2)
    points: list[dict] = []
    failures: list[dict] = []
    skipped = 0

    for k in range(args.n_ell):
        ell = lmin + 2 * k
        seeds = seed_table.get(str(ell), [])
        for N in range(min(args.overtones, len(seeds))):
            seed = complex(seeds[N][0], seeds[N][1])
            label = BranchLabel(m1=m1, m2=m2, ell=ell, overtone=N)

            # ---- stage 1: mu at (s=0, delta=0) --------------------------
            path1 = [(0.0, 0.0, mu) for mu in MU_GRID]
            p1, f1 = continue_along(label, path1, seed, commit=commit,
                                    depth=args.depth, angular_N=args.angular_N,
                                    parent="static")
            points += [p.to_dict() for p in p1]
            failures += f1
            by_mu = {p.mu: p.omega for p in p1}

            for mu, w_mu in by_mu.items():
                # ---- stage 2: s at delta=0 ------------------------------
                path2 = [(s, 0.0, mu) for s in S_GRID
                         if in_domain(s, 0.0, args.z_minus_max,
                                      args.extremality_min)]
                if not path2:
                    continue
                p2, f2 = continue_along(label, path2, w_mu, commit=commit,
                                        depth=args.depth,
                                        angular_N=args.angular_N,
                                        parent=f"spine_mu{mu}")
                points += [p.to_dict() for p in p2]
                failures += f2

                for pt in p2:
                    s = pt.s
                    if s == 0.0:
                        continue  # delta at s=0 gives a=-b, covered by |delta| walk
                    # ---- stage 3: delta outward, both signs -------------
                    for sign in (+1, -1):
                        path3 = []
                        for d in DELTA_POS:
                            dd = sign * d
                            if in_domain(s, dd, args.z_minus_max,
                                         args.extremality_min):
                                path3.append((s, dd, mu))
                            else:
                                break
                        if not path3:
                            skipped += 1
                            continue
                        p3, f3 = continue_along(
                            label, path3, pt.omega, commit=commit,
                            depth=args.depth, angular_N=args.angular_N,
                            parent=f"s{s}_mu{mu}")
                        points += [p.to_dict() for p in p3]
                        failures += f3

    return {"sector": [m1, m2], "n_points": len(points),
            "n_failures": len(failures), "n_skipped_domain": skipped,
            "points": points, "failures": failures}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", default=None,
                    help='e.g. "1,1"; default builds every sector')
    ap.add_argument("--checkpoint-dir", default="data/full_branch_atlas")
    ap.add_argument("--seed-table", default="results/static_seed_table.json")
    ap.add_argument("--resume", action="store_true",
                    help="skip sectors whose shard already exists")
    ap.add_argument("--n-ell", type=int, default=3)
    ap.add_argument("--overtones", type=int, default=4)
    ap.add_argument("--depth", type=int, default=200)
    ap.add_argument("--angular-N", type=int, default=40)
    ap.add_argument("--z-minus-max", type=float, default=0.20)
    ap.add_argument("--extremality-min", type=float, default=0.05)
    args = ap.parse_args()

    seed_table = json.loads(pathlib.Path(args.seed_table).read_text())["table"]
    ckpt = pathlib.Path(args.checkpoint_dir)
    ckpt.mkdir(parents=True, exist_ok=True)
    commit = commit_hash()

    sectors = SECTORS
    if args.sector:
        m1, m2 = (int(x) for x in args.sector.split(","))
        sectors = [(m1, m2)]

    for (m1, m2) in sectors:
        shard = ckpt / f"sector_{m1}_{m2}.json"
        if args.resume and shard.exists():
            print(f"[skip] {shard} exists", flush=True)
            continue
        t0 = time.time()
        out = build_sector(m1, m2, seed_table, args, commit)
        out["provenance"] = {
            "commit": commit, "solver": "leaver-cf", "depth": args.depth,
            "angular_N": args.angular_N, "precision": "double",
            "z_minus_max": args.z_minus_max,
            "extremality_min": args.extremality_min,
            "mu_grid": MU_GRID, "s_grid": S_GRID, "delta_pos": DELTA_POS,
            "runtime_s": round(time.time() - t0, 1),
        }
        tmp = shard.with_suffix(".tmp")
        tmp.write_text(json.dumps(out, indent=1) + "\n")
        tmp.replace(shard)   # atomic
        print(f"[done] sector ({m1},{m2}): {out['n_points']} pts, "
              f"{out['n_failures']} failures, {time.time()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

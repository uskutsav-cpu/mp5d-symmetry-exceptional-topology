"""Adaptive collision search over the branch atlas.

Two diagnostics, because they fail in different directions:

**Tier 1 -- pairwise branch gap.**  At each atlas parameter point, the minimum
``|omega_i - omega_j|`` over distinct branches *within the same sector*.  Free,
given the atlas.  An EP2 requires two branches to coalesce, so a lower bound on
this gap over a box is a genuine exclusion.  Its weakness: it only sees branches
that were actually computed.

**Tier 2 -- scale-free root separation.**  ``|a1/a2|`` from a Cauchy expansion of
the spectral condition about the converged root (see
``mp5d.exceptional.diagnostics``).  This sees *any* nearby zero, including
branches absent from the atlas, and it is invariant under rescaling of the
spectral condition -- unlike ``|dF/domega|``, which can be made arbitrarily small
by rescaling and therefore cannot support a bound at all.

Tier 2 runs on every point whose tier-1 gap is small, plus a deterministic
random subsample of the rest, which bounds the risk that an unmodelled partner
branch was missed everywhere.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.exceptional.diagnostics import root_separation  # noqa: E402
from mp5d.exceptional.ep_solver import spectral_condition  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def load_atlas(d: pathlib.Path) -> list[dict]:
    pts: list[dict] = []
    for f in sorted(d.glob("sector_*.json")):
        pts += json.loads(f.read_text())["points"]
    return pts


def tier1_gaps(points: list[dict]) -> list[dict]:
    """Minimum within-sector, same-parameter branch gap at every atlas point."""
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for p in points:
        key = (p["m1"], p["m2"], round(p["s"], 9), round(p["delta"], 9),
               round(p["mu"], 9))
        groups[key].append(p)

    out = []
    for key, members in groups.items():
        if len(members) < 2:
            continue
        best = None
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                pi, pj = members[i], members[j]
                if pi["branch"] == pj["branch"]:
                    continue
                gap = abs(complex(pi["omega_re"], pi["omega_im"])
                          - complex(pj["omega_re"], pj["omega_im"]))
                if best is None or gap < best["gap"]:
                    best = {"gap": gap, "branch_a": pi["branch"],
                            "branch_b": pj["branch"],
                            "omega_a": [pi["omega_re"], pi["omega_im"]],
                            "omega_b": [pj["omega_re"], pj["omega_im"]],
                            "ell_a": pi["ell"], "ell_b": pj["ell"],
                            "N_a": pi["overtone"], "N_b": pj["overtone"]}
        if best is None:
            continue
        m1, m2, s, delta, mu = key
        best.update({"m1": m1, "m2": m2, "s": s, "delta": delta, "mu": mu,
                     "z_minus": members[0]["z_minus"],
                     "extremality": members[0]["extremality"]})
        out.append(best)
    out.sort(key=lambda r: r["gap"])
    return out


def tier2(rec: dict, radius_factor: float = 0.15, n: int = 32) -> dict:
    """Scale-free separation at one candidate, via a Cauchy expansion."""
    w = complex(*rec["omega_a"])
    ell = rec["ell_a"]

    def F(z: complex) -> complex:
        return spectral_condition(z, rec["s"], rec["delta"], rec["mu"],
                                  rec["m1"], rec["m2"], ell, 0, depth=200)

    radius = max(radius_factor * rec["gap"], 1e-4)
    try:
        out = root_separation(F, w, radius=radius, n=n)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": str(exc)[:200]}
    return {"ok": True, "separation": out["separation"],
            "residual_relative": out["residual_relative"],
            "radius": out["radius"]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas-dir", default="data/full_branch_atlas")
    ap.add_argument("--out", default="results/all_collision_candidates.json")
    ap.add_argument("--negative-out", default="results/negative_regions.json")
    ap.add_argument("--top", type=int, default=40,
                    help="tier-2 on this many smallest-gap points")
    ap.add_argument("--subsample", type=int, default=40,
                    help="tier-2 on this many additional random points")
    ap.add_argument("--seed", type=int, default=20260802)
    args = ap.parse_args()

    t0 = time.time()
    points = load_atlas(pathlib.Path(args.atlas_dir))
    print(f"loaded {len(points)} atlas points")
    gaps = tier1_gaps(points)
    print(f"{len(gaps)} parameter points with >=2 branches")
    if not gaps:
        print("no multi-branch points; nothing to do")
        return 1

    chosen = list(range(min(args.top, len(gaps))))
    rng = random.Random(args.seed)
    rest = list(range(len(chosen), len(gaps)))
    rng.shuffle(rest)
    chosen += rest[:args.subsample]

    for idx in chosen:
        gaps[idx]["tier2"] = tier2(gaps[idx])
        print(f"  tier2 [{idx}] gap={gaps[idx]['gap']:.4e} "
              f"sep={gaps[idx]['tier2'].get('separation')}", flush=True)

    # ---- branch-label collapse guard ---------------------------------------
    # A tier-1 gap at machine precision together with a tier-2 separation that
    # is orders of magnitude LARGER means the two labelled branches are the same
    # mode: continuation lost one of them.  That is a mode-label error, not a
    # degeneracy, and it must never be counted as an EP candidate.  Without this
    # guard the tier-1 minimum reads ~1e-15 and looks like a perfect coalescence.
    for g in gaps:
        t2 = g.get("tier2")
        g["collapsed_label"] = bool(
            t2 and t2.get("ok") and g["gap"] < 1e-10
            and t2["separation"] > 1e3 * max(g["gap"], 1e-300)
        )
    n_collapsed = sum(1 for g in gaps if g.get("collapsed_label"))

    genuine = [g for g in gaps if not g.get("collapsed_label")]
    # untested small gaps are suspect too: treat any sub-1e-10 tier-1 gap that
    # was not tier-2 tested as unverified rather than as a candidate
    genuine_tested = [g for g in genuine
                      if g["gap"] > 1e-10 or g.get("tier2", {}).get("ok")]

    seps = [g["tier2"]["separation"] for g in gaps
            if g.get("tier2", {}).get("ok") and not g.get("collapsed_label")]
    result = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "runtime_s": round(time.time() - t0, 1),
                       "solver": "leaver-cf", "precision": "double"},
        "n_atlas_points": len(points),
        "n_parameter_points": len(gaps),
        "n_tier2": len(seps),
        "n_collapsed_labels": n_collapsed,
        "collapse_note": "gap < 1e-10 with tier-2 separation > 1000x the gap means "
                         "two labelled branches are the SAME mode (continuation lost "
                         "one). Excluded from candidacy: it is a mode-label error, "
                         "not a degeneracy.",
        "min_branch_gap_raw": gaps[0]["gap"],
        "min_branch_gap": genuine_tested[0]["gap"] if genuine_tested else None,
        "min_root_separation": min(seps) if seps else None,
        "candidates": genuine_tested[:200],
        # Every tier-2 evaluation, including the random subsample, whose ranks
        # generally fall outside the top 200.  Without this the quoted
        # min_root_separation would not be recoverable from the saved records,
        # and a figure drawn from `candidates` alone would disagree with the
        # reported bound.
        "tier2_evaluated": [g for g in gaps
                            if g.get("tier2", {}).get("ok")
                            and not g.get("collapsed_label")],
        "collapsed_examples": [g for g in gaps if g.get("collapsed_label")][:10],
    }
    pathlib.Path(args.out).write_text(json.dumps(result, indent=1) + "\n")

    neg = {
        "status": "complete",
        "description": "Bounded exclusion of EP2 over the searched atlas domain. "
                       "The bound is on the scale-free root separation |a1/a2| and "
                       "on the pairwise branch gap; it is NOT based on |dF/domega|, "
                       "which is not invariant under rescaling of the spectral "
                       "condition and therefore cannot bound anything.",
        "provenance": result["provenance"],
        "domain": {
            "sectors": sorted({(g["m1"], g["m2"]) for g in gaps}.__iter__(),
                              key=lambda t: (t[0], t[1])),
            "z_minus_max": 0.20, "extremality_min": 0.05,
            "s_range": [min(g["s"] for g in gaps), max(g["s"] for g in gaps)],
            "delta_range": [min(g["delta"] for g in gaps),
                            max(g["delta"] for g in gaps)],
            "mu_range": [min(g["mu"] for g in gaps), max(g["mu"] for g in gaps)],
        },
        "min_branch_gap": genuine_tested[0]["gap"] if genuine_tested else None,
        "min_root_separation": min(seps) if seps else None,
        "n_atlas_points": len(points),
        "n_parameter_points": len(gaps),
        "n_collapsed_labels_excluded": n_collapsed,
        "verdict": None,  # filled in by the reporting step after inspection
    }
    pathlib.Path(args.negative_out).write_text(json.dumps(neg, indent=1) + "\n")
    print(f"collapsed labels excluded = {n_collapsed}")
    print(f"min branch gap        = {result['min_branch_gap']}")
    print(f"min root separation   = {result['min_root_separation']}")
    print(f"wrote {args.out} and {args.negative_out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Apply the EP2 verification gate to the ranked collision candidates.

Implements `docs/EP_VERIFICATION_GATE.md`.  Every candidate is either verified
against all applicable conditions or **rejected with the condition that failed
named**.  Nothing is left as "promising".

Gates run here (G1-G6); G8-G10 are solver/domain conditions handled elsewhere:

  G1  scale-free root separation |a1/a2|, at two contour radii
  G2  argument-principle zero count inside a small contour
  G3  both roots located by contour moments, without a starting guess
  G4  monodromy: one loop in a parameter 2-plane transposes the branches
  G6  Puiseux: splitting ~ t^{1/2} preferred over t^1

A candidate that is really an avoided crossing fails G1, G2 and G4 together.
That is the expected outcome and it is recorded as evidence, not as a
non-result.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.exceptional.diagnostics import root_separation, winding_zero_count  # noqa: E402
from mp5d.exceptional.ep_solver import spectral_condition  # noqa: E402
from mp5d.exceptional.verification import (  # noqa: E402
    fit_puiseux,
    monodromy,
    roots_in_disc,
)


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def make_F(c: dict, depth: int = 200):
    """Spectral condition for this candidate, inverted at its own overtone.

    Using inversion 0 for a branch with N > 0 yields a function that does not
    vanish at the tracked root, so every gate built on it is void.
    """
    ell = c["ell_a"]
    inv = c["N_a"]

    def F(z: complex, s=None, delta=None, mu=None) -> complex:
        return spectral_condition(
            z,
            c["s"] if s is None else s,
            c["delta"] if delta is None else delta,
            c["mu"] if mu is None else mu,
            c["m1"], c["m2"], ell, inv, depth=depth)

    return F


def verify(c: dict, depth: int, monodromy_points: int) -> dict:
    F = make_F(c, depth)
    w = complex(*c["omega_a"])
    gap = c["gap"]
    out: dict = {"candidate": {k: c[k] for k in
                              ("m1", "m2", "s", "delta", "mu", "gap",
                               "branch_a", "branch_b", "ell_a", "ell_b",
                               "N_a", "N_b")},
                 "gates": {}}

    # ---- G1: scale-free separation at two radii -------------------------
    try:
        r1 = max(0.2 * gap, 1e-4)
        r2 = max(0.4 * gap, 2e-4)
        s1 = root_separation(F, w, radius=r1, n=32)
        s2 = root_separation(F, w, radius=r2, n=32)
        out["gates"]["G1"] = {
            "separation_r1": s1["separation"], "separation_r2": s2["separation"],
            "radius_r1": r1, "radius_r2": r2,
            "residual_relative": s1["residual_relative"],
            "pass": bool(min(s1["separation"], s2["separation"]) < 1e-3),
        }
    except Exception as exc:  # noqa: BLE001
        out["gates"]["G1"] = {"error": str(exc)[:200], "pass": False}

    # ---- G2: argument-principle zero count ------------------------------
    try:
        rad = max(2.0 * gap, 1e-3)
        wind = winding_zero_count(F, w, rad, n=192)
        out["gates"]["G2"] = {"winding": wind, "radius": rad,
                              "pass": bool(abs(wind - 2.0) < 0.15)}
    except Exception as exc:  # noqa: BLE001
        out["gates"]["G2"] = {"error": str(exc)[:200], "pass": False}

    # ---- G3: both roots, no starting guess ------------------------------
    try:
        rad = max(2.0 * gap, 1e-3)
        rts = roots_in_disc(F, w, rad, n=192, max_roots=3)
        out["gates"]["G3"] = {"n_roots": len(rts),
                              "roots": [[z.real, z.imag] for z in rts],
                              "pass": bool(len(rts) >= 2)}
    except Exception as exc:  # noqa: BLE001
        out["gates"]["G3"] = {"error": str(exc)[:200], "pass": False}

    # ---- G4: monodromy in the (delta, mu) plane -------------------------
    try:
        rad_p = 0.02

        def root_fn(p):
            def G(z):
                return F(z, s=c["s"], delta=p[0], mu=p[1])
            return roots_in_disc(G, w, max(4.0 * gap, 5e-3), n=128, max_roots=3)

        mono = monodromy(root_fn, [c["delta"], c["mu"]], rad_p,
                         n_points=monodromy_points, n_roots=2, plane=(0, 1))
        out["gates"]["G4"] = {"permutation": mono.permutation,
                              "is_swap": mono.is_swap,
                              "closed_error": mono.closed_error,
                              "max_step": mono.max_step,
                              "pass": bool(mono.is_swap)}
    except Exception as exc:  # noqa: BLE001
        out["gates"]["G4"] = {"error": str(exc)[:200], "pass": False}

    # ---- G6: Puiseux exponent along mu ----------------------------------
    try:
        ts, split = [], []
        for dt in np.geomspace(1e-4, 2e-2, 7):
            def G(z, dd=dt):
                return F(z, s=c["s"], delta=c["delta"], mu=c["mu"] + dd)
            rts = roots_in_disc(G, w, max(6.0 * gap, 8e-3), n=128, max_roots=3)
            if len(rts) >= 2:
                rts = sorted(rts, key=lambda z: abs(z - w))[:2]
                ts.append(float(dt))
                split.append(rts[0] - rts[1])
        if len(ts) >= 3:
            fit = fit_puiseux(ts, split)
            out["gates"]["G6"] = {"exponent": fit.exponent,
                                  "residual_sqrt": fit.residual_sqrt,
                                  "residual_linear": fit.residual_linear,
                                  "prefers_sqrt": fit.prefers_sqrt,
                                  "n_points": fit.n_points,
                                  "pass": bool(fit.prefers_sqrt
                                               and abs(fit.exponent - 0.5) < 0.15)}
        else:
            out["gates"]["G6"] = {"error": "fewer than 3 usable samples",
                                  "pass": False}
    except Exception as exc:  # noqa: BLE001
        out["gates"]["G6"] = {"error": str(exc)[:200], "pass": False}

    passed = [g for g, v in out["gates"].items() if v.get("pass")]
    failed = [g for g, v in out["gates"].items() if not v.get("pass")]
    out["gates_passed"] = passed
    out["gates_failed"] = failed

    # G2 and G3 are PRECONDITIONS -- they establish only that two distinct roots
    # sit near each other, which is what makes something a candidate in the
    # first place.  The EP-specific evidence is G1 (coalescence), G4 (monodromy)
    # and G6 (square-root splitting).  Counting G2/G3 toward an EP verdict would
    # let every ordinary crossing look half-verified.
    decisive = ("G1", "G4", "G6")
    dec_pass = [g for g in decisive if out["gates"].get(g, {}).get("pass")]
    dec_fail = [g for g in decisive if not out["gates"].get(g, {}).get("pass")]
    preconditions_met = all(out["gates"].get(g, {}).get("pass") for g in ("G2", "G3"))

    if not dec_fail:
        verdict, classification = "VERIFIED_EP2", "exceptional point of order 2"
    elif not dec_pass:
        verdict = "REJECTED"
        classification = ("avoided or ordinary crossing: two distinct roots "
                          "present (G2, G3) but no coalescence, no monodromy "
                          "exchange and no square-root splitting"
                          if preconditions_met else
                          "not a coalescence; preconditions not met either")
    else:
        verdict, classification = "UNRESOLVED", "mixed decisive evidence"

    out["decisive_gates_passed"] = dec_pass
    out["decisive_gates_failed"] = dec_fail
    out["verdict"] = verdict
    out["classification"] = classification
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="results/all_collision_candidates.json")
    ap.add_argument("--out", default="results/verified_exceptional_points.json")
    ap.add_argument("--rejected-out", default="results/all_rejected_candidates.json")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--depth", type=int, default=200)
    ap.add_argument("--monodromy-points", type=int, default=32)
    args = ap.parse_args()

    t0 = time.time()
    data = json.loads(pathlib.Path(args.candidates).read_text())
    cands = data["candidates"][:args.top]
    print(f"verifying {len(cands)} candidates")

    results = []
    for i, c in enumerate(cands):
        r = verify(c, args.depth, args.monodromy_points)
        results.append(r)
        print(f"  [{i}] ({c['m1']},{c['m2']}) gap={c['gap']:.3e} -> {r['verdict']} "
              f"passed={r['gates_passed']} failed={r['gates_failed']}", flush=True)

    verified = [r for r in results if r["verdict"] == "VERIFIED_EP2"]
    rejected = [r for r in results if r["verdict"] == "REJECTED"]
    unresolved = [r for r in results if r["verdict"] == "UNRESOLVED"]

    prov = {"commit": commit_hash(), "solver": "leaver-cf", "precision": "double",
            "depth": args.depth, "runtime_s": round(time.time() - t0, 1)}

    pathlib.Path(args.out).write_text(json.dumps({
        "status": "complete",
        "provenance": prov,
        "gate": "docs/EP_VERIFICATION_GATE.md (G1-G6 here; G8-G10 elsewhere)",
        "n_verified": len(verified),
        "n_rejected": len(rejected),
        "n_unresolved": len(unresolved),
        "verified": verified,
    }, indent=1) + "\n")

    pathlib.Path(args.rejected_out).write_text(json.dumps({
        "status": "complete",
        "provenance": prov,
        "note": "Rejections are evidence. Each entry names the gate that failed.",
        "rejected": rejected,
        "unresolved": unresolved,
    }, indent=1) + "\n")

    print(f"verified={len(verified)} rejected={len(rejected)} "
          f"unresolved={len(unresolved)} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

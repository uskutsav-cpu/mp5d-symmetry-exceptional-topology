"""Track massive-scalar branches into the long-lived / quasi-bound regime.

Physics being resolved
----------------------
As the scalar mass ``mu`` grows at fixed background, a branch migrates toward
the threshold ``omega^2 = mu^2``.  Past it, ``Omega = sqrt(omega^2 - mu^2)``
turns imaginary and the outgoing condition at infinity becomes a *decaying*
condition -- the mode is quasi-bound, and the damping ``-Im omega`` collapses by
orders of magnitude.  This is the long-lived behaviour reported by Huang-Huang.

Why a fixed step fails
----------------------
Near the transition ``domega/dmu`` becomes large, so a uniform ``mu`` grid steps
straight over it and the corrector lands on a different branch.  That is exactly
what the fixed-grid atlas run reported as ``predictor_miss`` at ``mu = 1.8``.
Here the step is halved whenever the predictor misses or the frequency moves too
far, down to ``--min-step``, and the walk stops -- with the reason recorded --
rather than accepting a jump.

Classification recorded per point:

* ``ordinary``    -- ``-Im omega > 0.05``
* ``long_lived``  -- ``-Im omega <= 0.05``
* ``quasi_bound`` -- ``Re omega < mu`` (below threshold)
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.continuation.branch import extrapolate  # noqa: E402
from mp5d.geometry import sdelta_to_ab  # noqa: E402
from mp5d.radial.qnm import solve_qnm_cf  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def classify(omega: complex, mu: float) -> str:
    if omega.real < mu:
        return "quasi_bound"
    if -omega.imag <= 0.05:
        return "long_lived"
    return "ordinary"


def walk_mu(m1: int, m2: int, ell: int, N: int, s: float, delta: float,
            seed: complex, mu_max: float, step0: float, min_step: float,
            depth: int, angular_N: int) -> tuple[list[dict], dict]:
    a, b = sdelta_to_ab(s, delta)
    pts: list[dict] = []
    hist: list[complex] = []
    mus: list[float] = []
    mu = 0.0
    step = step0
    stop = {"reason": "reached_mu_max"}

    while mu <= mu_max + 1e-12:
        guess = seed if not hist else extrapolate(hist, mus, mu)
        try:
            sol = solve_qnm_cf(a, b, mu, m1, m2, ell, N,
                               initial_frequency=guess, depth=depth,
                               angular_N=angular_N, depth_schedule=(depth // 2, depth))
        except Exception as exc:  # noqa: BLE001
            sol = None
            err = str(exc)[:150]

        ok = sol is not None and sol.converged and sol.cf_residual < 1e-8
        if ok:
            w = sol.omega
            moved = abs(w - guess)
            # reject an implausible jump; a real branch moves smoothly in mu
            if hist and moved > 0.15 * max(abs(w), 1.0):
                ok = False

        if not ok:
            if step > min_step:
                step *= 0.5
                mu = (mus[-1] + step) if mus else step
                continue
            stop = {"reason": "stalled" if sol is not None else "solver_error",
                    "mu": mu, "min_step": min_step}
            if sol is None:
                stop["detail"] = err
            break

        w = sol.omega
        pts.append({
            "mu": mu, "s": s, "delta": delta, "a": a, "b": b,
            "m1": m1, "m2": m2, "ell": ell, "overtone": N,
            "omega": [w.real, w.imag],
            "Lambda": [sol.Lambda.real, sol.Lambda.imag],
            "damping": -w.imag,
            "threshold_distance": abs(w.real) - mu,
            "regime": classify(w, mu),
            "residual": sol.cf_residual,
            "depth": depth, "angular_N": angular_N,
        })
        hist.append(w)
        mus.append(mu)
        if mu >= mu_max - 1e-12:
            break
        # grow the step back once things are smooth again
        step = min(step * 1.3, step0)
        mu = min(mu + step, mu_max)

    return pts, stop


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/long_lived_branches_final.json")
    ap.add_argument("--mu-max", type=float, default=3.0)
    ap.add_argument("--step", type=float, default=0.15)
    ap.add_argument("--min-step", type=float, default=0.002)
    ap.add_argument("--depth", type=int, default=200)
    ap.add_argument("--angular-N", type=int, default=40)
    ap.add_argument("--seed-table", default="results/static_seed_table.json")
    args = ap.parse_args()

    seeds = json.loads(pathlib.Path(args.seed_table).read_text())["table"]
    t0 = time.time()
    branches = []

    configs = [
        # (m1, m2, ell, N, s, delta)
        (1, 1, 2, 0, 0.0, 0.0), (1, 1, 2, 0, 0.24, 0.0), (1, 1, 2, 0, 0.24, 0.10),
        (1, 1, 2, 1, 0.24, 0.0),
        (0, 0, 0, 0, 0.0, 0.0), (0, 0, 0, 0, 0.24, 0.0),
        (0, 0, 2, 0, 0.24, 0.0),
        (2, 2, 4, 0, 0.24, 0.0),
        (1, 0, 1, 0, 0.24, 0.0),
        (2, 0, 2, 0, 0.24, 0.0),
        (2, 1, 3, 0, 0.24, 0.0),
    ]
    for (m1, m2, ell, N, s, delta) in configs:
        tbl = seeds.get(str(ell), [])
        if N >= len(tbl):
            continue
        seed = complex(tbl[N][0], tbl[N][1])
        pts, stop = walk_mu(m1, m2, ell, N, s, delta, seed, args.mu_max,
                            args.step, args.min_step, args.depth, args.angular_N)
        regimes = {p["regime"] for p in pts}
        lowest = min((p["damping"] for p in pts), default=None)
        branches.append({
            "labels": {"m1": m1, "m2": m2, "ell": ell, "overtone": N},
            "background": {"s": s, "delta": delta},
            "n_points": len(pts), "stop": stop,
            "regimes_visited": sorted(regimes),
            "lowest_damping": lowest,
            "mu_at_lowest_damping": min(pts, key=lambda p: p["damping"])["mu"]
            if pts else None,
            "points": pts,
        })
        print(f"({m1},{m2}) l={ell} N={N} s={s} d={delta}: {len(pts)} pts, "
              f"regimes={sorted(regimes)}, min damping={lowest}, stop={stop['reason']}",
              flush=True)

    out = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "solver": "leaver-cf",
                       "precision": "double", "runtime_s": round(time.time() - t0, 1)},
        "classification": {
            "ordinary": "-Im omega > 0.05",
            "long_lived": "-Im omega <= 0.05",
            "quasi_bound": "Re omega < mu (below the omega^2=mu^2 threshold)",
        },
        "lowest_damping_overall": min(
            (b["lowest_damping"] for b in branches if b["lowest_damping"] is not None),
            default=None),
        "branches": branches,
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

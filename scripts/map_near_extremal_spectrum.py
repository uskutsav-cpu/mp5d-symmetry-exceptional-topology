"""Classify the near-extremal difficulty: spectral, or merely numerical?

Standing question (`docs/MULTIDOMAIN_NEAR_HORIZON.md`): near extremality no
solver converges, and two hypotheses -- unresolved near-horizon scale, and
double-precision conditioning -- were each measured and **refuted**.  The
remaining possibility was that the difficulty is *spectral*: the object being
tracked may not be an isolated resonance at all.

Discriminating test (exterior complex scaling)
----------------------------------------------
Under exterior complex scaling by angle ``theta``:

* a **genuine resonance** is an eigenvalue of the scaled operator and does
  **not** move as ``theta`` varies -- it is exposed by the rotation, not created
  by it;
* the **rotated continuum** sweeps through the lower half plane *with* ``theta``;
* a **branch point** or accumulation locus shows up as an arc or cluster of
  small singular values rather than an isolated minimum.

So: map ``sigma_min`` over a complex-frequency window at several ``theta`` and
ask what stays put.  This is the standard resonance/continuum discriminator and
it needs no root finder, so it is immune to the Newton stagnation that produced
the refuted candidate.

Output feeds `docs/NEAR_EXTREMAL_SPECTRAL_CLASSIFICATION.md`.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.geometry import MPGeometry  # noqa: E402
from mp5d.radial.solver_c import SolverCProblem, solve_qnm_c  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def ab_for_z_minus(target: float, ratio: float = 1.6) -> tuple[float, float]:
    """Find ``(a, b)`` with ``b = a/ratio`` giving the requested inner horizon.

    ``z_+ z_- = a^2 b^2`` and ``z_+ + z_- = M - a^2 - b^2``, so ``z_-`` grows
    monotonically as the spins grow toward extremality; a bisection is enough.
    """
    lo, hi = 0.01, 0.999
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        g = MPGeometry(a=mid, b=mid / ratio, M=1.0)
        if not g.has_horizon:
            hi = mid
            continue
        if g.z_minus < target:
            lo = mid
        else:
            hi = mid
    a = 0.5 * (lo + hi)
    return a, a / ratio


def map_window(prob: SolverCProblem, centre: complex, half_re: float,
               half_im: float, n_re: int, n_im: int) -> dict:
    res = np.linspace(centre.real - half_re, centre.real + half_re, n_re)
    ims = np.linspace(centre.imag - half_im, centre.imag + half_im, n_im)
    grid = np.zeros((n_im, n_re))
    for i, im in enumerate(ims):
        for j, re in enumerate(res):
            try:
                s, _ = prob.sigma_min(complex(re, im))
            except Exception:  # noqa: BLE001
                s = float("nan")
            grid[i, j] = s
    k = int(np.nanargmin(grid))
    i0, j0 = divmod(k, n_re)
    return {
        "re": res.tolist(), "im": ims.tolist(),
        "log10_sigma_min": np.log10(np.maximum(grid, 1e-300)).round(4).tolist(),
        "min_sigma": float(np.nanmin(grid)),
        "argmin": [float(res[j0]), float(ims[i0])],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/near_extremal_spectral_map.json")
    ap.add_argument("--z-minus", type=float, nargs="+",
                    default=[0.20, 0.30, 0.40, 0.44])
    ap.add_argument("--thetas-deg", type=float, nargs="+", default=[50.0, 60.0, 70.0])
    ap.add_argument("--resolution", type=int, default=180)
    ap.add_argument("--n-re", type=int, default=25)
    ap.add_argument("--n-im", type=int, default=25)
    ap.add_argument("--mu", type=float, default=1.9)
    ap.add_argument("--ell", type=int, default=2)
    ap.add_argument("--m1", type=int, default=1)
    ap.add_argument("--m2", type=int, default=1)
    ap.add_argument("--centre", type=float, nargs=2, default=[2.14, -0.49],
                help="seed for the LOWEST z_minus; later windows are re-centred "
                     "by continuation")
    ap.add_argument("--half-re", type=float, default=0.45)
    ap.add_argument("--half-im", type=float, default=0.45)
    args = ap.parse_args()

    t0 = time.time()
    entries = []
    # Centre each window on the actual mode at that background, obtained by
    # continuation in z_minus from a converged low-z_minus mode.  Using one
    # fixed window for every z_minus put the minimum on the window edge and made
    # the theta comparison meaningless -- the drift measured the boundary, not
    # the spectrum.
    centre = complex(*args.centre)
    for z in sorted(args.z_minus):
        a, b = ab_for_z_minus(z)
        geo = MPGeometry(a=a, b=b, M=1.0)
        try:
            seeded = solve_qnm_c(a, b, args.mu, args.m1, args.m2, args.ell,
                                 initial_frequency=centre,
                                 theta=math.radians(args.thetas_deg[0]),
                                 L=60.0, resolution=args.resolution)
            if seeded.converged:
                centre = seeded.omega
        except Exception as exc:  # noqa: BLE001
            print(f"  [warn] seeding failed at z-={z}: {str(exc)[:80]}", flush=True)
        print(f"z-={z:.3f}: window centred on {centre:.6f}", flush=True)
        per_theta = {}
        for th_deg in args.thetas_deg:
            th = math.radians(th_deg)
            prob = SolverCProblem(geo, args.mu, args.m1, args.m2, args.ell,
                                  theta=th, L=60.0, resolution=args.resolution)
            m = map_window(prob, centre, args.half_re, args.half_im,
                           args.n_re, args.n_im)
            per_theta[f"{th_deg:g}"] = m
            print(f"z-={z:.3f} (a={a:.4f},b={b:.4f}) theta={th_deg:g}deg "
                  f"min_sigma={m['min_sigma']:.3e} at "
                  f"{m['argmin'][0]:.4f}{m['argmin'][1]:+.4f}i", flush=True)

        # how far the deepest point moves as theta varies:
        # small  -> resonance (theta-independent)
        # large  -> rotated continuum / not an isolated mode
        pts = [complex(*per_theta[k]["argmin"]) for k in per_theta]
        drift = max(abs(p - q) for p in pts for q in pts)
        entries.append({
            "z_minus": z, "a": a, "b": b,
            "extremality": geo.extremality,
            "argmin_drift_over_theta": float(drift),
            "min_sigma_by_theta": {k: per_theta[k]["min_sigma"] for k in per_theta},
            "maps": per_theta,
        })

    out = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "solver": "solver-C "
                       "(exterior complex scaling)", "precision": "double",
                       "resolution": args.resolution,
                       "runtime_s": round(time.time() - t0, 1)},
        "test": "A genuine resonance is theta-independent under exterior complex "
                "scaling; the rotated continuum moves with theta. "
                "argmin_drift_over_theta is the discriminator.",
        "parameters": {"mu": args.mu, "ell": args.ell,
                       "sector": [args.m1, args.m2], "centre": args.centre},
        "entries": [{k: v for k, v in e.items() if k != "maps"} for e in entries],
        "maps": {str(e["z_minus"]): e["maps"] for e in entries},
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Targeted pseudospectral analysis at selected spectral points.

Deliberately **not** a generic pseudospectrum survey.  The quantities below are
computed only where they discriminate: at an ordinary mode (control) and at the
strongest branch interaction found in the atlas.

Definitions actually used (stating these is mandatory -- a pseudospectrum is
meaningless without its norm and operator)
-------------------------------------------------------------------------
* **Operator**: Solver C's exterior-complex-scaled Chebyshev collocation matrix
  ``T(omega)`` for one ``(m1, m2, l)`` sector, with the angular eigenvalue
  ``Lambda(omega)`` recomputed at every ``omega`` (never frozen).
* **Norm**: the matrix 2-norm induced by the Euclidean vector norm on
  collocation coefficients.  This is a *discretization* norm, not an
  energy norm on the continuum solution space; the two differ, and no claim
  here is transferred to the continuum.
* **epsilon-pseudospectrum**: ``{omega : sigma_min(T(omega)) <= epsilon}``,
  reported as ``log10 sigma_min`` on a grid, normalized by ``sigma_max`` so the
  numbers do not drift with the arbitrary overall scale of ``T``.
* **Eigenvalue condition number** for the nonlinear pencil at a simple root:

      kappa = ||u|| ||v|| / |v^H T'(omega_*) u|

  with ``u``, ``v`` the right and left null vectors.  ``kappa -> infinity`` is
  the signature of an approaching defective coalescence, and it is the
  finite-dimensional shadow of the ``|a1/a2| -> 0`` diagnostic.
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

from mp5d.geometry import MPGeometry, sdelta_to_ab  # noqa: E402
from mp5d.radial.solver_c import SolverCProblem  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def condition_number(prob: SolverCProblem, omega: complex,
                     h: float = 1e-6) -> dict:
    T = prob.matrix(omega)
    U, S, Vh = np.linalg.svd(T)
    u = Vh[-1].conj()
    v = U[:, -1]
    dT = (prob.matrix(omega + h) - prob.matrix(omega - h)) / (2 * h)
    denom = complex(v.conj() @ (dT @ u))
    kappa = (np.linalg.norm(u) * np.linalg.norm(v) / abs(denom)) if denom != 0 else np.inf
    return {
        "sigma_min_rel": float(S[-1] / S[0]),
        "sigma_ratio_last_two": float(S[-2] / S[-1]),
        "abs_v_dT_u": float(abs(denom)),
        "condition_number": float(kappa),
    }


def perturbation_response(prob: SolverCProblem, omega: complex, eps: float,
                          n_samples: int, seed: int, structured: bool) -> dict:
    """Move ``sigma_min`` under random perturbations of size ``eps``.

    ``structured`` restricts the perturbation to the ODE-coefficient rows (a
    perturbation of the potential), whereas the unstructured case perturbs every
    entry including the boundary rows.  A large gap between the two means the
    naive unstructured pseudospectrum overstates the physical sensitivity.
    """
    rng = np.random.default_rng(seed)
    T = prob.matrix(omega)
    nrm = np.linalg.norm(T, 2)
    n = T.shape[0]
    shifts = []
    for _ in range(n_samples):
        E = rng.normal(size=T.shape) + 1j * rng.normal(size=T.shape)
        if structured:
            E[-2:, :] = 0.0          # leave the two boundary rows untouched
        E *= eps * nrm / np.linalg.norm(E, 2)
        S = np.linalg.svd(T + E, compute_uv=False)
        shifts.append(float(S[-1] / S[0]))
    base = float(np.linalg.svd(T, compute_uv=False)[-1] / np.linalg.svd(T, compute_uv=False)[0])
    return {"eps": eps, "n": n, "baseline_sigma_min_rel": base,
            "mean_sigma_min_rel": float(np.mean(shifts)),
            "max_sigma_min_rel": float(np.max(shifts)),
            "structured": structured}


def grid_map(prob: SolverCProblem, centre: complex, half: float, n: int) -> dict:
    res = np.linspace(centre.real - half, centre.real + half, n)
    ims = np.linspace(centre.imag - half, centre.imag + half, n)
    g = np.zeros((n, n))
    for i, im in enumerate(ims):
        for j, re in enumerate(res):
            try:
                s, _ = prob.sigma_min(complex(re, im))
            except Exception:  # noqa: BLE001
                s = float("nan")
            g[i, j] = s
    return {"re": res.tolist(), "im": ims.tolist(),
            "log10_sigma_min_rel": np.log10(np.maximum(g, 1e-300)).round(4).tolist(),
            "min": float(np.nanmin(g))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/pseudospectral_results.json")
    ap.add_argument("--candidates", default="results/all_collision_candidates.json")
    ap.add_argument("--resolution", type=int, default=180)
    ap.add_argument("--grid", type=int, default=21)
    ap.add_argument("--theta-deg", type=float, default=60.0)
    ap.add_argument("--seed", type=int, default=20260802)
    args = ap.parse_args()

    t0 = time.time()
    th = math.radians(args.theta_deg)
    points = []

    # -- control: an ordinary mode well inside the validated domain ---------
    control = {"label": "ordinary_control", "s": 0.24, "delta": 0.10, "mu": 0.6,
               "m1": 1, "m2": 1, "ell": 2}

    targets = [control]
    cand_path = pathlib.Path(args.candidates)
    if cand_path.exists():
        cands = json.loads(cand_path.read_text())["candidates"]
        if cands:
            c = cands[0]
            targets.append({"label": "tightest_interaction", "s": c["s"],
                            "delta": c["delta"], "mu": c["mu"], "m1": c["m1"],
                            "m2": c["m2"], "ell": c["ell_a"],
                            "omega_hint": c["omega_a"], "branch_gap": c["gap"]})

    for tgt in targets:
        a, b = sdelta_to_ab(tgt["s"], tgt["delta"])
        geo = MPGeometry(a=a, b=b, M=1.0)
        prob = SolverCProblem(geo, tgt["mu"], tgt["m1"], tgt["m2"], tgt["ell"],
                              theta=th, L=60.0, resolution=args.resolution)
        # locate the deepest point of sigma_min near the hint
        hint = complex(*tgt["omega_hint"]) if "omega_hint" in tgt else 1.6 - 0.35j
        gm = grid_map(prob, hint, 0.25, args.grid)
        k = int(np.nanargmin(np.array(gm["log10_sigma_min_rel"])))
        i0, j0 = divmod(k, args.grid)
        w = complex(gm["re"][j0], gm["im"][i0])

        entry = {
            **{k2: v for k2, v in tgt.items()},
            "a": a, "b": b,
            "omega_used": [w.real, w.imag],
            "conditioning": condition_number(prob, w),
            "perturbation_unstructured": perturbation_response(
                prob, w, 1e-8, 8, args.seed, structured=False),
            "perturbation_structured": perturbation_response(
                prob, w, 1e-8, 8, args.seed, structured=True),
            "sigma_map": gm,
        }
        points.append(entry)
        print(f"{tgt['label']}: kappa={entry['conditioning']['condition_number']:.4e} "
              f"sigma_min_rel={entry['conditioning']['sigma_min_rel']:.3e}", flush=True)

    out = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "solver": "solver-C",
                       "precision": "double", "resolution": args.resolution,
                       "theta_deg": args.theta_deg,
                       "runtime_s": round(time.time() - t0, 1)},
        "definitions": {
            "operator": "Solver C exterior-complex-scaled Chebyshev collocation "
                        "matrix; Lambda(omega) recomputed at every omega",
            "norm": "matrix 2-norm on collocation coefficients (a DISCRETIZATION "
                    "norm, not a continuum energy norm)",
            "pseudospectrum": "{omega : sigma_min(T)/sigma_max(T) <= epsilon}",
            "condition_number": "||u|| ||v|| / |v^H T'(omega) u| at a simple root",
        },
        "caveat": "All quantities are properties of the finite-dimensional "
                  "discretization. No continuum pseudospectral claim is made.",
        "points": [{k: v for k, v in p.items() if k != "sigma_map"} for p in points],
        "sigma_maps": {p["label"]: p["sigma_map"] for p in points},
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

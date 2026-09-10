"""Local two-branch effective model and its discriminant.

Any two interacting branches can be described locally by a 2x2 non-Hermitian
effective matrix

    H_eff(p) = [[E1, V], [W, E2]]

whose eigenvalues are ``(T +- sqrt(D))/2`` with ``T = E1+E2`` and
``D = (E1-E2)^2 + 4 V W``.  The individual entries are not observable -- any
similarity transformation changes them -- but the two symmetric functions are:

    T(p) = omega_+ + omega_-          (trace)
    D(p) = (omega_+ - omega_-)^2      (discriminant)

**An EP2 is exactly ``D = 0``**, and ``D`` is built from eigenvalue differences,
so unlike ``|dF/domega|`` it is invariant under any rescaling of the spectral
condition.  Fitting ``D`` is therefore both the effective theory and the
quantitative EP diagnostic, and ``min |D|`` over a box is a bound with meaning.

Symmetry constraint used
------------------------
In a diagonal sector (``m1 = m2``) every branch is exactly even in ``delta``
(``docs/SYMMETRY_GROUP.md``), so ``D`` is an analytic function of

    t = delta^2

and the fit is performed in ``(s, t, mu)``.  Imposing this rather than fitting a
general dependence on ``delta`` is what makes the model predictive with few
coefficients -- and the residual of the fit is itself a test of the symmetry.

The fit is validated **out of sample**: coefficients are determined on one
subset of atlas points and scored on held-out points.
"""

from __future__ import annotations

import argparse
import itertools
import json
import pathlib
import subprocess
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))


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


def monomials(s, t, mu, degree: int):
    """All monomials ``s^i t^j mu^k`` with ``i+j+k <= degree``."""
    cols, names = [], []
    for i, j, k in itertools.product(range(degree + 1), repeat=3):
        if i + j + k <= degree:
            cols.append((s ** i) * (t ** j) * (mu ** k))
            names.append(f"s^{i} t^{j} mu^{k}")
    return np.array(cols).T, names


def fit_pair(records: list[dict], degree: int, holdout_frac: float,
             seed: int) -> dict | None:
    """Fit ``T`` and ``D`` for one branch pair over its atlas points."""
    if len(records) < 15:
        return None
    s = np.array([r["s"] for r in records])
    t = np.array([r["delta"] for r in records]) ** 2
    mu = np.array([r["mu"] for r in records])
    wp = np.array([complex(*r["omega_a"]) for r in records])
    wm = np.array([complex(*r["omega_b"]) for r in records])

    T = wp + wm
    D = (wp - wm) ** 2

    A, names = monomials(s, t, mu, degree)
    if A.shape[0] <= A.shape[1] + 3:
        return None

    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(records))
    n_hold = max(3, int(holdout_frac * len(records)))
    hold, train = idx[:n_hold], idx[n_hold:]

    out = {}
    for name, y in (("T", T), ("D", D)):
        coef, *_ = np.linalg.lstsq(A[train], y[train], rcond=None)
        pred_h = A[hold] @ coef
        pred_t = A[train] @ coef
        scale = float(np.mean(np.abs(y))) or 1.0
        out[name] = {
            "coefficients": [[c.real, c.imag] for c in coef],
            "monomials": names,
            "train_rel_rms": float(np.sqrt(np.mean(np.abs(pred_t - y[train]) ** 2)) / scale),
            "holdout_rel_rms": float(np.sqrt(np.mean(np.abs(pred_h - y[hold]) ** 2)) / scale),
            "n_train": int(len(train)), "n_holdout": int(len(hold)),
        }
    out["min_abs_D_observed"] = float(np.min(np.abs(D)))
    out["max_abs_D_observed"] = float(np.max(np.abs(D)))
    out["argmin_D"] = {
        k: records[int(np.argmin(np.abs(D)))][k]
        for k in ("s", "delta", "mu", "m1", "m2", "ell_a", "ell_b", "N_a", "N_b")
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="results/all_collision_candidates.json")
    ap.add_argument("--out", default="results/effective_models.json")
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--holdout-frac", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=20260802)
    ap.add_argument("--max-pairs", type=int, default=12)
    args = ap.parse_args()

    data = json.loads(pathlib.Path(args.candidates).read_text())
    cands = data["candidates"]

    groups: dict[tuple, list[dict]] = defaultdict(list)
    for c in cands:
        groups[(c["m1"], c["m2"], c["branch_a"], c["branch_b"])].append(c)

    ranked = sorted(groups.items(), key=lambda kv: min(r["gap"] for r in kv[1]))
    models = []
    for key, recs in ranked[:args.max_pairs]:
        fit = fit_pair(recs, args.degree, args.holdout_frac, args.seed)
        if fit is None:
            continue
        m1, m2, ba, bb = key
        fit["pair"] = {"m1": m1, "m2": m2, "branch_a": ba, "branch_b": bb,
                       "n_points": len(recs),
                       "diagonal_sector": (m1 == m2)}
        models.append(fit)
        print(f"({m1},{m2}) {ba} x {bb}: n={len(recs)} "
              f"min|D|={fit['min_abs_D_observed']:.4e} "
              f"D holdout rel rms={fit['D']['holdout_rel_rms']:.3e}", flush=True)

    out = {
        "status": "complete" if models else "no_pairs_with_enough_points",
        "provenance": {"commit": commit_hash(), "degree": args.degree,
                       "holdout_frac": args.holdout_frac, "seed": args.seed},
        "model": "H_eff 2x2; observables are T = w+ + w- and D = (w+ - w-)^2; "
                 "EP2 <=> D = 0; D is invariant under rescaling of the spectral "
                 "condition, unlike |dF/domega|",
        "symmetry_constraint": "diagonal sectors are even in delta, so D is fitted "
                               "as an analytic function of t = delta^2",
        "n_models": len(models),
        "min_abs_D_overall": min((m["min_abs_D_observed"] for m in models),
                                 default=None),
        "models": models,
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out}: {len(models)} models")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

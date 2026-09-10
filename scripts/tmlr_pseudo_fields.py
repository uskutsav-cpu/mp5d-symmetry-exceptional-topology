"""High-resolution *display* fields for the Figure 2 pseudospectral panels.

This script recomputes nothing that the paper reports.  The spectral points
``omega_used``, their condition numbers and their ``sigma_min`` values are read
from the ``results/pseudospectral_results.json`` checkpoint and are re-evaluated
only to assert that this run reproduces them; if they do not match, the script
fails rather than quietly writing a different number into the figure.

What it does add is the field that gets drawn.  The checkpoint stores a 21x21
map on a window centred on the *hint* frequency, which for the tightest
interaction put the tracked mode on the boundary of its own panel.  Here each
window is centred on ``omega_used`` itself and sampled finely, so the panel
shows the mode it is about, and both panels use the same window half-width so
that a visual comparison of their sharpness is a fair one.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.geometry import MPGeometry, sdelta_to_ab  # noqa: E402
from mp5d.radial.solver_c import SolverCProblem  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
CHECKPOINT = ROOT / "results" / "pseudospectral_results.json"
OUTFILE = ROOT / "results" / "pseudospectral_display_fields.json"

# Same solver settings as scripts/pseudospectrum.py; changing either of these
# would make the field a different object from the checkpointed one.
THETA_DEG = 60.0
L_CONTOUR = 60.0
RESOLUTION = 180

RTOL = 1e-6         # tolerance on reproducing the checkpointed diagnostics


def build(point: dict) -> SolverCProblem:
    a, b = sdelta_to_ab(point["s"], point["delta"])
    if abs(a - point["a"]) > 1e-12 or abs(b - point["b"]) > 1e-12:
        raise RuntimeError(f"{point['label']}: (a,b) disagrees with checkpoint")
    return SolverCProblem(MPGeometry(a=a, b=b, M=1.0), point["mu"],
                          point["m1"], point["m2"], point["ell"],
                          theta=math.radians(THETA_DEG), L=L_CONTOUR,
                          resolution=RESOLUTION)


def verify(prob: SolverCProblem, point: dict) -> None:
    """Refuse to draw a field whose solver disagrees with the checkpoint."""
    w = complex(*point["omega_used"])
    got, _ = prob.sigma_min(w)
    want = point["conditioning"]["sigma_min_rel"]
    if abs(got - want) > RTOL * abs(want):
        raise RuntimeError(f"{point['label']}: sigma_min {got:.6e} does not "
                           f"reproduce checkpoint {want:.6e}")
    print(f"  verified sigma_min_rel = {got:.6e} at omega_used")


def field(prob: SolverCProblem, centre: complex, half: float, n: int) -> dict:
    res = np.linspace(centre.real - half, centre.real + half, n)
    ims = np.linspace(centre.imag - half, centre.imag + half, n)
    g = np.empty((n, n))
    t0 = time.time()
    for i, y in enumerate(ims):
        for j, x in enumerate(res):
            try:
                s, _ = prob.sigma_min(complex(x, y))
            except Exception:  # noqa: BLE001
                s = float("nan")
            g[i, j] = s
        if (i + 1) % 10 == 0:
            print(f"    row {i + 1}/{n}  ({time.time() - t0:.0f}s)", flush=True)
    return {"re": res.tolist(), "im": ims.tolist(),
            "log10_sigma_min_rel": np.log10(np.maximum(g, 1e-300)).round(5).tolist()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=61)
    ap.add_argument("--half", type=float, default=0.30)
    ap.add_argument("--out", default=str(OUTFILE))
    args = ap.parse_args()

    chk = json.loads(CHECKPOINT.read_text())
    fields = {}
    t0 = time.time()
    for point in chk["points"]:
        label = point["label"]
        print(f"{label}: centring on omega_used = "
              f"{point['omega_used'][0]:.7f}{point['omega_used'][1]:+.7f}i")
        prob = build(point)
        verify(prob, point)
        fields[label] = {
            "centre": point["omega_used"],
            "half": args.half,
            "grid": args.grid,
            **field(prob, complex(*point["omega_used"]), args.half, args.grid),
        }

    out = {
        "status": "complete",
        "purpose": "display fields for the Figure 2 panels; reported "
                   "diagnostics come from pseudospectral_results.json",
        "provenance": {"solver": "solver-C", "precision": "double",
                       "resolution": RESOLUTION, "theta_deg": THETA_DEG,
                       "L": L_CONTOUR, "grid": args.grid, "half": args.half,
                       "runtime_s": round(time.time() - t0, 1)},
        "caveat": chk["caveat"],
        "fields": fields,
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out} in {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

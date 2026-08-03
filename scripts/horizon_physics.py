"""Horizon thermodynamics and the superradiance condition over the atlas.

Supplies the physical quantities the manuscript quotes: Hawking temperature,
horizon angular velocities, horizon area and entropy across the searched
domain, and the superradiance factor

    S = Re(omega) - (m1 Omega_a + m2 Omega_b)

whose sign decides whether a mode extracts rotational energy.  ``S < 0`` is the
superradiant regime, in which the horizon exponent
``sigma_+ = (omega - m1 Omega_a - m2 Omega_b)/(2 kappa)`` changes sign and the
group velocity at the horizon reverses relative to the co-rotating frame.

For a five-dimensional Myers-Perry black hole the horizon is a squashed
three-sphere of area

    A = 2 pi^2 (z_+ + a^2)(z_+ + b^2) / r_+ ,

which follows from the induced metric on the horizon; we verify it against the
Schwarzschild-Tangherlini limit ``A = 2 pi^2 r_+^3``.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.geometry import MPGeometry, sdelta_to_ab  # noqa: E402


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def horizon_area(g: MPGeometry) -> float:
    """Area of the (squashed) horizon three-sphere."""
    return 2.0 * math.pi**2 * (g.z_plus + g.a**2) * (g.z_plus + g.b**2) / g.r_plus


def entry(s: float, delta: float) -> dict:
    a, b = sdelta_to_ab(s, delta)
    g = MPGeometry(a=a, b=b, M=1.0)
    A = horizon_area(g)
    return {"s": s, "delta": delta, "a": a, "b": b,
            "r_plus": g.r_plus, "r_minus": g.r_minus,
            "extremality": g.extremality,
            "Omega_a": g.Omega_a, "Omega_b": g.Omega_b,
            "kappa": g.kappa, "T_H": g.T_H,
            "area": A, "entropy": A / 4.0}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/horizon_physics.json")
    ap.add_argument("--long-lived",
                    default="results/long_lived_branches_final.json")
    args = ap.parse_args()

    # --- sanity: Schwarzschild-Tangherlini limit ---------------------------
    g0 = MPGeometry(a=0.0, b=0.0, M=1.0)
    A0, A0_exact = horizon_area(g0), 2 * math.pi**2 * g0.r_plus**3
    assert abs(A0 - A0_exact) < 1e-12, (A0, A0_exact)
    assert abs(g0.T_H - 1.0 / (2 * math.pi * g0.r_plus)) < 1e-12

    corners = [(0.0, 0.0), (0.20, 0.0), (0.42, 0.0),
               (0.42, 0.20), (0.42, -0.20), (0.34, 0.20)]
    grid = [entry(s, d) for (s, d) in corners]

    # --- superradiance across the tracked long-lived branches -------------
    superradiant = []
    p = pathlib.Path(args.long_lived)
    if p.exists():
        data = json.loads(p.read_text())
        for br in data["branches"]:
            lab, bg = br["labels"], br["background"]
            a, b = sdelta_to_ab(bg["s"], bg["delta"])
            g = MPGeometry(a=a, b=b, M=1.0)
            crit = lab["m1"] * g.Omega_a + lab["m2"] * g.Omega_b
            pts = br["points"]
            if not pts:
                continue
            worst = min(pts, key=lambda q: q["omega"][0] - crit)
            superradiant.append({
                "labels": lab, "background": bg,
                "Omega_a": g.Omega_a, "Omega_b": g.Omega_b,
                "critical_frequency": crit,
                "min_superradiance_factor": worst["omega"][0] - crit,
                "at_mu": worst["mu"],
                "omega_re": worst["omega"][0],
                "is_superradiant": bool(worst["omega"][0] < crit),
                "min_damping": min(q["damping"] for q in pts),
            })

    n_sr = sum(1 for x in superradiant if x["is_superradiant"])
    out = {
        "status": "complete",
        "provenance": {"commit": commit_hash()},
        "conventions": {
            "area": "A = 2 pi^2 (z_+ + a^2)(z_+ + b^2)/r_+ ; "
                    "checked against A = 2 pi^2 r_+^3 at a=b=0",
            "entropy": "S = A/4 in units G = c = hbar = 1",
            "superradiance": "S_factor = Re(omega) - (m1 Omega_a + m2 Omega_b); "
                             "negative means superradiant",
        },
        "domain_corners": grid,
        "T_H_range": [min(e["T_H"] for e in grid), max(e["T_H"] for e in grid)],
        "Omega_a_range": [min(e["Omega_a"] for e in grid),
                          max(e["Omega_a"] for e in grid)],
        "n_superradiant_branches": n_sr,
        "superradiance": superradiant,
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")

    print("horizon data at the domain corners:")
    for e in grid:
        print(f"  s={e['s']:.2f} d={e['delta']:+.2f}  r+={e['r_plus']:.4f} "
              f"T_H={e['T_H']:.5f}  Om_a={e['Omega_a']:.4f} "
              f"Om_b={e['Omega_b']:.4f}  S={e['entropy']:.4f}")
    print(f"\nsuperradiant branches: {n_sr}/{len(superradiant)}")
    for x in superradiant[:6]:
        lab = x["labels"]
        print(f"  ({lab['m1']},{lab['m2']}) l={lab['ell']}: "
              f"crit={x['critical_frequency']:.4f} "
              f"min(Re w - crit)={x['min_superradiance_factor']:+.4f} "
              f"-> {'SUPERRADIANT' if x['is_superradiant'] else 'damped regime'}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

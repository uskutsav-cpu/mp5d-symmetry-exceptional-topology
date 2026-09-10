"""Enumerate static (Schwarzschild-Tangherlini) QNM roots, indexed by ``l`` only.

At ``a = b = 0`` the horizon angular velocities vanish and the angular
eigenvalue collapses to ``Lambda = l(l+2)``, so the radial problem depends on
the mode labels **only through ``l``**.  One root table therefore seeds every
``(m1, m2)`` sector, which is what makes the atlas affordable.

Overtones are labelled by increasing ``|Im omega|`` at this static point and
then carried by continuation.  They are *not* taken from the continued-fraction
inversion index: the inversion index changes the conditioning of the CF but not
its zero set, so it does not by itself select an overtone.  Relying on it
produces duplicate "overtones", which is exactly what a first attempt here did.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.radial.qnm import solve_qnm_cf  # noqa: E402


def labels_for(ell: int) -> tuple[int, int]:
    """Any ``(m1, m2)`` consistent with ``l = 2n + |m1| + |m2|``, ``n >= 0``.

    The static root depends only on ``l``, but the solver still needs labels
    obeying the parity constraint, so even ``l`` uses ``(0, 0)`` and odd ``l``
    uses ``(1, 0)``.  Using ``(0, 0)`` for odd ``l`` silently yields no roots.
    """
    return (0, 0) if ell % 2 == 0 else (1, 0)


def find_static_roots(ell: int, n_overtones: int, depth: int = 200,
                      dedupe: float = 1e-5) -> list[complex]:
    m1, m2 = labels_for(ell)
    found: list[complex] = []
    base = 0.34 * (ell + 2)
    reals = np.linspace(0.4 * base, 1.35 * base, 6)
    imags = np.linspace(-0.3, -2.8, 8)
    for gi in imags:
        for gr in reals:
            try:
                sol = solve_qnm_cf(0.0, 0.0, 0.0, m1, m2, ell, 0,
                                   initial_frequency=complex(gr, gi),
                                   depth=depth, depth_schedule=(depth,),
                                   angular_N=24, maxiter=40)
            except Exception:  # noqa: BLE001 - a failed seed is not an error
                continue
            if not sol.converged or sol.cf_residual > 1e-8:
                continue
            w = sol.omega
            if w.imag >= -1e-6 or w.real <= 1e-3:
                continue
            if all(abs(w - q) > dedupe for q in found):
                found.append(w)
    found.sort(key=lambda z: -z.imag)  # N = 0 is the least damped
    return found[:n_overtones]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--l-max", type=int, default=8)
    ap.add_argument("--overtones", type=int, default=4)
    ap.add_argument("--out", default="results/static_seed_table.json")
    args = ap.parse_args()

    t0 = time.time()
    table: dict[str, list[list[float]]] = {}
    for ell in range(args.l_max + 1):
        roots = find_static_roots(ell, args.overtones)
        table[str(ell)] = [[z.real, z.imag] for z in roots]
        print(f"l={ell}: " + ", ".join(f"{z.real:.6f}{z.imag:+.6f}i" for z in roots),
              flush=True)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "status": "complete",
        "description": "static ST5D roots; depends on l only because Omega_a=Omega_b=0 "
                       "and Lambda=l(l+2) at a=b=0",
        "overtone_convention": "N ordered by increasing |Im omega| at a=b=mu=0",
        "runtime_s": round(time.time() - t0, 1),
        "table": table,
    }, indent=1) + "\n")
    print(f"wrote {out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

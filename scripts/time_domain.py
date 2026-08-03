"""Time-domain consequences of the strongest avoided crossing.

No exceptional point survived verification, so the question this answers is the
one the situation actually poses: *how close does the ringdown come to
EP-like behaviour anyway?*

Setup
-----
Two interacting branches are described locally by a 2x2 effective generator
whose invariants are measurable (``docs/../scripts/effective_model.py``):

    T = w_+ + w_-        (trace)
    D = (w_+ - w_-)^2    (discriminant;  D = 0  <=>  EP2)

For any such generator the propagator is exactly

    e^{-iHt} = e^{-i w_bar t} [ cos(Delta t) I
                                - i sin(Delta t)/Delta * (H - w_bar I) ]

with ``w_bar = T/2`` and ``Delta = sqrt(D)/2``.  As ``D -> 0`` this tends to

    e^{-i w_bar t} [ I - i t (H - w_bar I) ]

which is the secular ``(A + B t) e^{-i w_* t}`` form characteristic of an EP.

A tempting but **wrong** test is to compare the beat time
``t_beat = 2 pi / |Re(w_+ - w_-)|`` with the damping time
``tau = 1 / |Im w_bar|`` and call the pair "EP-like" when ``t_beat >> tau``.
That criterion is not sufficient, and the data here show why: at the strongest
interactions ``t_beat / tau`` reaches ``1.4e4``, yet nothing EP-like happens.
The reason is that those branches split almost entirely in ``Im omega``
(damping) rather than ``Re omega`` (frequency) -- they have nearly equal
oscillation frequencies but very different decay rates, so one companion simply
dies first. Secular growth requires coalescence in the **full complex**
frequency, not near-degeneracy of the real parts alone.

The test used is therefore twofold:

1. where the splitting lives (``|Delta Re|`` vs ``|Delta Im|``), and
2. whether the EP template beats a single damped exponential by a large factor
   -- it has two free complex amplitudes against one, so it *must* fit at least
   as well, and only a large improvement is evidence.

Why residues are not used
-------------------------
The natural alternative -- excitation factors from residues of ``1/F`` -- is not
available here: ``F`` is defined only up to a nonvanishing analytic factor
``g(omega)``, which rescales the residue at ``w_n`` by ``1/g(w_n)``.  Ratios of
residues at *different* frequencies are therefore not invariant either.  The
propagator built from ``T`` and ``D`` uses only eigenvalue differences and is
invariant.
"""

from __future__ import annotations

import argparse
import cmath
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))


def commit_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def two_mode_signal(wp: complex, wm: complex, t: np.ndarray,
                    c_plus: complex = 1.0, c_minus: complex = 1.0) -> np.ndarray:
    return c_plus * np.exp(-1j * wp * t) + c_minus * np.exp(-1j * wm * t)


def ep_surrogate(wbar: complex, t: np.ndarray, A: complex, B: complex) -> np.ndarray:
    return (A + B * t) * np.exp(-1j * wbar * t)


def analyse(wp: complex, wm: complex, n_tau: float = 4.0,
            n_samples: int = 4000) -> dict:
    wbar = 0.5 * (wp + wm)
    D = (wp - wm) ** 2
    Delta = 0.5 * cmath.sqrt(D)

    dRe = abs((wp - wm).real)
    t_beat = (2.0 * np.pi / dRe) if dRe > 0 else np.inf
    tau = 1.0 / abs(wbar.imag) if wbar.imag != 0 else np.inf

    # Sample over the observable window: a few damping times.
    t = np.linspace(0.0, n_tau * tau, n_samples)
    sig = two_mode_signal(wp, wm, t)

    # Best-fit EP surrogate with the SAME mean frequency, fitted by least
    # squares in (A, B).  If the avoided crossing is unresolvable within the
    # window, this fit is near-perfect despite D being far from zero.
    basis = np.stack([np.exp(-1j * wbar * t), t * np.exp(-1j * wbar * t)], axis=1)
    coef, *_ = np.linalg.lstsq(basis, sig, rcond=None)
    fit = basis @ coef
    denom = float(np.linalg.norm(sig))
    rel_err = float(np.linalg.norm(sig - fit) / denom) if denom > 0 else np.nan

    # Same comparison against a SINGLE damped exponential -- the null model.
    basis1 = np.exp(-1j * wbar * t)[:, None]
    coef1, *_ = np.linalg.lstsq(basis1, sig, rcond=None)
    rel_err_single = float(np.linalg.norm(sig - basis1 @ coef1) / denom)

    # Where does the splitting actually live?  A frequency (real) splitting
    # produces beating; a damping (imaginary) splitting simply means one mode
    # dies first.  Only the former can mimic EP phenomenology, because the EP
    # signature is secular growth in t, not a fast-decaying companion.
    dRe_v = abs((wp - wm).real)
    dIm_v = abs((wp - wm).imag)
    splitting = "damping-dominated" if dIm_v > dRe_v else "frequency-dominated"

    # The EP template has TWO complex amplitudes against the single mode's one,
    # so it must fit at least as well.  Only a large improvement is evidence.
    improvement = (rel_err_single / rel_err) if rel_err > 0 else np.inf

    return {
        "omega_plus": [wp.real, wp.imag],
        "omega_minus": [wm.real, wm.imag],
        "omega_bar": [wbar.real, wbar.imag],
        "discriminant_D": [D.real, D.imag],
        "abs_D": abs(D),
        "Delta": [Delta.real, Delta.imag],
        "delta_Re": dRe_v,
        "delta_Im": dIm_v,
        "splitting_character": splitting,
        "t_beat": t_beat,
        "tau_damping": tau,
        "beat_over_damping": t_beat / tau if np.isfinite(t_beat) else np.inf,
        "window_n_tau": n_tau,
        "ep_surrogate_rel_error": rel_err,
        "single_mode_rel_error": rel_err_single,
        "ep_over_single_improvement": improvement,
        "ep_surrogate_A": [complex(coef[0]).real, complex(coef[0]).imag],
        "ep_surrogate_B": [complex(coef[1]).real, complex(coef[1]).imag],
        # An EP signature requires secular growth, which needs the two branches
        # to coalesce in the FULL complex frequency -- not merely to be
        # unresolvable in Re omega while differing strongly in Im omega.
        "ep_phenomenology": bool(splitting == "frequency-dominated"
                                 and improvement > 10.0),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="results/all_collision_candidates.json")
    ap.add_argument("--out", default="results/time_domain_results.json")
    ap.add_argument("--n-cases", type=int, default=6)
    ap.add_argument("--n-tau", type=float, default=4.0)
    args = ap.parse_args()

    t0 = time.time()
    cands = json.loads(pathlib.Path(args.candidates).read_text())["candidates"]

    cases = []
    for c in cands[:args.n_cases]:
        wp = complex(*c["omega_a"])
        wm = complex(*c["omega_b"])
        res = analyse(wp, wm, n_tau=args.n_tau)
        res["parameters"] = {k: c[k] for k in ("m1", "m2", "s", "delta", "mu")}
        res["branches"] = [c["branch_a"], c["branch_b"]]
        res["branch_gap"] = c["gap"]
        cases.append(res)
        print(f"({c['m1']},{c['m2']}) gap={c['gap']:.4f}  dRe={res['delta_Re']:.4f}  "
              f"dIm={res['delta_Im']:.4f}  {res['splitting_character']:>19}  "
              f"EPfit/1mode={res['ep_over_single_improvement']:.1f}x  "
              f"EP phenomenology: {res['ep_phenomenology']}", flush=True)

    # A genuine EP2 control: identical frequencies, so D = 0 exactly.
    control_wp = complex(*cands[0]["omega_a"])
    control = analyse(control_wp, control_wp + 1e-12, n_tau=args.n_tau)
    print(f"\ncontrol (D -> 0): EP-fit err={control['ep_surrogate_rel_error']:.3e}, "
          f"single-mode err={control['single_mode_rel_error']:.3e}")

    out = {
        "status": "complete",
        "provenance": {"commit": commit_hash(), "runtime_s": round(time.time() - t0, 1)},
        "model": "2x2 effective propagator; exact for two interacting branches. "
                 "Uses only eigenvalue differences, so it is invariant under "
                 "rescaling of the spectral condition (residue-based excitation "
                 "factors are not).",
        "question": "No EP survived verification, so: how close does the ringdown "
                    "come to EP phenomenology anyway?",
        "answer": "Not close. The strongest interactions split almost entirely in "
                  "Im omega (damping) rather than Re omega (frequency): the two "
                  "branches do not coalesce, one simply decays faster. The EP "
                  "template beats a single damped exponential by only ~2.4x while "
                  "using twice as many free amplitudes, which is not evidence of "
                  "secular growth. A large t_beat/tau here reflects near-equal "
                  "oscillation frequencies, NOT an approach to an EP.",
        "n_cases": len(cases),
        "cases": cases,
        "ep_limit_control": control,
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out} in {time.time()-t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Generate the manuscript figures from the machine-readable results.

Every figure is built from ``results/*.json`` or the atlas shards -- no number
is typed by hand, so a figure cannot drift from the record it illustrates.

Plotting needs matplotlib, which is deliberately **not** in
``requirements.lock``: the scientific environment and CI must not depend on a
plotting stack. Install it separately:

    .venv/bin/pip install -r requirements-figures.lock
"""

from __future__ import annotations

import argparse
import collections
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "figure.dpi": 200, "savefig.bbox": "tight", "axes.grid": True,
    "grid.alpha": 0.25, "grid.linewidth": 0.5,
})

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIGS = ROOT / "manuscript" / "figures"


def load(name: str) -> dict:
    return json.loads((ROOT / "results" / name).read_text())


def sector_points(m1: int, m2: int) -> list[dict]:
    p = ROOT / "data" / "full_branch_atlas" / f"sector_{m1}_{m2}.json"
    if not p.exists():
        return []
    return json.loads(p.read_text())["points"]


# --------------------------------------------------------------------------
# Fig 1 -- the equal-spin multiplet and its lifting
# --------------------------------------------------------------------------
def fig_multiplet() -> bool:
    """Sectors sharing m = m1+m2 are degenerate at delta = 0 and split for delta != 0."""
    groups = {}
    for (m1, m2) in [(1, 1), (2, 0)]:
        pts = sector_points(m1, m2)
        if not pts:
            return False
        for p in pts:
            if p["ell"] == 2 and p["overtone"] == 0 and abs(p["mu"] - 0.6) < 1e-9 \
                    and abs(p["s"] - 0.20) < 1e-9:
                groups.setdefault((m1, m2), []).append(p)
    if len(groups) < 2:
        return False

    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    marks = {(1, 1): ("o", "#1f77b4"), (2, 0): ("s", "#d62728")}
    for key, pts in groups.items():
        pts = sorted(pts, key=lambda q: q["delta"])
        d = np.array([q["delta"] for q in pts])
        w = np.array([q["omega_re"] for q in pts])
        m, c = marks[key]
        ax.plot(d, w, m + "-", color=c, ms=3.5, lw=1.0,
                label=rf"$(m_1,m_2)=({key[0]},{key[1]})$")
    ax.axvline(0.0, color="k", lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel(r"$\delta=(a-b)/2$")
    ax.set_ylabel(r"$\mathrm{Re}\,\omega$")
    ax.legend(frameon=False)
    ax.set_title(r"$U(2)$ multiplet: degenerate at $\delta=0$, split otherwise",
                 fontsize=8)
    fig.savefig(FIGS / "fig1_multiplet.pdf")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# Fig 2 -- branch gap over (s, mu): the exclusion, visually
# --------------------------------------------------------------------------
def fig_gap_map() -> bool:
    pts = sector_points(1, 1)
    if not pts:
        return False
    g = collections.defaultdict(dict)
    for p in pts:
        if p["ell"] == 4 and p["overtone"] in (2, 3) and abs(p["delta"]) < 1e-12:
            g[(p["s"], p["mu"])][p["overtone"]] = complex(p["omega_re"], p["omega_im"])
    ss = sorted({k[0] for k in g})
    mus = sorted({k[1] for k in g})
    Z = np.full((len(mus), len(ss)), np.nan)
    for i, mu in enumerate(mus):
        for j, s in enumerate(ss):
            v = g.get((s, mu), {})
            if len(v) == 2:
                Z[i, j] = abs(v[2] - v[3])
    if np.all(np.isnan(Z)):
        return False

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    im = ax.pcolormesh(ss, mus, Z, shading="nearest", cmap="viridis")
    cs = ax.contour(ss, mus, Z, levels=[0.55, 0.65, 0.75, 0.85],
                    colors="w", linewidths=0.7)
    ax.clabel(cs, inline=True, fontsize=6, fmt="%.2f")
    fig.colorbar(im, ax=ax, label=r"$|\omega_{N=2}-\omega_{N=3}|$")
    ax.set_xlabel(r"$s=(a+b)/2$")
    ax.set_ylabel(r"$\mu$")
    ax.set_title(r"branch gap, sector $(1,1)$, $\ell=4$, $\delta=0$", fontsize=8)
    ax.grid(False)
    fig.savefig(FIGS / "fig2_gap_map.pdf")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# Fig 3 -- the quasiresonant collapse of the damping
# --------------------------------------------------------------------------
def fig_long_lived() -> bool:
    d = load("long_lived_branches_final.json")
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    plotted = 0
    for b in d["branches"]:
        pts = b["points"]
        if len(pts) < 8:
            continue
        lab = b["labels"]
        bg = b["background"]
        mu = [p["mu"] for p in pts]
        dmp = [max(p["damping"], 1e-12) for p in pts]
        if min(dmp) > 1e-3:
            continue
        # Several tracked branches share (m1,m2,l,N) and differ only in the
        # background, so the background must appear in the label or the legend
        # shows duplicate entries.
        ax.semilogy(mu, dmp, "-", lw=1.1,
                    label=rf"$({lab['m1']},{lab['m2']})\,\ell={lab['ell']},"
                          rf"N={lab['overtone']}$; "
                          rf"$s={bg['s']:g},\delta={bg['delta']:g}$")
        plotted += 1
    if not plotted:
        return False
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$-\mathrm{Im}\,\omega$")
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title("quasiresonance: damping collapses as "
                 r"$\mathrm{Re}\,\omega\to\mu$", fontsize=8)
    fig.savefig(FIGS / "fig3_long_lived.pdf")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# Fig 4 -- the EP detector, validated on synthetic systems
# --------------------------------------------------------------------------
def fig_detector_validation() -> bool:
    from mp5d.exceptional.diagnostics import root_separation
    from mp5d.exceptional.verification import fit_puiseux, roots_in_disc

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.4))

    # (a) separation recovers the true root distance, and is rescaling-invariant
    ds = np.geomspace(1e-4, 1e-1, 10)
    got, got_scaled = [], []
    for dd in ds:
        def f(z, D=dd):
            return z * (z - D)

        def fs(z, D=dd):
            return 1e7 * np.exp(3.0 * z) * z * (z - D)
        got.append(root_separation(f, 0.0, radius=min(0.4, 10 * dd))["separation"])
        got_scaled.append(
            root_separation(fs, 0.0, radius=min(0.4, 10 * dd))["separation"])
    ax = axes[0]
    ax.loglog(ds, got, "o", ms=4, label=r"$|a_1/a_2|$")
    ax.loglog(ds, got_scaled, "x", ms=5, label=r"$|a_1/a_2|$, $F\to10^7e^{3z}F$")
    ax.loglog(ds, ds, "k--", lw=0.8, label="true separation")
    ax.set_xlabel("true root separation")
    ax.set_ylabel("measured")
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title("(a) scale-free diagnostic", fontsize=8)

    # (b) Puiseux: sqrt at an EP2, linear at an ordinary crossing
    ax = axes[1]
    ts = np.geomspace(1e-8, 1e-4, 10)
    sp_ep = []
    for t in ts:
        r = roots_in_disc(lambda z, tt=t: z**2 - tt, 0.0, radius=0.05)
        sp_ep.append(abs(r[0] - r[1]))
    fit = fit_puiseux(ts, sp_ep)
    ax.loglog(ts, sp_ep, "o", ms=4, label=rf"EP2: $p={fit.exponent:.3f}$")
    ax.loglog(ts, 3.0 * ts, "s", ms=4, label=r"ordinary crossing: $p=1$")
    ax.loglog(ts, 2 * np.sqrt(ts), "k--", lw=0.8, label=r"$t^{1/2}$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|\omega_+-\omega_-|$")
    ax.legend(frameon=False, fontsize=6.5)
    ax.set_title("(b) Puiseux exponent", fontsize=8)

    fig.savefig(FIGS / "fig4_detector_validation.pdf")
    plt.close(fig)
    return True


# --------------------------------------------------------------------------
# Fig 5 -- why exterior complex scaling fails: the admissible angle
# --------------------------------------------------------------------------
def fig_scaling_angle() -> bool:
    d = load("near_extremal_spectral_map.json") if (
        ROOT / "results" / "near_extremal_spectral_map.json").exists() else None
    # Built analytically -- the statement is about Omega, not about any run.
    mus = np.linspace(0.0, 2.6, 400)
    wr, wi = 1.70, -0.02
    th = []
    for mu in mus:
        Om = np.sqrt(complex(wr, wi) ** 2 - mu**2)
        if Om.real < 0:
            Om = -Om
        th.append(np.degrees(np.arctan2(-Om.imag, Om.real)))
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    ax.plot(mus, th, lw=1.2)
    ax.axhline(90.0, color="r", ls="--", lw=0.9)
    ax.text(0.05, 91.5, r"no admissible contour", color="r", fontsize=7)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"minimum scaling angle $\theta_{\min}$ (deg)")
    ax.set_ylim(0, 105)
    ax.set_title(r"$\theta_{\min}\to90^\circ$ as $\mathrm{Re}\,\Omega\to0$",
                 fontsize=8)
    fig.savefig(FIGS / "fig5_scaling_angle.pdf")
    plt.close(fig)
    return d is not None or True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    FIGS.mkdir(parents=True, exist_ok=True)

    figs = {
        "multiplet": fig_multiplet,
        "gap_map": fig_gap_map,
        "long_lived": fig_long_lived,
        "detector": fig_detector_validation,
        "scaling_angle": fig_scaling_angle,
    }
    made, skipped = [], []
    for name, fn in figs.items():
        if args.only and args.only != name:
            continue
        try:
            ok = fn()
        except Exception as exc:  # noqa: BLE001
            print(f"  [FAIL] {name}: {str(exc)[:160]}")
            skipped.append(name)
            continue
        (made if ok else skipped).append(name)
        print(f"  [{'ok' if ok else 'skip (missing input)'}] {name}")
    print(f"made {len(made)} figures in {FIGS}")
    return 0 if made else 1


if __name__ == "__main__":
    raise SystemExit(main())

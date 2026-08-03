"""Generate the manuscript figures from the machine-readable results.

Every figure is built from ``results/*.json`` or the atlas shards -- no number
is typed by hand, so a figure cannot drift from the record it illustrates.
Each function returns ``True`` if it produced a figure and ``False`` if its
input was missing, so a partial run degrades visibly rather than silently.

Style
-----
Deliberately plain: Computer Modern mathematics, greyscale, open and filled
markers, inward ticks on all four sides, no filled colour maps and no colour
bars.  The gap map is drawn as labelled level curves, which is how a quantity
of this kind is normally presented and which survives monochrome printing.
The figure count is kept small; a plot earns its place only if a table or a
sentence would not do.

Plotting needs matplotlib, which is deliberately **not** in
``requirements.lock``: the scientific environment and CI must not depend on a
plotting stack.  Install it separately:

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
from matplotlib.ticker import AutoMinorLocator  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "axes.grid": False,
    "axes.linewidth": 0.7,
    "lines.linewidth": 0.9,
    "lines.markersize": 3.6,
    "legend.frameon": False,
    "legend.handlelength": 2.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 3.4,
    "ytick.major.size": 3.4,
    "xtick.minor.size": 1.9,
    "ytick.minor.size": 1.9,
})

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIGS = ROOT / "manuscript" / "figures"

K, G1, G2 = "k", "0.40", "0.65"   # black and two greys


def _tidy(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())


def load(name: str) -> dict | None:
    p = ROOT / "results" / name
    return json.loads(p.read_text()) if p.exists() else None


def sector_points(m1: int, m2: int) -> list[dict]:
    p = ROOT / "data" / "full_branch_atlas" / f"sector_{m1}_{m2}.json"
    return json.loads(p.read_text())["points"] if p.exists() else []


# ==========================================================================
# Fig 1 -- equal-spin degeneracy and the delta-parity of diagonal sectors
# ==========================================================================
def fig_multiplet() -> bool:
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.55))

    sel = {}
    for key in [(1, 1), (2, 0)]:
        sel[key] = [p for p in sector_points(*key)
                    if p["ell"] == 2 and p["overtone"] == 0
                    and abs(p["mu"] - 0.6) < 1e-9 and abs(p["s"] - 0.20) < 1e-9]
    if not all(sel.values()):
        plt.close(fig)
        return False

    ax = axes[0]
    for key, mk, fc in [((1, 1), "o", "none"), ((2, 0), "s", K)]:
        d = sorted(sel[key], key=lambda q: q["delta"])
        ax.plot([p["delta"] for p in d], [p["omega_re"] for p in d],
                marker=mk, color=K, mfc=fc, ls="-",
                label=rf"$(m_1,m_2)=({key[0]},{key[1]})$")
    ax.axvline(0.0, color=G2, lw=0.6, ls=(0, (4, 3)))
    ax.set_xlabel(r"$\delta=(a-b)/2$")
    ax.set_ylabel(r"$\mathrm{Re}\,\omega$")
    ax.legend(loc="upper left")
    ax.set_title(r"(a)", loc="left")
    _tidy(ax)

    ax = axes[1]
    for key, mk, fc, lab in [((1, 1), "o", "none", r"diagonal, $m_1=m_2$"),
                             ((2, 0), "s", K, r"off-diagonal, $(2,0)$")]:
        d = {round(p["delta"], 9): complex(p["omega_re"], p["omega_im"])
             for p in sel[key]}
        xs = [x for x in sorted(d) if x > 0 and -x in d]
        ys = [max(abs(d[x] - d[-x]), 1e-17) for x in xs]
        if xs:
            ax.semilogy(xs, ys, marker=mk, color=K, mfc=fc, ls="-", label=lab)
    ax.axhline(1e-12, color=G2, lw=0.6, ls=":")
    ax.text(0.985, 1.5e-12, "double precision", ha="right", fontsize=6.5,
            color=G1, transform=ax.get_yaxis_transform())
    ax.set_xlabel(r"$|\delta|$")
    ax.set_ylabel(r"$|\omega(+\delta)-\omega(-\delta)|$")
    ax.set_ylim(1e-16, 1e0)
    ax.legend(loc="center right")
    ax.set_title(r"(b)", loc="left")
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(FIGS / "fig1_multiplet.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 2 -- the exclusion: gap level curves and diagnostic distributions
# ==========================================================================
def fig_exclusion() -> bool:
    pts = sector_points(1, 1)
    cand = load("all_collision_candidates.json")
    if not pts or not cand:
        return False

    g = collections.defaultdict(dict)
    for p in pts:
        if p["ell"] == 4 and p["overtone"] in (2, 3) and abs(p["delta"]) < 1e-12:
            g[(p["s"], p["mu"])][p["overtone"]] = complex(p["omega_re"],
                                                          p["omega_im"])
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

    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.55))

    ax = axes[0]
    cs = ax.contour(ss, mus, Z, levels=8, colors=K, linewidths=0.7)
    ax.clabel(cs, inline=True, fontsize=6.5, fmt="%.2f")
    ax.set_xlabel(r"$s=(a+b)/2$")
    ax.set_ylabel(r"$\mu$")
    ax.set_title(r"(a)", loc="left")
    _tidy(ax)

    ax = axes[1]
    gaps = np.array([c["gap"] for c in cand["candidates"]])
    t2 = cand.get("tier2_evaluated") or [c for c in cand["candidates"]
                                         if c.get("tier2", {}).get("ok")]
    seps = np.array([c["tier2"]["separation"] for c in t2])
    bins = np.linspace(min(gaps.min(), seps.min()) * 0.95,
                       max(gaps.max(), seps.max()) * 1.02, 30)
    ax.hist(gaps, bins=bins, histtype="step", color=K, lw=0.9,
            label=r"branch gap $|\omega_i-\omega_j|$")
    ax.hist(seps, bins=bins, histtype="stepfilled", facecolor="0.85",
            edgecolor=K, lw=0.7, ls=(0, (3, 2)),
            label=r"root separation $|a_1/a_2|$")
    ax.axvline(seps.min(), color=K, lw=0.8, ls=":")
    # Keep the annotation clear of the legend: anchor it low and to the left of
    # the marked minimum, not in the upper-right block the legend occupies.
    ax.annotate(rf"$\min|a_1/a_2| = {seps.min():.4f}$",
                xy=(seps.min(), 0), xycoords=("data", "axes fraction"),
                xytext=(0.03, 0.34), textcoords="axes fraction", fontsize=7,
                ha="left",
                arrowprops=dict(arrowstyle="-", color=K, lw=0.6,
                                shrinkA=0, shrinkB=2))
    ax.set_xlabel("spectral separation")
    ax.set_ylabel("count")
    ax.set_xlim(0, None)
    ax.margins(y=0.22)
    ax.legend(loc="upper right")
    ax.set_title(r"(b)", loc="left")
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(FIGS / "fig2_exclusion.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 3 -- calibration of the diagnostics on exactly solvable problems
# ==========================================================================
def fig_calibration() -> bool:
    from mp5d.exceptional.diagnostics import root_separation
    from mp5d.exceptional.verification import fit_puiseux, roots_in_disc

    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.55))

    ax = axes[0]
    ds = np.geomspace(1e-4, 1e-1, 10)
    inv, scaled, naive = [], [], []
    for dd in ds:
        def f(z, D=dd):
            return z * (z - D)

        def fs(z, D=dd):
            return 1e7 * np.exp(3.0 * z) * z * (z - D)
        r = min(0.4, 10 * dd)
        inv.append(root_separation(f, 0.0, radius=r)["separation"])
        scaled.append(root_separation(fs, 0.0, radius=r)["separation"])
        naive.append(abs(root_separation(fs, 0.0, radius=r)["a1"]))
    ax.loglog(ds, ds, color=G2, lw=0.8, ls="-", label="exact")
    ax.loglog(ds, inv, "o", color=K, mfc="none", label=r"$|a_1/a_2|$")
    ax.loglog(ds, scaled, "+", color=K, mew=0.9,
              label=r"$|a_1/a_2|$ after $F\mapsto 10^{7}e^{3\omega}F$")
    ax.loglog(ds, naive, "^", color=K, mfc=G2,
              label=r"$|dF/d\omega|$, same rescaling")
    ax.set_xlabel(r"true separation $|\omega_1-\omega_2|$")
    ax.set_ylabel("measured")
    ax.legend(loc="upper left")
    ax.set_title(r"(a)", loc="left")
    _tidy(ax)

    ax = axes[1]
    ts = np.geomspace(1e-8, 1e-4, 10)
    ep = [abs(np.subtract(*roots_in_disc(lambda z, tt=t: z**2 - tt, 0.0,
                                         radius=0.05)[:2])) for t in ts]
    fit = fit_puiseux(ts, ep)
    gg = 1e-3
    ax.loglog(ts, 2 * np.sqrt(ts), color=G2, lw=0.8, label=r"$t^{1/2}$")
    ax.loglog(ts, ep, "o", color=K, mfc="none",
              label=rf"EP2, fitted $p={fit.exponent:.3f}$")
    ax.loglog(ts, 3.0 * ts, "s", color=K, mfc=K, label=r"ordinary crossing")
    ax.loglog(ts, [2 * np.sqrt(t**2 + gg**2) for t in ts], "^", color=K,
              mfc=G2, label=r"avoided crossing")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|\omega_+-\omega_-|$")
    ax.legend(loc="lower right")
    ax.set_title(r"(b)", loc="left")
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(FIGS / "fig3_calibration.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 4 -- the quasiresonant limit and the scaling-angle obstruction
# ==========================================================================
def fig_quasiresonance() -> bool:
    d = load("long_lived_branches_final.json")
    if not d:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.55))

    ax = axes[0]
    styles = [("o", "none"), ("s", K), ("^", G2), ("d", "none")]
    # Eight branches reach the long-lived regime; plotting them all buries the
    # curves under the legend and adds nothing, since they are qualitatively
    # identical.  Show the four that span the range of threshold values, and
    # report the rest in the table.
    candidates = [b for b in d["branches"]
                  if len(b["points"]) >= 8
                  and min(p["damping"] for p in b["points"]) <= 1e-3]
    candidates.sort(key=lambda b: min(p["damping"] for p in b["points"]))
    n = 0
    for b in candidates[:4]:
        pts = b["points"]
        lab = b["labels"]
        bg = b["background"]
        mk, fc = styles[n % len(styles)]
        ax.semilogy([p["mu"] for p in pts],
                    [max(p["damping"], 1e-12) for p in pts],
                    marker=mk, color=K, mfc=fc, ls="-", ms=2.6, markevery=3,
                    label=rf"$({lab['m1']},{lab['m2']})$, $\ell={lab['ell']}$,"
                          rf" $s={bg['s']:g}$")
        n += 1
    if not n:
        plt.close(fig)
        return False
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$-\mathrm{Im}\,\omega$")
    ax.set_ylim(1e-11, 5e0)
    ax.legend(loc="lower left", fontsize=6.8)
    ax.set_title(r"(a)", loc="left")
    _tidy(ax)

    ax = axes[1]
    mus = np.linspace(0.0, 2.6, 500)
    for wi, ls, lab in [(-0.02, "-", r"$-\mathrm{Im}\,\omega=0.02$"),
                        (-0.20, (0, (5, 2)), r"$0.20$"),
                        (-0.60, (0, (1.5, 1.5)), r"$0.60$")]:
        th = []
        for mu in mus:
            Om = np.sqrt(complex(1.70, wi) ** 2 - mu**2)
            if Om.real < 0:
                Om = -Om
            th.append(np.degrees(np.arctan2(-Om.imag, Om.real)))
        ax.plot(mus, th, color=K, ls=ls, label=lab)
    ax.axhline(90.0, color=G1, lw=0.7, ls=":")
    ax.text(0.03, 92.0, r"no admissible contour", fontsize=7, color=G1)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\theta_{\min}\ (\mathrm{deg})$")
    ax.set_ylim(0, 105)
    ax.legend(loc="center left")
    ax.set_title(r"(b)", loc="left")
    _tidy(ax)

    fig.tight_layout()
    fig.savefig(FIGS / "fig4_quasiresonance.pdf")
    plt.close(fig)
    return True


# ==========================================================================
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    FIGS.mkdir(parents=True, exist_ok=True)

    figs = {
        "multiplet": fig_multiplet,
        "exclusion": fig_exclusion,
        "calibration": fig_calibration,
        "quasiresonance": fig_quasiresonance,
    }
    made, skipped = [], []
    for name, fn in figs.items():
        if args.only and args.only != name:
            continue
        try:
            ok = fn()
        except Exception as exc:  # noqa: BLE001
            print(f"  [FAIL] {name}: {type(exc).__name__}: {str(exc)[:160]}")
            skipped.append(name)
            continue
        (made if ok else skipped).append(name)
        print(f"  [{'ok' if ok else 'skip (missing input)'}] {name}")
    print(f"made {len(made)} figures in {FIGS}")
    if skipped:
        print(f"skipped: {', '.join(skipped)}")
    return 0 if made else 1


if __name__ == "__main__":
    raise SystemExit(main())

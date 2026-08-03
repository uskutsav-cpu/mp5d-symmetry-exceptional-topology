"""Generate the manuscript figures from the machine-readable results.

Every figure is built from ``results/*.json`` or the atlas shards -- no number
is typed by hand, so a figure cannot drift from the record it illustrates.
Each function returns ``True`` if it produced a figure and ``False`` if its
input was missing, so a partial run degrades visibly rather than silently.

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
from matplotlib.lines import Line2D  # noqa: E402

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "legend.fontsize": 7,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.titlesize": 8.5,
    "figure.dpi": 220, "savefig.bbox": "tight", "axes.grid": True,
    "grid.alpha": 0.22, "grid.linewidth": 0.5, "axes.linewidth": 0.7,
    "lines.linewidth": 1.1, "legend.frameon": False,
})

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIGS = ROOT / "manuscript" / "figures"
SECTORS = [(0, 0), (1, 1), (-1, -1), (1, 0), (2, 0), (2, 1), (2, 2)]
CMAP = "viridis"


def load(name: str) -> dict | None:
    p = ROOT / "results" / name
    return json.loads(p.read_text()) if p.exists() else None


def sector_points(m1: int, m2: int) -> list[dict]:
    p = ROOT / "data" / "full_branch_atlas" / f"sector_{m1}_{m2}.json"
    return json.loads(p.read_text())["points"] if p.exists() else []


def all_points() -> list[dict]:
    out = []
    for sec in SECTORS:
        out += sector_points(*sec)
    return out


# ==========================================================================
# Fig 1 -- the spectrum itself: branch ladder in the complex plane
# ==========================================================================
def fig_spectrum() -> bool:
    """Where the modes actually are, and how the atlas covers them."""
    pts = sector_points(1, 1)
    if not pts:
        return False
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9))

    # (a) complex plane, coloured by overtone, marker by l
    ax = axes[0]
    colors = plt.get_cmap("plasma")(np.linspace(0.05, 0.85, 4))
    marks = {2: "o", 4: "s", 6: "^"}
    for ell in (2, 4, 6):
        for N in range(4):
            sel = [p for p in pts if p["ell"] == ell and p["overtone"] == N]
            if not sel:
                continue
            ax.scatter([p["omega_re"] for p in sel], [p["omega_im"] for p in sel],
                       s=1.2, c=[colors[N]], marker=marks[ell], alpha=0.45,
                       linewidths=0, rasterized=True)
    ax.set_xlabel(r"$\mathrm{Re}\,\omega$")
    ax.set_ylabel(r"$\mathrm{Im}\,\omega$")
    ax.set_title(r"(a) atlas coverage, sector $(1,1)$", loc="left")
    handles = [Line2D([], [], marker="o", ls="", ms=4, color=colors[n],
                      label=rf"$N={n}$") for n in range(4)]
    handles += [Line2D([], [], marker=m, ls="", ms=4, color="0.35",
                       label=rf"$\ell={e}$") for e, m in marks.items()]
    ax.legend(handles=handles, ncol=2, loc="lower left", fontsize=6)

    # (b) the same branches as trajectories in mu, at fixed background
    ax = axes[1]
    for ell in (2, 4, 6):
        for N in range(4):
            sel = sorted([p for p in pts if p["ell"] == ell and p["overtone"] == N
                          and abs(p["delta"]) < 1e-12 and abs(p["s"] - 0.20) < 1e-9],
                         key=lambda q: q["mu"])
            if len(sel) < 3:
                continue
            ax.plot([p["omega_re"] for p in sel], [p["omega_im"] for p in sel],
                    "-", color=colors[N], marker=marks[ell], ms=2.6, lw=0.9)
    ax.set_xlabel(r"$\mathrm{Re}\,\omega$")
    ax.set_ylabel(r"$\mathrm{Im}\,\omega$")
    ax.set_title(r"(b) branch motion under $\mu$, $s=0.20$, $\delta=0$",
                 loc="left")

    fig.tight_layout()
    fig.savefig(FIGS / "fig1_spectrum.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 2 -- the U(2) multiplet and the delta/s parity classification
# ==========================================================================
def fig_multiplet() -> bool:
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.4))

    # (a) two sectors sharing m = 2 meet exactly at delta = 0
    ax = axes[0]
    style = {(1, 1): ("o", "#1f77b4"), (2, 0): ("s", "#d62728")}
    ok = 0
    for key in [(1, 1), (2, 0)]:
        sel = [p for p in sector_points(*key)
               if p["ell"] == 2 and p["overtone"] == 0
               and abs(p["mu"] - 0.6) < 1e-9 and abs(p["s"] - 0.20) < 1e-9]
        if not sel:
            continue
        ok += 1
        sel.sort(key=lambda q: q["delta"])
        m, c = style[key]
        ax.plot([p["delta"] for p in sel], [p["omega_re"] for p in sel],
                m + "-", color=c, ms=3.4, label=rf"$({key[0]},{key[1]})$")
    if ok < 2:
        plt.close(fig)
        return False
    ax.axvline(0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xlabel(r"$\delta$")
    ax.set_ylabel(r"$\mathrm{Re}\,\omega$")
    ax.set_title(r"(a) $m_1{+}m_2=2$ multiplet", loc="left")
    ax.legend()

    # (b) delta-parity: diagonal even, off-diagonal not
    ax = axes[1]
    for key, lab in [((1, 1), "diagonal $(1,1)$"), ((2, 0), "off-diag $(2,0)$")]:
        sel = [p for p in sector_points(*key)
               if p["ell"] == 2 and p["overtone"] == 0
               and abs(p["mu"] - 0.6) < 1e-9 and abs(p["s"] - 0.20) < 1e-9]
        d = {round(p["delta"], 9): complex(p["omega_re"], p["omega_im"])
             for p in sel}
        xs, ys = [], []
        for dd in sorted(x for x in d if x > 0):
            if -dd in d:
                xs.append(dd)
                ys.append(max(abs(d[dd] - d[-dd]), 1e-17))
        if xs:
            ax.semilogy(xs, ys, "o-", ms=3.4, color=style[key][1], label=lab)
    ax.axhline(1e-12, color="0.4", lw=0.7, ls=":")
    ax.text(0.02, 1.6e-12, "machine precision", fontsize=6, color="0.35",
            transform=ax.get_yaxis_transform())
    ax.set_xlabel(r"$|\delta|$")
    ax.set_ylabel(r"$|\omega(+\delta)-\omega(-\delta)|$")
    ax.set_title(r"(b) $\delta$-parity", loc="left")
    ax.legend(loc="center right")

    # (c) multiplet dimension l+1, exact combinatorics
    ax = axes[2]
    ells = range(0, 13)
    dims = []
    for ell in ells:
        ms = [m for m in range(-ell, ell + 1) if (ell - abs(m)) % 2 == 0]
        d = []
        for m in ms:
            cnt = sum(1 for m1 in range(-ell - 1, ell + 2)
                      if abs(m1) + abs(m - m1) <= ell
                      and (ell - abs(m1) - abs(m - m1)) % 2 == 0)
            d.append(cnt)
        dims.append((ell, d))
    for ell, d in dims:
        ax.plot([ell] * len(d), d, "o", ms=2.6, color="#1f77b4", alpha=0.6)
    ax.plot(list(ells), [e + 1 for e in ells], "k--", lw=0.9,
            label=r"$\ell+1$")
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel("multiplet dimension")
    ax.set_title(r"(c) dimension is $\ell+1$ for every $m$", loc="left")
    ax.legend()

    fig.tight_layout()
    fig.savefig(FIGS / "fig2_multiplet.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 3 -- the exclusion: gap landscape across every sector
# ==========================================================================
def _gap_grid(pts, ell, na, nb):
    g = collections.defaultdict(dict)
    for p in pts:
        if p["ell"] == ell and p["overtone"] in (na, nb) and abs(p["delta"]) < 1e-12:
            g[(p["s"], p["mu"])][p["overtone"]] = complex(p["omega_re"],
                                                          p["omega_im"])
    ss = sorted({k[0] for k in g})
    mus = sorted({k[1] for k in g})
    Z = np.full((len(mus), len(ss)), np.nan)
    for i, mu in enumerate(mus):
        for j, s in enumerate(ss):
            v = g.get((s, mu), {})
            if len(v) == 2:
                Z[i, j] = abs(v[na] - v[nb])
    return ss, mus, Z


def fig_gap_landscape() -> bool:
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.6), sharex=True, sharey=True)
    vmin, vmax = np.inf, -np.inf
    panels = []
    for sec in SECTORS:
        pts = sector_points(*sec)
        if not pts:
            continue
        lmin = abs(sec[0]) + abs(sec[1])
        ell = lmin + 4 if lmin + 4 <= 8 else lmin
        ss, mus, Z = _gap_grid(pts, ell, 2, 3)
        if Z.size == 0 or np.all(np.isnan(Z)):
            continue
        panels.append((sec, ell, ss, mus, Z))
        vmin = min(vmin, np.nanmin(Z))
        vmax = max(vmax, np.nanmax(Z))
    if not panels:
        plt.close(fig)
        return False

    for ax, (sec, ell, ss, mus, Z) in zip(axes.ravel(), panels,
                                          strict=False):
        im = ax.pcolormesh(ss, mus, Z, shading="nearest", cmap=CMAP,
                           vmin=vmin, vmax=vmax)
        ax.set_title(rf"$({sec[0]},{sec[1]})$, $\ell={ell}$", fontsize=7.5)
        ax.grid(False)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel(r"$s$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\mu$")
    cb = fig.colorbar(im, ax=axes, fraction=0.022, pad=0.015)
    cb.set_label(r"$|\omega_{N=2}-\omega_{N=3}|$ at $\delta=0$")
    fig.savefig(FIGS / "fig3_gap_landscape.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 4 -- the bounds: distributions of the three invariant diagnostics
# ==========================================================================
def fig_bounds() -> bool:
    cand = load("all_collision_candidates.json")
    eff = load("effective_models.json")
    if not cand:
        return False
    gaps = np.array([c["gap"] for c in cand["candidates"]])
    # Use the full tier-2 population, not just the top-200 `candidates`: the
    # random subsample lies outside that slice, and drawing from it alone made
    # the panel disagree with the reported bound.
    t2 = cand.get("tier2_evaluated") or [c for c in cand["candidates"]
                                         if c.get("tier2", {}).get("ok")]
    seps = np.array([c["tier2"]["separation"] for c in t2])

    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.2))

    ax = axes[0]
    ax.hist(gaps, bins=28, color="#1f77b4", alpha=0.85)
    ax.axvline(gaps.min(), color="r", lw=1.1)
    ax.annotate(rf"$\min={gaps.min():.4f}$", xy=(gaps.min(), 0),
                xytext=(0.30, 0.82), textcoords="axes fraction", color="r",
                fontsize=6.5, arrowprops=dict(arrowstyle="->", color="r", lw=0.7))
    ax.set_xlabel(r"branch gap $|\omega_i-\omega_j|$")
    ax.set_ylabel("count")
    ax.set_title("(a) tier 1", loc="left")

    ax = axes[1]
    if seps.size:
        ax.hist(seps, bins=20, color="#2ca02c", alpha=0.85)
        ax.axvline(seps.min(), color="r", lw=1.1)
        ax.annotate(rf"$\min={seps.min():.4f}$", xy=(seps.min(), 0),
                    xytext=(0.30, 0.82), textcoords="axes fraction", color="r",
                    fontsize=6.5,
                    arrowprops=dict(arrowstyle="->", color="r", lw=0.7))
    ax.set_xlabel(r"root separation $|a_1/a_2|$")
    ax.set_title("(b) tier 2 (scale-free)", loc="left")

    ax = axes[2]
    if eff and eff.get("models"):
        vals = np.array([m["min_abs_D_observed"] for m in eff["models"]])
        labs = [rf"$({m['pair']['m1']},{m['pair']['m2']})$"
                for m in eff["models"]]
        ax.barh(range(len(vals)), vals, color="#9467bd", alpha=0.85)
        ax.set_yticks(range(len(vals)))
        ax.set_yticklabels(labs, fontsize=6.5)
        ax.axvline(vals.min(), color="r", lw=1.1)
        ax.set_xlabel(r"$\min|D|$ per branch pair")
        ax.set_title("(c) discriminant", loc="left")
    fig.tight_layout()
    fig.savefig(FIGS / "fig4_bounds.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 5 -- detector validation on synthetic systems
# ==========================================================================
def fig_detector_validation() -> bool:
    from mp5d.exceptional.diagnostics import root_separation
    from mp5d.exceptional.verification import fit_puiseux, roots_in_disc

    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.2))

    ax = axes[0]
    ds = np.geomspace(1e-4, 1e-1, 10)
    got, got_scaled, naive = [], [], []
    for dd in ds:
        def f(z, D=dd):
            return z * (z - D)

        def fs(z, D=dd):
            return 1e7 * np.exp(3.0 * z) * z * (z - D)
        r = min(0.4, 10 * dd)
        got.append(root_separation(f, 0.0, radius=r)["separation"])
        got_scaled.append(root_separation(fs, 0.0, radius=r)["separation"])
        naive.append(abs(root_separation(fs, 0.0, radius=r)["a1"]))
    ax.loglog(ds, got, "o", ms=4, label=r"$|a_1/a_2|$")
    ax.loglog(ds, got_scaled, "x", ms=5, mew=1.1,
              label=r"same, $F\!\to\!10^7e^{3z}F$")
    ax.loglog(ds, naive, "^", ms=4, color="0.55",
              label=r"$|dF/d\omega|$ (not invariant)")
    ax.loglog(ds, ds, "k--", lw=0.8, label="truth")
    ax.set_xlabel("true root separation")
    ax.set_ylabel("measured")
    ax.set_title("(a) invariance", loc="left")
    ax.legend(fontsize=6)

    ax = axes[1]
    ts = np.geomspace(1e-8, 1e-4, 10)
    sp = [abs(np.subtract(*roots_in_disc(lambda z, tt=t: z**2 - tt, 0.0,
                                         radius=0.05)[:2])) for t in ts]
    fit = fit_puiseux(ts, sp)
    ax.loglog(ts, sp, "o", ms=4, label=rf"EP2, fit $p={fit.exponent:.3f}$")
    ax.loglog(ts, 3.0 * ts, "s", ms=4, label=r"ordinary crossing, $p=1$")
    g = 1e-3
    ax.loglog(ts, [2 * np.sqrt(t**2 + g**2) for t in ts], "d", ms=3.6,
              label="avoided crossing")
    ax.loglog(ts, 2 * np.sqrt(ts), "k--", lw=0.8, label=r"$t^{1/2}$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|\omega_+-\omega_-|$")
    ax.set_title("(b) Puiseux exponent", loc="left")
    ax.legend(fontsize=6, loc="lower right")

    # (c) monodromy: EP2 swaps, avoided crossing does not
    ax = axes[2]
    for lab, fn, col in [
        ("EP2", lambda z, t: z**2 - t, "#d62728"),
        ("avoided", lambda z, t: z**2 - (t**2 + 0.01**2), "#1f77b4"),
    ]:
        def root_fn(p, F=fn):
            return roots_in_disc(lambda z: F(z, complex(p[0], p[1])), 0.0,
                                 radius=0.6, max_roots=3)
        # transport the two roots around a small loop and plot their phase
        from mp5d.exceptional.verification import track_roots_on_loop
        loop = [[1e-3 * np.cos(a), 1e-3 * np.sin(a)]
                for a in np.linspace(0, 2 * np.pi, 48)]
        hist = track_roots_on_loop(root_fn, loop, n_roots=2)
        for k in range(2):
            ax.plot(np.linspace(0, 1, len(hist)),
                    [np.angle(h[k]) for h in hist],
                    "-" if k == 0 else "--", color=col, lw=1.1,
                    label=f"{lab}" if k == 0 else None)
    ax.set_xlabel("fraction of loop")
    ax.set_ylabel(r"$\arg\omega$")
    ax.set_title("(c) monodromy", loc="left")
    ax.legend(fontsize=6)

    fig.tight_layout()
    fig.savefig(FIGS / "fig5_detector_validation.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 6 -- quasiresonance and the scaling-angle obstruction
# ==========================================================================
def fig_quasiresonance() -> bool:
    d = load("long_lived_branches_final.json")
    if not d:
        return False
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.3))

    ax = axes[0]
    n = 0
    for b in d["branches"]:
        pts = b["points"]
        if len(pts) < 8 or min(p["damping"] for p in pts) > 1e-3:
            continue
        lab = b["labels"]
        bg = b["background"]
        ax.semilogy([p["mu"] for p in pts],
                    [max(p["damping"], 1e-12) for p in pts], "-",
                    label=rf"$({lab['m1']},{lab['m2']})\ell{lab['ell']}$"
                          rf"$,s{bg['s']:g}$")
        n += 1
    if not n:
        plt.close(fig)
        return False
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$-\mathrm{Im}\,\omega$")
    ax.set_title("(a) damping collapse", loc="left")
    ax.legend(fontsize=5.6)

    ax = axes[1]
    for b in d["branches"]:
        pts = b["points"]
        if len(pts) < 8:
            continue
        ax.plot([p["mu"] for p in pts],
                [p["omega"][0] - p["mu"] for p in pts], "-", lw=0.9)
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\mathrm{Re}\,\omega-\mu$")
    ax.set_title(r"(b) approach to threshold", loc="left")

    ax = axes[2]
    mus = np.linspace(0.0, 2.6, 400)
    for wr, wi, lab in [(1.70, -0.02, r"$-\mathrm{Im}\,\omega=0.02$"),
                        (1.70, -0.20, r"$0.20$"),
                        (1.70, -0.60, r"$0.60$")]:
        th = []
        for mu in mus:
            Om = np.sqrt(complex(wr, wi) ** 2 - mu**2)
            if Om.real < 0:
                Om = -Om
            th.append(np.degrees(np.arctan2(-Om.imag, Om.real)))
        ax.plot(mus, th, lw=1.1, label=lab)
    ax.axhline(90.0, color="r", ls="--", lw=0.9)
    ax.text(0.03, 0.93, "no admissible contour", color="r", fontsize=6,
            transform=ax.transAxes)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\theta_{\min}$ (deg)")
    ax.set_ylim(0, 105)
    ax.set_title(r"(c) scaling-angle limit", loc="left")
    ax.legend(fontsize=6, loc="center left")

    fig.tight_layout()
    fig.savefig(FIGS / "fig6_quasiresonance.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 7 -- pseudospectrum, and the effective model's out-of-sample skill
# ==========================================================================
def fig_pseudospectrum() -> bool:
    p = load("pseudospectral_results.json")
    eff = load("effective_models.json")
    if not p:
        return False
    maps = p["sigma_maps"]
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.3))

    for ax, key, title in [(axes[0], "ordinary_control", "(a) ordinary mode"),
                           (axes[1], "tightest_interaction",
                            "(b) strongest interaction")]:
        if key not in maps:
            continue
        m = maps[key]
        Z = np.array(m["log10_sigma_min_rel"])
        im = ax.contourf(m["re"], m["im"], Z, levels=18, cmap="magma")
        ax.contour(m["re"], m["im"], Z, levels=8, colors="w", linewidths=0.4,
                   alpha=0.6)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02,
                     label=r"$\log_{10}\sigma_{\min}/\sigma_{\max}$")
        ax.set_xlabel(r"$\mathrm{Re}\,\omega$")
        ax.set_ylabel(r"$\mathrm{Im}\,\omega$")
        ax.set_title(title, loc="left")
        ax.grid(False)

    ax = axes[2]
    if eff and eff.get("models"):
        labs, tr, ho, diag = [], [], [], []
        for m in eff["models"]:
            labs.append(rf"$({m['pair']['m1']},{m['pair']['m2']})$")
            tr.append(m["D"]["train_rel_rms"])
            ho.append(m["D"]["holdout_rel_rms"])
            diag.append(m["pair"]["diagonal_sector"])
        x = np.arange(len(labs))
        ax.semilogy(x, tr, "o", ms=4, label="train")
        for xx, yy, dd in zip(x, ho, diag, strict=True):
            ax.semilogy([xx], [yy], "s", ms=5,
                        color="#d62728" if dd else "#7f7f7f")
        ax.set_xticks(x)
        ax.set_xticklabels(labs, fontsize=6.5)
        ax.set_ylabel(r"relative RMS in $D$")
        ax.set_title(r"(c) $D$ fitted in $t=\delta^2$", loc="left")
        handles = [Line2D([], [], marker="o", ls="", color="#1f77b4",
                          label="train"),
                   Line2D([], [], marker="s", ls="", color="#d62728",
                          label="hold-out, diagonal"),
                   Line2D([], [], marker="s", ls="", color="#7f7f7f",
                          label="hold-out, off-diagonal")]
        ax.legend(handles=handles, fontsize=6)

    fig.tight_layout()
    fig.savefig(FIGS / "fig7_pseudospectrum.pdf")
    plt.close(fig)
    return True


# ==========================================================================
# Fig 8 -- time domain: no EP phenomenology
# ==========================================================================
def fig_time_domain() -> bool:
    d = load("time_domain_results.json")
    if not d or not d.get("cases"):
        return False
    c = d["cases"][0]
    wp = complex(*c["omega_plus"])
    wm = complex(*c["omega_minus"])
    wbar = 0.5 * (wp + wm)
    tau = c["tau_damping"]
    t = np.linspace(0, 4 * tau, 3000)
    sig = np.exp(-1j * wp * t) + np.exp(-1j * wm * t)
    basis = np.stack([np.exp(-1j * wbar * t), t * np.exp(-1j * wbar * t)], 1)
    coef, *_ = np.linalg.lstsq(basis, sig, rcond=None)
    ep = basis @ coef
    b1 = np.exp(-1j * wbar * t)[:, None]
    c1, *_ = np.linalg.lstsq(b1, sig, rcond=None)
    one = (b1 @ c1).ravel()

    fig, axes = plt.subplots(1, 3, figsize=(6.9, 2.2))
    ax = axes[0]
    ax.plot(t, sig.real, lw=1.2, label="two-mode signal")
    ax.plot(t, ep.real, "--", lw=1.0, label=r"EP template $(A{+}Bt)e^{-i\bar\omega t}$")
    ax.plot(t, one.real, ":", lw=1.0, label="single mode")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\mathrm{Re}\,\psi$")
    ax.set_title("(a) waveform", loc="left")
    ax.legend(fontsize=6)

    ax = axes[1]
    ax.semilogy(t, np.abs(sig - ep) / np.abs(sig).max(), lw=1.0,
                label="EP template")
    ax.semilogy(t, np.abs(sig - one) / np.abs(sig).max(), lw=1.0,
                label="single mode")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel("relative residual")
    ax.set_title("(b) fit residual", loc="left")
    ax.legend(fontsize=6)

    ax = axes[2]
    dre = [x["delta_Re"] for x in d["cases"]]
    dim = [x["delta_Im"] for x in d["cases"]]
    ax.scatter(dre, dim, s=22, c="#d62728", zorder=3)
    lim = max(max(dre), max(dim)) * 1.15
    ax.plot([0, lim], [0, lim], "k--", lw=0.8, label=r"$|\Delta\mathrm{Re}|=|\Delta\mathrm{Im}|$")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(r"$|\mathrm{Re}(\omega_+-\omega_-)|$")
    ax.set_ylabel(r"$|\mathrm{Im}(\omega_+-\omega_-)|$")
    ax.set_title("(c) splitting is in the damping", loc="left")
    ax.legend(fontsize=6)

    fig.tight_layout()
    fig.savefig(FIGS / "fig8_time_domain.pdf")
    plt.close(fig)
    return True


# ==========================================================================
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    FIGS.mkdir(parents=True, exist_ok=True)

    figs = {
        "spectrum": fig_spectrum,
        "multiplet": fig_multiplet,
        "gap_landscape": fig_gap_landscape,
        "bounds": fig_bounds,
        "detector": fig_detector_validation,
        "quasiresonance": fig_quasiresonance,
        "pseudospectrum": fig_pseudospectrum,
        "time_domain": fig_time_domain,
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

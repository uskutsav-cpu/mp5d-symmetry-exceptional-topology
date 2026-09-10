"""Export the two TMLR figures' data as plain tables for PGFPlots.

The figures themselves are drawn by TikZ/PGFPlots inside the manuscript, so
they are typeset in the document's own font and are true vector art.  This
script only moves numbers out of the repository record and into ``.dat`` files;
it does no smoothing, gridding or interpolation.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "tmlr" / "data"

FLOOR = 1.58        # base plane of the Fig. 1 axis, in Re omega
NLEVELS = 9         # contour levels drawn over each pseudospectral field


def _contour_generator(X, Y, Z):
    """Level-set extractor: contourpy if present, else matplotlib's."""
    try:
        from contourpy import LineType, contour_generator
        # LineType.Separate makes lines() return a plain list of (N,2) arrays.
        cg = contour_generator(X, Y, Z, line_type=LineType.Separate)
        return cg.lines
    except ImportError:  # pragma: no cover - fallback path
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()

        def gen(lv):
            cs = ax.contour(X, Y, Z, levels=[lv])
            return [p.vertices for c in cs.collections for p in c.get_paths()]
        return gen


def sector_points(m1: int, m2: int) -> list[dict]:
    p = ROOT / "data" / "full_branch_atlas" / f"sector_{m1}_{m2}.json"
    if not p.exists():
        raise FileNotFoundError(f"missing atlas shard {p}; do not substitute data")
    return json.loads(p.read_text())["points"]


# ------------------------------------------------------- Fig 1: continuation --
def export_continuation() -> None:
    """The measured equal-spin crossing, on its true line in the (a,b) plane."""
    for m1, m2 in [(1, 1), (2, 0)]:
        rows = [p for p in sector_points(m1, m2)
                if p["ell"] == 2 and p["overtone"] == 0
                and abs(p["mu"] - 0.6) < 1e-9 and abs(p["s"] - 0.20) < 1e-9]
        if not rows:
            raise RuntimeError(f"no atlas rows for sector ({m1},{m2})")
        rows.sort(key=lambda q: q["delta"])
        lines = ["a b reom imom delta"]
        for p in rows:
            lines.append(f"{p['a']:.10f} {p['b']:.10f} {p['omega_re']:.10f} "
                         f"{p['omega_im']:.10f} {p['delta']:.10f}")
        (OUT / f"cont_{m1}_{m2}.dat").write_text("\n".join(lines) + "\n")

        # A two-scanline mesh sweeping each curve down to the floor: the
        # "ribbon on its true line" of the figure brief.  It adds no data --
        # the lower edge is the floor, the upper edge is the measurement.
        curt = ["a b z"]
        for p in rows:
            curt.append(f"{p['a']:.10f} {p['b']:.10f} {FLOOR:.10f}")
        curt.append("")
        for p in rows:
            curt.append(f"{p['a']:.10f} {p['b']:.10f} {p['omega_re']:.10f}")
        (OUT / f"curtain_{m1}_{m2}.dat").write_text("\n".join(curt) + "\n")

        print(f"cont_{m1}_{m2}.dat  {len(rows)} rows  "
              f"Re w in [{min(p['omega_re'] for p in rows):.4f}, "
              f"{max(p['omega_re'] for p in rows):.4f}]  "
              f"Im w in [{min(p['omega_im'] for p in rows):.4f}, "
              f"{max(p['omega_im'] for p in rows):.4f}]")


# ----------------------------------------------------- Fig 2: pseudospectra --
def _display_maps(d: dict) -> tuple[dict, str]:
    """The fields Figure 2 draws.

    ``scripts/tmlr_pseudo_fields.py`` recomputes each panel's field on a fine
    grid centred on the panel's own ``omega_used``; the checkpoint's coarse
    maps are centred on the *hint* instead, which leaves the tracked mode on
    the boundary of its own panel.  Prefer the display fields when they exist,
    and say which set was used so the provenance is never ambiguous.
    """
    p = ROOT / "results" / "pseudospectral_display_fields.json"
    if p.exists():
        f = json.loads(p.read_text())
        return f["fields"], f"display fields, {f['provenance']['grid']}^2"
    return d["sigma_maps"], "checkpoint maps, 21^2"


def export_pseudospectra() -> None:
    """log10 sigma_min/sigma_max on the stored complex-frequency grids."""
    d = json.loads((ROOT / "results" / "pseudospectral_results.json").read_text())
    pts = {p["label"]: p for p in d["points"]}

    maps, source = _display_maps(d)
    print(f"Fig. 2 fields: {source}")
    facts: list[str] = []
    for label, sm in maps.items():
        re = np.array(sm["re"])
        im = np.array(sm["im"])
        Z = np.array(sm["log10_sigma_min_rel"])
        if Z.shape != (im.size, re.size):
            raise RuntimeError(f"{label}: map shape {Z.shape} does not match axes")

        # No blank lines here: the mesh dimensions are declared in the figure
        # via mesh/rows and mesh/cols, which is the unambiguous form.
        lines = ["re im logsig"]
        for i, y in enumerate(im):
            for j, x in enumerate(re):
                lines.append(f"{x:.10f} {y:.10f} {Z[i, j]:.6f}")
        (OUT / f"pseudo_{label}.dat").write_text("\n".join(lines) + "\n")

        # Contour polylines, extracted by marching squares on the same grid.
        # pgfplots' own contour engines need gnuplot or LuaTeX; tectonic runs
        # XeTeX, so the level sets are precomputed here and drawn as paths.
        levels = np.linspace(Z.min(), Z.max(), NLEVELS + 2)[1:-1]
        X, Y = np.meshgrid(re, im)
        cg = _contour_generator(X, Y, Z)
        paths: list[str] = []
        for lv in levels:
            for seg in cg(lv):
                if len(seg) < 2:
                    continue
                # A "nan nan" row separates polylines; the figure reads this
                # file with "unbounded coords=jump" so the paths do not join
                # up.  (The blank-line form pgfplots also accepts makes it
                # typeset a stray character in a null font on every break.)
                if paths:
                    paths.append("nan nan")
                paths.extend(f"{x:.8f} {y:.8f}" for x, y in seg)
        (OUT / f"pseudocont_{label}.dat").write_text(
            "re im\n" + "\n".join(paths) + "\n")
        (OUT / f"pseudolevels_{label}.txt").write_text(
            " ".join(f"{lv:.4f}" for lv in levels) + "\n")

        p = pts[label]
        meta = {
            "label": label,
            "s": p["s"], "delta": p["delta"], "mu": p["mu"],
            "m1": p["m1"], "m2": p["m2"], "ell": p["ell"],
            "omega_used": p["omega_used"],
            "condition_number": p["conditioning"]["condition_number"],
            "sigma_min_rel": p["conditioning"]["sigma_min_rel"],
            "re_range": [float(re.min()), float(re.max())],
            "im_range": [float(im.min()), float(im.max())],
            "z_range": [float(Z.min()), float(Z.max())],
            "branch_gap": p.get("branch_gap"),
        }
        (OUT / f"pseudo_{label}.json").write_text(json.dumps(meta, indent=1))
        print(f"pseudo_{label}.dat  {re.size}x{im.size}  "
              f"log10 sig in [{Z.min():.3f}, {Z.max():.3f}]  "
              f"kappa = {meta['condition_number']:.4g}")

        # Axis limits, colour range and mesh size as macros, so the figure
        # source cannot fall out of step with the field it draws.
        tag = "A" if label == "ordinary_control" else "B"
        facts += [
            rf"\newcommand{{\ps{tag}xmin}}{{{re.min():.6f}}}",
            rf"\newcommand{{\ps{tag}xmax}}{{{re.max():.6f}}}",
            rf"\newcommand{{\ps{tag}ymin}}{{{im.min():.6f}}}",
            rf"\newcommand{{\ps{tag}ymax}}{{{im.max():.6f}}}",
            rf"\newcommand{{\ps{tag}zmin}}{{{Z.min():.4f}}}",
            rf"\newcommand{{\ps{tag}zmax}}{{{Z.max():.4f}}}",
        ]
        # The tracked mode ships as a one-row table rather than as inline
        # coordinates, so the figure never hard-codes a frequency.
        (OUT / f"pseudo_mark_{label}.dat").write_text(
            "re im\n"
            f"{p['omega_used'][0]:.10f} {p['omega_used'][1]:.10f}\n")
        if not any("psGrid" in f for f in facts):
            facts.append(rf"\newcommand{{\psGrid}}{{{re.size}}}")

    # The tracked pair at the tightest interaction.  Only the members that
    # actually fall inside the stored grid are written: a marker outside the
    # plotted window would silently stretch the axis and misrepresent the
    # field's extent.  Members left out are reported here and in the caption.
    cand = json.loads((ROOT / "results" / "all_collision_candidates.json").read_text())
    c = min(cand["candidates"], key=lambda x: x["gap"])
    sm = maps["tightest_interaction"]
    rlo, rhi = min(sm["re"]), max(sm["re"])
    ilo, ihi = min(sm["im"]), max(sm["im"])

    inside, outside = [], []
    for name in ("omega_a", "omega_b"):
        x, y = c[name]
        (inside if (rlo <= x <= rhi and ilo <= y <= ihi) else outside).append(
            (name, x, y))
    (OUT / "pseudo_pair.dat").write_text(
        "re im\n" + "\n".join(f"{x:.10f} {y:.10f}" for _, x, y in inside) + "\n")
    (OUT / "pseudo_pair.json").write_text(json.dumps(
        {"gap": c["gap"], "inside": inside, "outside": outside}, indent=1))
    print(f"pseudo_pair.dat  gap = {c['gap']:.6f}  "
          f"in window: {[n for n, _, _ in inside]}  "
          f"outside: {[(n, round(y, 4)) for n, _, y in outside]}")

    facts.append(rf"\newcommand{{\psGap}}{{{c['gap']:.4f}}}")
    facts.append(rf"\newcommand{{\psPartnerIm}}{{{outside[0][2]:.3f}}}"
                 if outside else r"\newcommand{\psPartnerIm}{}")
    (OUT / "pseudo_facts.tex").write_text("\n".join(facts) + "\n")


SECTORS = [(0, 0), (1, 1), (-1, -1), (1, 0), (2, 0), (2, 1), (2, 2)]
NRE, NIM = 180, 136
RE_LO, RE_HI, IM_LO, IM_HI = 0.0, 5.4, -3.4, 0.0
# Integer block replication: adds print resolution, never data.  At 18 the
# 180 x 136 field ships as 3240 x 2448 px, so at the 13.9 cm axis width of
# Fig. 1 every cell edge lands on a pixel edge at just under 600 dpi.
UPSCALE = 18

# magma control points, matching the colormap declared in fig_atlas.tex
EMBER = [(0.001, 0.000, 0.014), (0.232, 0.060, 0.438), (0.554, 0.161, 0.506),
         (0.868, 0.288, 0.409), (0.988, 0.553, 0.348), (0.987, 0.991, 0.750)]


def export_atlas() -> None:
    """Occupancy of the whole stored atlas in the complex-frequency plane.

    Every stored entry is counted; nothing is smoothed or interpolated.  The
    cell value is log10(count+1), so an empty cell is exactly zero and the
    filamentary continuation traces stay visible next to the dense ladders.
    """
    re_all, im_all, over = [], [], []
    for m1, m2 in SECTORS:
        p = ROOT / "data" / "full_branch_atlas" / f"sector_{m1}_{m2}.json"
        for q in json.loads(p.read_text())["points"]:
            re_all.append(q["omega_re"])
            im_all.append(q["omega_im"])
            over.append(q["overtone"])
    re_all = np.array(re_all)
    im_all = np.array(im_all)
    over = np.array(over)

    H, xe, ye = np.histogram2d(re_all, im_all, bins=[NRE, NIM],
                               range=[[RE_LO, RE_HI], [IM_LO, IM_HI]])
    Z = np.log10(H.T + 1.0)
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])

    # The field is genuinely pixel data, so it ships as an image: a grid this
    # size exhausts TeX's memory as vector cells, and every axis, label and
    # annotation around it is still typeset by PGFPlots in the document font.
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.image as mpimg

    cmap = LinearSegmentedColormap.from_list("emberdark", EMBER)
    norm = Z / Z.max()
    rgb = (cmap(norm)[:, :, :3] * 255).astype(np.uint8)
    rgb = np.repeat(np.repeat(rgb, UPSCALE, axis=0), UPSCALE, axis=1)
    rgb = rgb[::-1]                      # PNG rows run top-down
    mpimg.imsave(OUT / "atlas_density.png", rgb)

    bands = {int(n): float(np.median(im_all[over == n])) for n in np.unique(over)}
    meta = {
        "n_entries": int(re_all.size),
        "n_branches": len({(m1, m2, q["ell"], q["overtone"])
                           for m1, m2 in SECTORS
                           for q in json.loads(
                               (ROOT / "data" / "full_branch_atlas" /
                                f"sector_{m1}_{m2}.json").read_text())["points"]}),
        "grid": [NRE, NIM],
        "re_range": [RE_LO, RE_HI], "im_range": [IM_LO, IM_HI],
        "logn_max": float(Z.max()), "max_count": int(H.max()), "upscale": UPSCALE,
        "occupied_fraction": float((H > 0).mean()),
        "overtone_band_median_im": bands,
    }
    (OUT / "atlas_density.json").write_text(json.dumps(meta, indent=1))

    # Caption numbers as macros, so the prose cannot drift from the grid.
    n = meta["n_entries"]
    (OUT / "atlas_facts.tex").write_text("\n".join([
        rf"\newcommand{{\atlasN}}{{{n // 1000}{{,}}{n % 1000:03d}}}",
        rf"\newcommand{{\atlasBranches}}{{{meta['n_branches']}}}",
        rf"\newcommand{{\atlasGridRe}}{{{NRE}}}",
        rf"\newcommand{{\atlasGridIm}}{{{NIM}}}",
        rf"\newcommand{{\atlasOccupied}}{{{100 * meta['occupied_fraction']:.0f}}}",
        rf"\newcommand{{\atlasMetaMax}}{{{meta['logn_max']:.4f}}}",
    ]) + "\n")
    print(f"atlas_density.dat  {NRE}x{NIM} cells  {meta['n_entries']:,} entries  "
          f"{meta['n_branches']} branches  max count {meta['max_count']}")
    print(f"  overtone band medians (Im w): "
          f"{ {k: round(v, 3) for k, v in bands.items()} }")


def report_crossing() -> None:
    """Exact numbers for the equal-spin coincidence quoted in the text."""
    out = {}
    for m1, m2 in [(1, 1), (2, 0)]:
        rows = [p for p in sector_points(m1, m2)
                if p["ell"] == 2 and p["overtone"] == 0
                and abs(p["mu"] - 0.6) < 1e-9 and abs(p["s"] - 0.20) < 1e-9]
        at = {round(p["delta"], 6): p["omega_re"] for p in rows}
        out[f"({m1},{m2})"] = {"delta=0": at[0.0], "delta=0.2": at[0.2],
                               "delta=-0.2": at[-0.2]}
    (OUT / "crossing.json").write_text(json.dumps(out, indent=1))
    print("crossing:", json.dumps(out))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    export_continuation()
    export_pseudospectra()
    export_atlas()
    report_crossing()


if __name__ == "__main__":
    main()

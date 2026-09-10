"""Emit the numeric TMLR tables as LaTeX, computed from the repository record.

Table 1 comes from ``mp5d.geometry.MPGeometry``; Table 4 from the stored
collision-search result.  Nothing is typed by hand, so a table cannot drift
away from the code or the atlas it describes.
"""

from __future__ import annotations

import collections
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from mp5d.geometry.metric import MPGeometry, sdelta_to_ab  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "tmlr"

ROWS = [(0.00, 0.00), (0.20, 0.00), (0.42, 0.00),
        (0.34, 0.20), (0.42, 0.20), (0.42, -0.20)]


HEAD_H = r"""\begin{tabular}{
  S[table-format=1.2]
  S[table-format=+1.2]
  S[table-format=1.4]
  S[table-format=1.5]
  S[table-format=1.4]
  S[table-format=1.4]
}
\toprule
\thead{$s$} & \thead{$\delta$} & \thead{$r_+$} & \thead{$T_H$} &
\thead{$\Omega_a$} & \thead{$\Omega_b$} \\
\midrule"""

HEAD_C = r"""\begin{tabular}{lc}
\toprule
\thead{Quantity} & \thead{Value} \\
\midrule"""

FOOT = "\\bottomrule\n\\end{tabular}"


def table_horizons() -> str:
    body = []
    for s, d in ROWS:
        a, b = sdelta_to_ab(s, d)
        g = MPGeometry(a=a, b=b)
        body.append(
            f"  {s:.2f} & {d:+.2f} & {g.r_plus:.4f} & {g.T_H:.5f} & "
            f"{g.Omega_a:.4f} & {g.Omega_b:.4f} \\\\"
        )
    return "\n".join([HEAD_H, *body, FOOT])


def table_census() -> str:
    cand = json.loads((ROOT / "results" / "all_collision_candidates.json").read_text())
    atlas = {}
    for m1, m2 in [(0, 0), (1, 1), (-1, -1), (1, 0), (2, 0), (2, 1), (2, 2)]:
        p = ROOT / "data" / "full_branch_atlas" / f"sector_{m1}_{m2}.json"
        d = json.loads(p.read_text())
        atlas[(m1, m2)] = d
    n_pts = sum(d["n_points"] for d in atlas.values())
    max_res = max(p["residual"] for d in atlas.values() for p in d["points"])

    # A continuation failure is recorded with its reason; the two reasons mean
    # different things and are reported separately, never pooled.
    reasons = collections.Counter(f["reason"] for d in atlas.values()
                                  for f in d["failures"])

    mant, expo = f"{max_res:.1e}".split("e")
    res = f"${mant}\\times10^{{{int(expo)}}}$"

    rows = [
        ("Stored mode entries", f"{n_pts:,}"),
        ("Multi-branch parameter points", f"{cand['n_parameter_points']:,}"),
        ("Independent azimuthal sectors", f"{len(atlas)}"),
        ("Angular branches per sector", "3"),
        ("Overtones", "$N=0,1,2,3$"),
        ("Maximum stored continued-fraction residual", res),
        ("Predictor misses", f"{reasons['predictor_miss']}"),
        ("Nonconvergences", f"{reasons['not_converged']}"),
        ("Branch-label collapses excluded", f"{cand['n_collapsed_labels']}"),
    ]
    body = [f"  {k} & {v} \\\\" for k, v in rows]
    return "\n".join([HEAD_C, *body, FOOT])


def table_calibration() -> str:
    """Diagnostic calibration, run through the project's own routines.

    This replaces the former calibration figure panels.  Both blocks are
    computed here at build time on exactly solvable models, so the numbers
    cannot drift from the code that produced the atlas diagnostics.
    """
    import numpy as np
    from mp5d.exceptional.diagnostics import root_separation
    from mp5d.exceptional.verification import fit_puiseux, roots_in_disc

    ts = np.geomspace(1e-8, 1e-4, 12)

    def split(fn):
        return np.array([abs(np.subtract(*roots_in_disc(fn(t), 0.0,
                                                        radius=0.05)[:2]))
                         for t in ts])

    ep = split(lambda t: (lambda z, tt=t: z**2 - tt))
    lin = 3.0 * ts
    floor = 1e-3
    avoid = 2.0 * np.sqrt(ts**2 + floor**2)

    laws = [
        (r"square root, $\omega_\pm=\pm\sqrt{t}$", fit_puiseux(ts, ep).exponent,
         "EP2"),
        (r"linear, $\omega_\pm=\pm 3t$", fit_puiseux(ts, lin).exponent,
         "ordinary crossing"),
        (r"floored, $2\sqrt{t^2+g^2}$", fit_puiseux(ts, avoid).exponent,
         "avoided crossing"),
    ]

    kk = 3.0
    rows_n = []
    for dd in (1e-4, 1e-2, 3e-1):
        def f(z, D=dd):
            return z * (z - D)

        def fs(z, D=dd):
            return 1e7 * np.exp(kk * z) * z * (z - D)

        rad = min(0.4, 10 * dd)
        meas = (root_separation(fs, 0.0, radius=rad)["separation"]
                / root_separation(f, 0.0, radius=rad)["separation"])
        pred = abs(1.0 / (1.0 - kk * dd))
        rr = min(0.4, 3 * dd)
        n0 = len(roots_in_disc(f, 0.0, radius=rr))
        n1 = len(roots_in_disc(fs, 0.0, radius=rr))
        assert n0 == n1 == 2, (dd, n0, n1)
        mant, expo = f"{dd:.0e}".split("e")
        rows_n.append((f"${mant}\\times10^{{{int(expo)}}}$",
                       f"{meas:.3f}", f"{pred:.3f}", f"{n0}"))

    # Two stacked panels.  Both are set in tabular* at the same explicit
    # width so their rules line up: two booktabs blocks of different natural
    # widths sitting one above the other read as a typesetting accident.
    width = r"\begin{tabular*}{0.88\linewidth}{@{\extracolsep{\fill}}"

    out = [width + r"l c l@{}}",
           r"\toprule",
           r"\multicolumn{3}{@{}l}{\emph{Local splitting laws:} fitted Puiseux "
           r"exponent $p$ over $t\in[10^{-8},10^{-4}]$}\\",
           r"\midrule",
           r"\thead{Model} & \thead{$p$} & \thead{Classification}\\",
           r"\midrule"]
    out += [f"  {n} & {p:.3f} & {c} \\\\" for n, p, c in laws]
    out += [r"\bottomrule", r"\end{tabular*}", r"\\[1.4ex]",
            width + r"l c c c@{}}",
            r"\toprule",
            r"\multicolumn{4}{@{}l}{\emph{Normalization sensitivity} under "
            r"$F\mapsto gF$, $g=10^{7}e^{3\omega}$}\\",
            r"\midrule",
            r"\thead{$d$} & \thead{measured distortion} & "
            r"\thead{$|1+(g'/g)d|^{-1}$} & \thead{zero count}\\",
            r"\midrule"]
    out += [f"  {a} & {b} & {c} & {e} \\\\" for a, b, c, e in rows_n]
    out += [r"\bottomrule", r"\end{tabular*}"]
    return "\n".join(out)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tab_horizons.tex").write_text(table_horizons() + "\n")
    (OUT / "tab_census.tex").write_text(table_census() + "\n")
    (OUT / "tab_calibration.tex").write_text(table_calibration() + "\n")
    print(table_horizons())
    print()
    print(table_census())
    print()
    print(table_calibration())


if __name__ == "__main__":
    main()

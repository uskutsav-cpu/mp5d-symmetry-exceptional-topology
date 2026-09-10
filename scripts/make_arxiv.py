"""Derive the arXiv preprint tree from the TMLR submission tree.

The TMLR sources stay the single point of truth.  This script copies them into
``manuscript/arxiv/`` and applies one explicit, enumerated list of edits: turn
on the style file's ``preprint`` option, put the real author block in, and
reword the three places that only make sense inside an anonymized review.
Nothing is hand-edited on the arXiv side, so the two versions cannot drift.

Every substitution below is checked: if a pattern stops matching because the
TMLR source changed under it, the script fails instead of silently emitting a
preprint that still says "Anonymous authors".
"""

from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "manuscript" / "tmlr"
DST = ROOT / "manuscript" / "arxiv"

AUTHOR = "Utsav Sunil Kumar"
EMAIL = "uskutsav@gmail.com"
AFFIL = "Heritage High School, Frisco, Texas, USA"

# (description, pattern, replacement).  Patterns are literal strings.
EDITS: list[tuple[str, str, str]] = [
    (
        "style option: anonymous review -> preprint",
        "\\usepackage{tmlr}\n",
        # The preprint option empties the running head but leaves fancyhdr's
        # rule behind it, which on arXiv reads as a rule drawn by mistake.
        "\\usepackage[preprint]{tmlr}\n"
        "\\renewcommand{\\headrulewidth}{0pt}\n",
    ),
    (
        "author block",
        "\\author{\\name Anonymous Authors \\email anonymous@example.org \\\\\n"
        "        \\addr Anonymous Institution}",
        f"\\author{{\\name {AUTHOR} \\email {EMAIL} \\\\\n"
        f"        \\addr {AFFIL}}}",
    ),
    (
        "PDF author metadata",
        "  pdfauthor={},",
        f"  pdfauthor={{{AUTHOR}}},",
    ),
    (
        "disclosure footnote: single author",
        "consistency checking, and figure-layout refinement. The authors independently\n"
        "reviewed and approved all scientific claims, mathematical arguments, numerical\n"
        "results, references, figures, and text and take full responsibility for the\n"
        "submitted work.",
        "consistency checking, and figure-layout refinement. The author independently\n"
        "reviewed and approved all scientific claims, mathematical arguments, numerical\n"
        "results, references, figures, and text and takes full responsibility for this\n"
        "work.",
    ),
    (
        "data availability: no anonymized review process",
        "The numerical continuation atlas and solver materials supporting the\n"
        "finite-atlas claims are available from the authors upon reasonable request\n"
        "through the anonymized review process. They are not publicly archived at the\n"
        "time of submission. The anonymized Supplementary Material contains the solver\n"
        "inventory, atlas census, verification gate, symmetry checks, and threshold\n"
        "summary.",
        "The numerical continuation atlas and solver materials supporting the\n"
        "finite-atlas claims are available from the author upon reasonable request.\n"
        "They are not publicly archived at the time of posting. The Supplementary\n"
        "Material contains the solver inventory, atlas census, verification gate,\n"
        "symmetry checks, and threshold summary.",
    ),
]

# Everything the document actually reads.  Listing them explicitly keeps stale
# by-products of retired figures out of a public source package.
CARRY = [
    "main.tex", "refs.bib", "tmlr.sty", "tmlr.bst",
    "fig_atlas.tex", "fig_pseudospectrum.tex",
    "tab_horizons.tex", "tab_census.tex", "tab_calibration.tex",
]
CARRY_DATA = [
    "atlas_facts.tex", "pseudo_facts.tex", "atlas_density.png",
    "pseudo_ordinary_control.dat", "pseudo_tightest_interaction.dat",
    "pseudocont_ordinary_control.dat", "pseudocont_tightest_interaction.dat",
    "pseudo_mark_ordinary_control.dat", "pseudo_mark_tightest_interaction.dat",
    "pseudo_pair.dat",
    # provenance for the numbers the figures show; tiny, and worth publishing
    "atlas_density.json", "pseudo_ordinary_control.json",
    "pseudo_tightest_interaction.json", "pseudo_pair.json",
]


def check_inputs(tex: str) -> None:
    """Fail if the document reads a file the package would not ship."""
    read = set(re.findall(r"\\(?:input|include)\{([^}]+)\}", tex))
    read |= set(re.findall(r"\{(data/[A-Za-z0-9_./-]+)\}", tex))
    shipped = {p.removesuffix(".tex") for p in CARRY}
    shipped |= {f"data/{p}" for p in CARRY_DATA}
    shipped |= {f"data/{p.removesuffix('.tex')}" for p in CARRY_DATA}
    missing = {r for r in read if r not in shipped and not r.startswith("refs")}
    if missing:
        raise SystemExit(f"document reads files the package omits: {sorted(missing)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DST))
    args = ap.parse_args()
    dst = pathlib.Path(args.out)

    tex = (SRC / "main.tex").read_text()
    for what, pat, rep in EDITS:
        if tex.count(pat) != 1:
            raise SystemExit(
                f"cannot apply edit '{what}': pattern matched "
                f"{tex.count(pat)} times, expected exactly 1"
            )
        tex = tex.replace(pat, rep)
        print(f"  applied: {what}")

    for bad in ("Anonymous", "anonymized", "anonymous"):
        if bad in tex:
            raise SystemExit(f"'{bad}' still present in the preprint source")

    # Figures read data/... relative to the source they sit in, so the check
    # runs over the whole set of files that get concatenated at build time.
    whole = tex + "".join((SRC / f).read_text() for f in CARRY
                          if f.startswith(("fig_", "tab_")))
    check_inputs(whole)

    if dst.exists():
        shutil.rmtree(dst)
    (dst / "data").mkdir(parents=True)
    for f in CARRY:
        if f != "main.tex":
            shutil.copy2(SRC / f, dst / f)
    for f in CARRY_DATA:
        shutil.copy2(SRC / "data" / f, dst / "data" / f)
    (dst / "main.tex").write_text(tex)

    n = sum(1 for _ in dst.rglob("*") if _.is_file())
    size = sum(p.stat().st_size for p in dst.rglob("*") if p.is_file())
    print(f"wrote {dst} ({n} files, {size/1024:.0f} KiB)")
    print("next: build once here to produce main.bbl (arXiv does not run "
          "BibTeX), then `python scripts/pack_arxiv.py`")
    return 0


if __name__ == "__main__":
    sys.exit(main())

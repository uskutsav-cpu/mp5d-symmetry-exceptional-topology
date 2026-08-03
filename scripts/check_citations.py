"""Verify every \\cite key in the manuscript resolves to a bibliography entry.

LaTeX allows a citation list to be broken across lines with a trailing ``%``,
so the comment-continuations must be stripped before splitting on commas --
otherwise the checker invents missing keys that compile perfectly well.
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1] / "manuscript"


def cited_keys(tex: str) -> set[str]:
    tex = re.sub(r"%\s*\n\s*", "", tex)          # join %-continued lines
    tex = re.sub(r"(?m)^\s*%.*$", "", tex)        # drop whole-line comments
    keys: set[str] = set()
    for group in re.findall(r"\\cite[a-zA-Z]*\s*(?:\[[^\]]*\]\s*)*\{([^}]*)\}", tex):
        keys |= {k.strip() for k in group.split(",") if k.strip()}
    return keys


def main() -> int:
    bib = (ROOT / "references.bib").read_text()
    have = set(re.findall(r"@\w+\{([^,]+),", bib))
    bad = 0
    for name in ("main.tex", "supplement.tex"):
        p = ROOT / name
        if not p.exists():
            continue
        cited = cited_keys(p.read_text())
        missing = sorted(cited - have)
        print(f"{name}: {len(cited)} keys cited, {len(missing)} unresolved")
        for k in missing:
            print(f"  UNRESOLVED: {k}")
        bad += len(missing)
    print(f"bibliography: {len(have)} entries")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

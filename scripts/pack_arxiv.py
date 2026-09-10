"""Pack manuscript/arxiv into the tarball arXiv wants uploaded.

Two arXiv rules drive what goes in and what stays out:

* arXiv does **not** run BibTeX, so ``main.bbl`` must be shipped.  Build the
  directory once (``tectonic -X compile main.tex --keep-intermediates``)
  before packing, or this script refuses to run.
* a ``.pdf`` sitting next to the sources makes arXiv treat the upload as a
  PDF-only submission and skip the TeX build entirely, so build products are
  excluded rather than merely ignored.

Paths are stored flat (no wrapping directory) except for ``data/``, which the
figures reference by relative path.
"""

from __future__ import annotations

import pathlib
import sys
import tarfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "manuscript" / "arxiv"
OUT = ROOT / "manuscript" / "arxiv_upload.tar.gz"

DROP_SUFFIX = {".pdf", ".log", ".aux", ".out", ".blg", ".synctex", ".gz"}
DROP_NAME = {".DS_Store"}


def main() -> int:
    if not SRC.exists():
        raise SystemExit(f"{SRC} does not exist; run scripts/make_arxiv.py first")
    if not (SRC / "main.bbl").exists():
        raise SystemExit(
            "main.bbl is missing. arXiv does not run BibTeX, so the bibliography "
            "would come out empty. Build once in manuscript/arxiv with\n"
            "  tectonic -X compile main.tex --keep-intermediates"
        )

    members = sorted(
        p for p in SRC.rglob("*")
        if p.is_file()
        and p.suffix not in DROP_SUFFIX
        and p.name not in DROP_NAME
    )
    if not any(m.name == "main.tex" for m in members):
        raise SystemExit("main.tex missing from the package")

    OUT.unlink(missing_ok=True)
    with tarfile.open(OUT, "w:gz") as tar:
        for m in members:
            tar.add(m, arcname=str(m.relative_to(SRC)))

    print(f"wrote {OUT}  ({len(members)} files, {OUT.stat().st_size/1024:.0f} KiB)")
    for m in members:
        print(f"  {m.relative_to(SRC)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

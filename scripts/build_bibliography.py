"""Generate ``manuscript/references.bib`` from verified metadata only.

Sources of truth, in order:

1. ``bibliography/verified_metadata.json`` -- the authoritative arXiv record for
   every entry carrying an arXiv identifier, fetched by
   ``scripts/verify_bibliography.py`` (exact title, full author list, year, and
   the DOI and journal reference where arXiv holds them).
2. ``NON_ARXIV`` below -- entries without an arXiv identifier, each resolved
   against CrossRef.

Nothing is emitted that was not obtained from one of these.  In particular the
generated file contains no provenance notes, no "metadata not verified"
disclaimers and no internal commentary: a bibliography is part of the
manuscript, not a laboratory notebook.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]

# Entries with no arXiv identifier, each resolved against the CrossRef API.
# Kato's monograph is the Springer "Classics in Mathematics" reprint; CrossRef
# indexes only its individual chapters, so it is cited as a book, which needs no
# DOI.
NON_ARXIV: dict[str, dict] = {
    "A01": dict(type="article", author=["Myers, R. C.", "Perry, M. J."],
                title="Black holes in higher dimensional space-times",
                journal="Annals of Physics", volume="172", pages="304--347",
                year="1986", doi="10.1016/0003-4916(86)90186-7"),
    "B32": dict(type="article", author=["Wang, Y.", "Wu, S."],
                title="Coexistence of spectrally stable and unstable modes in "
                      "black hole ringdowns",
                journal="Physical Review D", volume="114", year="2026",
                doi="10.1103/srxw-wsvq"),
    "C34": dict(type="article", author=["Nollert, H.-P.", "Price, R. H."],
                title="Quantifying excitations of quasinormal mode systems",
                journal="Journal of Mathematical Physics", volume="40",
                pages="980--1010", year="1999", doi="10.1063/1.532698"),
    "C40": dict(type="article", author=["Gasper\\'in, E.", "Jaramillo, J. L."],
                title="Energy scales and black hole pseudospectra: the "
                      "structural role of the scalar product",
                journal="Classical and Quantum Gravity", volume="39",
                pages="115010", year="2022",
                doi="10.1088/1361-6382/ac5054"),
    "D51": dict(type="article", author=["Leaver, E. W."],
                title="An analytic representation for the quasi-normal modes "
                      "of Kerr black holes",
                journal="Proceedings of the Royal Society A", volume="402",
                pages="285--298", year="1985",
                doi="10.1098/rspa.1985.0119"),
    "D52": dict(type="article", author=["Leaver, E. W."],
                title="Spectral decomposition of the perturbation response of "
                      "the Schwarzschild geometry",
                journal="Physical Review D", volume="34", pages="384--408",
                year="1986", doi="10.1103/PhysRevD.34.384"),
    "D53": dict(type="article", author=["Detweiler, S."],
                title="Klein-Gordon equation and rotating black holes",
                journal="Physical Review D", volume="22", pages="2323--2326",
                year="1980", doi="10.1103/PhysRevD.22.2323"),
    "E71": dict(type="book", author=["Kato, T."],
                title="Perturbation Theory for Linear Operators",
                publisher="Springer", series="Classics in Mathematics",
                year="1995"),
    "E72": dict(type="book", author=["Trefethen, L. N.", "Embree, M."],
                title="Spectra and Pseudospectra: The Behavior of Nonnormal "
                      "Matrices and Operators",
                publisher="Princeton University Press", year="2005",
                doi="10.1515/9780691213101"),
    "E73": dict(type="book",
                author=["Moore, R. E.", "Kearfott, R. B.", "Cloud, M. J."],
                title="Introduction to Interval Analysis",
                publisher="Society for Industrial and Applied Mathematics",
                year="2009", doi="10.1137/1.9780898717716"),
    "E74": dict(type="book", author=["Tucker, W."],
                title="Validated Numerics: A Short Introduction to Rigorous "
                      "Computations",
                publisher="Princeton University Press", year="2011",
                doi="10.1515/9781400838974"),
}

# Cited but absent from the supplied corpus; both verified.
EXTRA: dict[str, dict] = {
    "Matyjasek2021": dict(type="article", author=["Matyjasek, J."],
                          title="Quasinormal modes of black holes in "
                                "Einstein-Gauss-Bonnet and other theories: "
                                "the Schwarzschild-Tangherlini case",
                          eprint="2107.04815", year="2021"),
    "KravanjaVanBarel2000": dict(
        type="book", author=["Kravanja, P.", "Van Barel, M."],
        title="Computing the Zeros of Analytic Functions",
        publisher="Springer", series="Lecture Notes in Mathematics 1727",
        year="2000", doi="10.1007/BFb0103927"),
}


def bibkey(row: dict) -> str:
    first = re.split(r"[;,]", row["authors"])[0].strip()
    first = re.sub(r"[^A-Za-z]", "", first) or "Anon"
    ident = row.get("identifier", "")
    m = re.search(r"\b(\d{2})(\d{2})\.\d{4,5}", ident)
    if m:
        yr = f"20{int(m.group(1)):02d}"
    else:
        m = re.search(r"/(\d{2})(\d{2})(\d{3})\b", ident)
        if m:
            yy = int(m.group(1))
            yr = f"{'19' if yy >= 90 else '20'}{yy:02d}"
        else:
            m2 = re.search(r"\b(19|20)\d{2}\b", ident)
            yr = m2.group(0) if m2 else ""
    return f"{first}{yr}{row['id']}"


def emit(key: str, rec: dict) -> str:
    lines = [f"@{rec.get('type', 'article')}{{{key},"]
    fields = [("author", " and ".join(rec["author"])),
              ("title", "{" + rec["title"] + "}")]
    for f in ("journal", "booktitle", "publisher", "series", "volume",
              "pages", "year"):
        if rec.get(f):
            fields.append((f, rec[f]))
    if rec.get("eprint"):
        fields.append(("eprint", rec["eprint"]))
        fields.append(("archivePrefix", "arXiv"))
    if rec.get("doi"):
        fields.append(("doi", rec["doi"]))
    body = ",\n".join(f"  {k:<14}= {{{v}}}" for k, v in fields)
    return lines[0] + "\n" + body + "\n}\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="bibliography/mandatory_sources.csv")
    ap.add_argument("--verified", default="bibliography/verified_metadata.json")
    ap.add_argument("--out", default="manuscript/references.bib")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(ROOT / args.csv)))
    ver = json.loads((ROOT / args.verified).read_text())["records"]

    out, keymap, incomplete = [], [], []
    for r in rows:
        key = bibkey(r)
        ident = (r.get("identifier") or "").strip()
        if ident.startswith("arXiv:"):
            aid = ident.split(":", 1)[1]
            v = ver.get(aid)
            if not v:
                incomplete.append(r["id"])
                continue
            jr = v.get("journal_ref") or ""
            m = re.search(r"([A-Za-z .]+?)\s*\.?\s*(\d+)[,:]\s*([0-9\-]+)", jr)
            rec = dict(type="article", author=v["authors"], title=v["title"],
                       year=v["year"], eprint=aid, doi=v.get("doi"))
            if m:
                rec["journal"] = m.group(1).strip().rstrip(".")
                rec["volume"] = m.group(2)
                rec["pages"] = m.group(3)
            out.append(emit(key, rec))
        elif r["id"] in NON_ARXIV:
            out.append(emit(key, NON_ARXIV[r["id"]]))
        else:
            incomplete.append(r["id"])
        keymap.append(f"{r['id']}\t{key}\t{r['short_title']}")

    for k, rec in EXTRA.items():
        out.append(emit(k, rec))

    (ROOT / args.out).write_text("\n".join(out))
    (ROOT / "manuscript" / "citekeys.txt").write_text("\n".join(keymap) + "\n")
    print(f"wrote {args.out}: {len(out)} entries")
    if incomplete:
        print(f"INCOMPLETE (not emitted): {incomplete}")
        return 1
    print("every emitted entry carries verified metadata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

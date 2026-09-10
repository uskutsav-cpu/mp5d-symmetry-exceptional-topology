"""Verify every arXiv entry in the corpus against the arXiv API.

The corpus supplies short titles and identifiers only.  A submitted manuscript
needs full author lists, exact titles, journal references and DOIs, and none of
those may be invented.  This script fetches the authoritative record for every
arXiv identifier in ``bibliography/mandatory_sources.csv`` and writes the result
to ``bibliography/verified_metadata.json``.

Entries without an arXiv identifier (books, and a few older journal articles)
cannot be verified this way and are listed in the output under
``unverifiable_here`` so that they can be completed by hand rather than silently
shipped with placeholder text.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parents[1]
ATOM = "{http://www.w3.org/2005/Atom}"
ARX = "{http://arxiv.org/schemas/atom}"


def fetch(ids: list[str], timeout: int = 40) -> str:
    url = ("https://export.arxiv.org/api/query?id_list="
           + ",".join(ids) + f"&max_results={len(ids)}")
    out = subprocess.run(["curl", "-sS", "--max-time", str(timeout), url],
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"curl failed: {out.stderr[:200]}")
    return out.stdout


def parse(xml: str) -> dict[str, dict]:
    root = ET.fromstring(xml)
    recs: dict[str, dict] = {}
    for e in root.findall(f"{ATOM}entry"):
        raw_id = e.findtext(f"{ATOM}id") or ""
        m = re.search(r"abs/(.+?)(v\d+)?$", raw_id)
        if not m:
            continue
        key = m.group(1)
        authors = [a.findtext(f"{ATOM}name").strip()
                   for a in e.findall(f"{ATOM}author")
                   if a.findtext(f"{ATOM}name")]
        title = " ".join((e.findtext(f"{ATOM}title") or "").split())
        published = (e.findtext(f"{ATOM}published") or "")[:4]
        recs[key] = {
            "arxiv_id": key,
            "title": title,
            "authors": authors,
            "year": published,
            "doi": (e.findtext(f"{ARX}doi") or "").strip() or None,
            "journal_ref": " ".join(
                (e.findtext(f"{ARX}journal_ref") or "").split()) or None,
        }
    return recs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="bibliography/mandatory_sources.csv")
    ap.add_argument("--out", default="bibliography/verified_metadata.json")
    ap.add_argument("--batch", type=int, default=25)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(ROOT / args.csv)))
    arxiv_ids, no_id = [], []
    for r in rows:
        ident = (r.get("identifier") or "").strip()
        if ident.startswith("arXiv:"):
            arxiv_ids.append(ident.split(":", 1)[1])
        else:
            no_id.append({"id": r["id"], "authors": r["authors"],
                          "short_title": r["short_title"],
                          "identifier": ident})

    got: dict[str, dict] = {}
    for i in range(0, len(arxiv_ids), args.batch):
        chunk = arxiv_ids[i:i + args.batch]
        for attempt in range(3):
            try:
                got.update(parse(fetch(chunk)))
                break
            except Exception as exc:  # noqa: BLE001
                print(f"  batch {i}: attempt {attempt+1} failed: "
                      f"{str(exc)[:120]}", file=sys.stderr)
                time.sleep(4)
        print(f"  fetched {min(i+args.batch, len(arxiv_ids))}/{len(arxiv_ids)}",
              flush=True)
        time.sleep(3.2)  # arXiv asks for >3 s between requests

    missing = [a for a in arxiv_ids if a not in got]
    out = {
        "source": "arXiv API (export.arxiv.org), authoritative record",
        "n_requested": len(arxiv_ids),
        "n_verified": len(got),
        "missing_from_api": missing,
        "unverifiable_here": no_id,
        "records": got,
    }
    (ROOT / args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"verified {len(got)}/{len(arxiv_ids)} arXiv entries; "
          f"{len(no_id)} have no arXiv id; {len(missing)} not returned")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

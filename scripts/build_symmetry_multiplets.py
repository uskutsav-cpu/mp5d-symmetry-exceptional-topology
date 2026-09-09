#!/usr/bin/env python3
"""Emit the exact equal-spin U(2) multiplet decomposition to results/.

At delta = 0 the separated problem depends on (m1, m2) only through
m = m1 + m2, so the degenerate sets are the fixed-(l, m) blocks.  Each is a
single irreducible SU(2) multiplet of dimension l + 1.

This is an *analytic* artifact: no solver is involved, nothing here is a
numerical estimate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mp5d.angular import s3_multiplet  # noqa: E402


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "uncommitted"


def lock_hash() -> str:
    p = pathlib.Path(__file__).resolve().parents[1] / "requirements.lock"
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16] if p.exists() else "none"


def build(l_max: int) -> dict:
    multiplets = []
    for l in range(l_max + 1):
        states = s3_multiplet(l)
        assert len(states) == (l + 1) ** 2
        for m in range(-l, l + 1):
            if (m - l) % 2 != 0:
                continue
            block = sorted(t for t in states if t[1] + t[2] == m)
            assert len(block) == l + 1
            multiplets.append(
                {
                    "l": l,
                    "m": m,
                    "u1_charge": m,
                    "su2_irrep_dim": l + 1,
                    "degeneracy": len(block),
                    "states_n_m1_m2": [list(t) for t in block],
                    "diagonal_sector_members": [list(t) for t in block if t[1] == t[2]],
                    "defective_possible": False,
                    "reason": (
                        "members carry distinct (m1,m2), which label blocks that "
                        "are decoupled at every delta; degeneracy is semisimple"
                    ),
                }
            )
    return {
        "status": "complete",
        "kind": "analytic",
        "provenance": {
            "commit": git_commit(),
            "lock_sha256_16": lock_hash(),
            "solver": "none (exact representation theory)",
            "precision": "exact",
        },
        "definitions": {
            "l": "2n + |m1| + |m2|, S^3 angular momentum",
            "m": "m1 + m2, the U(1) charge of U(2) = (SU(2) x U(1))/Z2",
            "validity": "delta = 0 only; for delta != 0 these sets split",
        },
        "l_max": l_max,
        "records": multiplets,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--l-max", type=int, default=8)
    ap.add_argument(
        "--out",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1]
        / "results"
        / "symmetry_multiplets.json",
    )
    args = ap.parse_args()
    data = build(args.l_max)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print(f"wrote {args.out} : {len(data['records'])} multiplets up to l={args.l_max}")


if __name__ == "__main__":
    main()

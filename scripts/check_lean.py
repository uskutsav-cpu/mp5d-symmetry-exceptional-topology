#!/usr/bin/env python3
"""Compile the exact layer and write a provenance-bearing receipt; no fake PASS."""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mp5d_science.provenance import atomic_json, envelope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--development", action="store_true")
    args = parser.parse_args()
    if not shutil.which("lake"):
        payload = {
            "status": "NOT_BUILT",
            "reason": "lake executable is unavailable",
            "continuum_qnm_theorem": False,
        }
    else:
        result = subprocess.run(
            ["lake", "build"],
            cwd=ROOT / "lean",
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        payload = {
            "status": "BUILD_PASS" if result.returncode == 0 else "BUILD_FAILED",
            "returncode": result.returncode,
            "build_log": result.stdout,
            "scope": "ALGEBRAIC_OPERATOR_SYMMETRY_AND_PARITY_STATEMENTS",
            "continuum_qnm_theorem": False,
        }
    atomic_json(
        args.output,
        envelope("lean-kernel-build", payload, ROOT, allow_development=args.development),
    )
    return 0 if payload["status"] == "BUILD_PASS" else 2


if __name__ == "__main__":
    sys.exit(main())

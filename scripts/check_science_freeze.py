#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mp5d_science.claims import validate_freeze
from mp5d_science.provenance import atomic_json, read_json


def main():
    parser = argparse.ArgumentParser(
        description="Do not tag or publish unless every science gate passes."
    )
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--expected-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = validate_freeze(
        ROOT,
        args.artifact_dir,
        args.expected_commit,
        read_json(ROOT / "config/science/claims.json"),
    )
    atomic_json(args.output, result)
    for error in result["errors"]:
        print(error, file=sys.stderr)
    return 0 if result["status"] == "SCIENCE_FREEZE_PASS" else 2


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mp5d_science.provenance import check_result, read_json, source_snapshot


def main():
    parser = argparse.ArgumentParser(
        description="Require exact committed clean source, never an ancestor commit."
    )
    parser.add_argument("files", type=Path, nargs="+")
    parser.add_argument("--expected-commit", required=True)
    args = parser.parse_args()
    errors = []
    for path in args.files:
        try:
            errors.extend(
                f"{path}: {e}"
                for e in check_result(
                    read_json(path),
                    expected_commit=args.expected_commit,
                    expected_source_digest=source_snapshot(ROOT)["sha256"],
                    authoritative=True,
                )
            )
        except (OSError, ValueError, TypeError) as exc:
            errors.append(f"{path}: {exc}")
    for error in errors:
        print(error, file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env bash
# Full reproduction: core science plus every validation and diagnostic run.
#
#   ./scripts/reproduce_full_science.sh            # everything (several hours)
#   ./scripts/reproduce_full_science.sh --resume   # skip completed shards
set -euo pipefail
cd "$(dirname "$0")/.."

RESUME=""
for arg in "$@"; do
  case "$arg" in
    --resume) RESUME="--resume" ;;
    --quick|--full) ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

PY=.venv/bin/python
export PYTHONPATH=src

echo "==> full test suite"
$PY -m pytest -q

echo "==> core science"
./scripts/reproduce_core_science.sh

echo "==> long-lived / quasi-bound branches"
$PY scripts/track_long_lived.py

echo "==> near-extremal spectral classification"
$PY scripts/map_near_extremal_spectrum.py

echo "==> targeted pseudospectral analysis"
$PY scripts/pseudospectrum.py

echo "==> published benchmarks"
$PY scripts/first_qnm.py
$PY scripts/rotating_validation.py

echo "done."

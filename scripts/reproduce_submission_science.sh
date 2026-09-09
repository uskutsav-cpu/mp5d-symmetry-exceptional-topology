#!/usr/bin/env bash
# Never discard an error, omit slow tests, relabel old results, or tag on failure.
# Usage: bash scripts/reproduce_submission_science.sh /ABSOLUTE/EXTERNAL/output-dir
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [[ $# -ne 1 ]]; then echo 'Supply one external output directory.' >&2; exit 2; fi
OUT="$1"
python3 - "$ROOT" "$OUT" <<'PY'
from pathlib import Path
import sys
root,out=map(lambda p:Path(p).resolve(),sys.argv[1:])
if out==root or out.is_relative_to(root):raise SystemExit('Results must be external to the source checkout.')
out.mkdir(parents=True,exist_ok=True)
PY
COMMIT="$(git rev-parse HEAD)"
if [[ -n "$(git status --porcelain --untracked-files=all)" ]]; then
  echo 'Commit the reviewed source changes before authoritative reproduction.' >&2; exit 2
fi
python3 -m venv "$OUT/environment"
PY="$OUT/environment/bin/python"
"$PY" -m pip install --disable-pip-version-check -r requirements-science.lock -r requirements-arb.lock
export PYTHONPATH="$ROOT/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$PY" -m pip check
# Original tests plus every new test; no marker exclusion and no xfail override.
"$PY" -m pytest tests -q --strict-markers --junitxml="$OUT/all-tests.xml"
"$PY" scripts/check_claim_consistency.py
"$PY" scripts/run_benchmark_audit.py --output "$OUT/benchmarks.json"
# This deliberately fails on the current invalid historical N2 roots.
"$PY" scripts/run_candidate_audit.py --output "$OUT/candidates.json"
# Additional required artifacts are NOT presumed present or blessed from history.
"$PY" scripts/check_science_freeze.py --artifact-dir "$OUT" --expected-commit "$COMMIT" --output "$OUT/freeze.json"
echo "Every configured science gate passed at $COMMIT. No tag was created automatically."

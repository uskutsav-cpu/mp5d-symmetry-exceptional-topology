#!/usr/bin/env bash
# Reproduce the core scientific result from a clean clone.
#
#   ./scripts/reproduce_core_science.sh --quick   # analytical + machinery only (~5 min)
#   ./scripts/reproduce_core_science.sh           # + atlas + collision search (~1 h)
#
# The core result is the multi-sector bounded exclusion of EP2 together with the
# symmetry classification.  Everything here is deterministic: seeds are fixed and
# no step depends on wall-clock time.
set -euo pipefail
cd "$(dirname "$0")/.."

QUICK=0
for arg in "$@"; do
  case "$arg" in
    --quick) QUICK=1 ;;
    --full)  QUICK=0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

if [ ! -d .venv ]; then
  echo "==> creating environment"
  python3 -m venv .venv
  .venv/bin/pip install --quiet -r requirements.lock
fi
PY=.venv/bin/python
export PYTHONPATH=src

echo "==> analytical + machinery tests (symmetry group, no-go, EP detectors)"
$PY -m pytest -q tests/unit/test_ep_machinery_synthetic.py tests/symmetry/

if [ "$QUICK" = "1" ]; then
  echo "==> quick mode: skipping the atlas"
  exit 0
fi

echo "==> static seed table"
$PY scripts/build_static_seed_table.py --l-max 8 --overtones 4

echo "==> branch atlas (7 sectors, resumable)"
for sec in "0,0" "1,1" "1,0" "2,0" "2,1" "2,2"; do
  $PY scripts/build_branch_atlas.py --sector "$sec" --resume
done
$PY scripts/build_branch_atlas.py --sector=-1,-1 --resume

echo "==> collision search + bounded exclusion"
$PY scripts/collision_search.py

echo "==> effective models"
$PY scripts/effective_model.py

echo
echo "core result:"
$PY - <<'PYEOF'
import json
d = json.load(open("results/negative_regions.json"))
print("  min branch gap      :", d["min_branch_gap"])
print("  min root separation :", d["min_root_separation"])
print("  atlas points        :", d["n_atlas_points"])
print("  sectors             :", d["domain"]["sectors"])
PYEOF

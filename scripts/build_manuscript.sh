#!/usr/bin/env bash
# Build the manuscript, supplement and letters from source.
#
#   ./scripts/build_manuscript.sh            # figures + all documents
#   ./scripts/build_manuscript.sh --no-figs  # documents only
#
# Requires `tectonic` (self-contained LaTeX). Figures additionally require the
# plotting stack, which is kept OUT of requirements.lock so that neither the
# scientific results nor CI depend on a graphics library:
#
#   .venv/bin/pip install -r requirements-figures.lock
set -euo pipefail
cd "$(dirname "$0")/.."

FIGS=1
for arg in "$@"; do
  case "$arg" in
    --no-figs) FIGS=0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

if ! command -v tectonic >/dev/null 2>&1; then
  echo "tectonic not found. Install it, or compile manuscript/*.tex with any" >&2
  echo "LaTeX distribution providing revtex4-2." >&2
  exit 1
fi

if [ "$FIGS" = "1" ]; then
  if .venv/bin/python -c "import matplotlib" >/dev/null 2>&1; then
    echo "==> figures"
    PYTHONPATH=src .venv/bin/python scripts/make_figures.py
  else
    echo "==> figures SKIPPED (matplotlib absent; see requirements-figures.lock)"
    echo "    Existing PDFs in manuscript/figures/ will be used if present."
  fi
fi

mkdir -p manuscript/build
cd manuscript
for doc in main supplement cover_letter_prd cover_letter_prl response_template; do
  echo "==> $doc"
  tectonic -X compile "$doc.tex" --outdir build
done

echo
echo "built:"
ls -1 build/*.pdf

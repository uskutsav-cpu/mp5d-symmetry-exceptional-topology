# MP5D science-hardening update bundle

**Read `DELIVERY_REPORT.md` before applying. This is not a science-freeze release.**

Audited repository: `uskutsav-cpu/mp5d-symmetry-exceptional-topology`.
Baseline: `research/absolute-final`, commit
`72abf0bf68d41d07de0fbd5056c40b7a713a6cb5`.
The checklist corresponds to this research branch, not the older `main` branch.

The update adds implemented/tested numerical and validation code, plus a guarded
installer for the original APIs. The physical audit found unresolved problems
with all ten saved N2 frequencies; it does not substantiate the old minimum-gap
claim. Some parts of the twelve-item research request remain unfinished.

## Contents

- `overlay/`: source, tests, scripts, regression manifests, configuration,
  documentation, workflow and uncompiled Lean drafts.
- `APPLY_UPDATE.py`: exact-baseline, clean-checkout, new-branch installer.
- `PAYLOAD_SHA256.json`: hashes of every payload file and installer.
- `evidence/`: actual numerical receipts, test logs and machine-readable results.
- `test_bundle_installer.py`: safety tests using isolated synthetic Git fixtures.
- `DELIVERY_REPORT.md`: checklist status, executed checks and limitations.

The ZIP is an update bundle, **not a complete clone of the original repository**.
Remote writes were denied by GitHub. Nothing has been pushed or merged.

## Install into a new review checkout

On macOS Terminal, after downloading the ZIP into Downloads, execute this block.
It deliberately uses a new directory and stops on errors. Existing directories
are not overwritten. Review the dry-run output before executing the separate
apply command.

```bash
set -euo pipefail
cd "$HOME/Downloads"
unzip mp5d-science-hardening.zip
git clone --branch research/absolute-final --single-branch \
  https://github.com/uskutsav-cpu/mp5d-symmetry-exceptional-topology.git \
  mp5d-hardening-review
cd mp5d-hardening-review
git switch -c research/science-hardening-review \
  72abf0bf68d41d07de0fbd5056c40b7a713a6cb5
python3 ../mp5d-science-hardening/APPLY_UPDATE.py "$PWD" --dry-run
```

Apply only after that preflight succeeds:

```bash
python3 ../mp5d-science-hardening/APPLY_UPDATE.py "$PWD"
git diff --stat
git diff -- src/mp5d results/status.json README.md pyproject.toml
```

The installer validates payload hashes, repository origin, exact source commit,
clean working tree and a non-protected review branch. It preflights every source
edit. It never fetches, resets, switches branches, commits, pushes or tags. A
write failure restores prior file contents. Source files from another revision
are refused instead of patched by guesswork.

This update intentionally retires unsafe legacy evidence APIs and changes
convergence decisions. Old scripts relying on those APIs may now stop with an
explicit error. The complete patched upstream suite and original Ruff lint have
not been run here. Inspect those failures rather than bypassing or weakening
scientific gates.

## Verify the standalone new code

A Python 3.13 virtual environment is recommended; the executed environment was
CPython 3.13.5 on Linux. The locked numerical packages are recorded in
`overlay/requirements-science.lock`. macOS and other environments have not been
tested by this delivery.

```bash
cd "$HOME/Downloads/mp5d-science-hardening"
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r overlay/requirements-science.lock
PYTHONPATH=overlay/src python -m pytest -c overlay/pytest-science.ini overlay/tests/submission -q
python -m pytest test_bundle_installer.py -q
```

To reproduce the physical audit, run the commands in
`overlay/docs/SCIENCE_HARDENING.md`. A failed scientific audit writes its receipt
and exits with code 2. A passing software test suite does not erase that failure.

## Before any source freeze

Repair/independently identify the failed historical branches, establish all
required convergence ladders, regenerate the atlas and boundary search, execute
the higher-mode probe, complete all benchmark/convention checks, compile/review
the exact Lean layer, and run the complete upstream tests and CI from the
reviewed source commit. The full massive MP5D hyperboloidal solver has not been
implemented; quasiresonant work is explicitly excluded in the present scope.

The strict submission script requires a clean committed source checkout and an
output directory outside it. It currently must fail on the unresolved science
and missing artifacts. Do not hand-edit receipts, set `science_ready: true` to
bypass failures, create a freeze tag, or regenerate manuscript numbers from the
old JSONs.

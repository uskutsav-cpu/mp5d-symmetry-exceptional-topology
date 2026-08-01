# mp5d-symmetry-exceptional-topology

Does the enhanced `U(2)` symmetry of the equal-spin surface organize the
exceptional-set topology of the massive-scalar quasinormal-mode spectrum of the
asymptotically flat, doubly rotating five-dimensional Myers-Perry black hole?

Parameter space: `(s, δ, μ)` with `s = (a+b)/2`, `δ = (a−b)/2`, `μ` the scalar
mass. `δ = 0` is the equal-spin surface.

**Private research repository. Do not make public.**

## Status

See `results/status.json` and `docs/SESSION_HANDOFF.md`.

## Reading order

1. `docs/CONVENTIONS.md` — metric, separation, sign conventions. All of it is
   re-derived symbolically on every CI run by `tests/unit/test_separation.py`.
2. `docs/SYMMETRY_STRUCTURE.md` — the equal-spin `U(2)` analysis and the
   structural no-go that follows from it.
3. `docs/CLAIM_LEDGER.md` — every claim with its evidence level.

## Install

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install -e . --no-deps
.venv/bin/python -m pytest -q
```

`uv` is not available on the development machine; `requirements.lock` is a
`pip freeze` of the working environment and plays the role of `uv.lock`.

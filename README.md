# mp5d-symmetry-exceptional-topology

Active research code. Nothing here is a finished result.

**Question.** Does the enhanced `U(2)` symmetry of the equal-spin surface
organize the exceptional-set topology of the massive-scalar quasinormal-mode
spectrum of the asymptotically flat, doubly rotating five-dimensional
Myers–Perry black hole?

Parameter space: `(s, δ, μ)` with `s = (a+b)/2`, `δ = (a−b)/2`, and `μ` the
scalar-field mass. `δ = 0` is the equal-spin surface. The Myers–Perry mass
parameter is `M`; it is never written `μ`.

## Honest status

Read this before reading anything else in the repository.

| Component | State |
| --- | --- |
| Metric, conventions, separation of the massive Klein–Gordon equation | **Derived symbolically**, re-checked in exact rational arithmetic on every test run |
| Geometry (horizons, `Ω_a`, `Ω_b`, `κ`, `T_H`) | Implemented, tested |
| Angular sector | Implemented; exact Jacobi construction cross-checked against an independent finite-difference solver |
| Equal-spin `U(2)` multiplet structure | Derived exactly |
| Radial singular structure and characteristic exponents | Derived exactly |
| **Radial QNM eigenvalue solver** | see `results/status.json` |
| **Computed QNM frequencies** | see `results/status.json` |
| **Exceptional points** | **none verified, none claimed** |

`results/status.json` is the authoritative machine-readable status and is
updated at every milestone. Where this README and that file disagree, the JSON
is correct.

### What is *not* claimed

* No quasinormal frequency here should be treated as publication-grade until it
  appears in `results/baseline_validation.json` having passed the benchmark gate
  in `.claude/skills/mp5d-bulk-research/SCIENTIFIC_GATES.md`.
* No exceptional point has been verified. A close pair of frequencies is not an
  exceptional point, and this repository will not call one that.
* The cross-sector no-go argument (`docs/CROSS_SECTOR_EP_NO_GO.md`) is
  **provisional**. Its proof obligations are listed there and not all are
  closed. It is not a theorem.
* Bibliography entries in `bibliography/mandatory_sources.csv` are **indexed**.
  The `metadata_verified` column records which have actually been checked
  against the source; most have not.
* Numerical convergence is not a proof about the continuum operator. Every
  certification is labelled with what it certifies — a matrix, a truncated
  recurrence, an algebraic system, or the continuum operator — and the first
  three never imply the fourth.

This project is not advertised as a discovery, and no claim of novelty is made
against the prior literature indexed in `docs/PRIOR_RESULTS.md`.

### Stability

APIs, file layouts, and JSON schemas change without notice. This is a research
repository, not a library.

## Reading order

1. `docs/CONVENTIONS.md` — metric, separation, sign conventions. Pinned by
   `tests/unit/test_separation.py`.
2. `docs/SYMMETRY_STRUCTURE.md` — the equal-spin analysis and what it implies
   about where exceptional points can live.
3. `docs/CROSS_SECTOR_EP_NO_GO.md` — the provisional no-go and its open
   obligations.
4. `docs/CLAIM_LEDGER.md` — every claim with its evidence level.
5. `docs/RECOVERY_AUDIT.md`, `docs/FAILED_APPROACHES.md` — state verification
   and preserved negative results.

## Install and run

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install -e . --no-deps
PYTHONPATH=src .venv/bin/python -m pytest -q
```

`uv` and Julia are not installed on the development machine; `requirements.lock`
is a `pip freeze` and plays the role of `uv.lock`.

## License

MIT, first-party code only. No third-party source is vendored. The bibliography
contains metadata only — no copyrighted article text or PDFs.

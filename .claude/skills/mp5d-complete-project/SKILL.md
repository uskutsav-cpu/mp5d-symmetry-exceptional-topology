---
name: mp5d-complete-project
description: Durable operating rules for the MP5D symmetry/exceptional-topology research repository. Load when continuing scientific work on doubly rotating five-dimensional Myers-Perry massive-scalar QNMs, exceptional points, branch atlases, or the associated solvers.
---

# MP5D research operating rules

## Scientific target

Classify exceptional spectral structures **within fixed `(m1,m2)` sectors** of
the massive-scalar QNM problem for asymptotically flat, doubly rotating 5D
Myers-Perry black holes, with attention to the equal-spin surface and to
diagonal sectors `m1 = m2`.

The cross-sector channel is **closed** (`docs/CROSS_SECTOR_EP_NO_GO.md`): Jordan
chains of a direct sum are the union of the summands' chains, so a coincidence
of eigenvalues across distinct sectors cannot create defectiveness. Do not
re-open it. All defectiveness in this system is within-sector defectiveness.

## Non-negotiables

* **Verify state before acting.** Run `pytest -q` first. Never trust a handoff.
* **Never call a close frequency pair an exceptional point.** The EP2 gate is in
  `docs/EP_VERIFICATION_GATE.md`; every applicable condition must pass.
* **Never use `|dF/domega|` as a bound.** It is not invariant under rescaling of
  the spectral condition. Use the scale-free root separation `|a1/a2|`
  (`mp5d.exceptional.diagnostics.root_separation`) or the pairwise branch gap.
* **Never claim continuum certification from a matrix or recurrence truncation.**
  Use the labels in `STATUS_SCHEMA.md`.
* **Never label a branch by sorting frequencies.** Identity comes from
  continuation history plus predictor agreement.
* **Stay inside the validated solver domain** (`z_minus <= 0.20`, extremality
  `>= 0.05`) unless independently validating. `z_minus` is the *inner* horizon
  radius squared. Stepping outside it silently produced the refuted
  near-extremal candidate.
* **Preserve failures** in `docs/FAILED_APPROACHES.md` with the mechanism, not
  just the symptom.
* Work on a research branch, never `main`. Push every validated milestone.

## Conventions

`docs/CONVENTIONS.md` is authoritative and is re-derived symbolically on every
test run. Code matches it, not a paper. Scalar mass is `mu`; the Myers-Perry
mass parameter is `M`. `Phi ~ e^{-i omega t}`, so damped modes have
`Im omega < 0`.

## Symmetry facts to use, not rederive

* `c^2 = 4 s delta (omega^2 - mu^2)` exactly; the angular operator is the round
  S^3 Laplacian on `delta = 0` **and** on `s = 0`.
* At `delta = 0`: `Lambda = l(l+2)`, and the radial problem depends on
  `(m1,m2)` only through `m = m1+m2`. Multiplet dimension is exactly `l+1`.
* Symmetry group is Klein four `{1, E, P, EP}` acting on
  (parameters x sector labels); see `docs/SYMMETRY_GROUP.md`. Diagonal sectors
  are even in `delta`; anti-diagonal sectors are even in `s`; sector `(0,0)` has
  both.
* Diagonal-sector exceptional lines must come in mirror pairs and can terminate
  quadratically on the equal-spin surface.

## Environment

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.lock
PYTHONPATH=src .venv/bin/python -m pytest -q
```

`uv` and Julia are not installed. `python-flint` (Arb) **is** available.

See `SOLVERS.md` and `STATUS_SCHEMA.md` in this directory.

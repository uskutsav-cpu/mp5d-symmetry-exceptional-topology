---
name: mp5d-bulk-research
description: Continue the MP5D exceptional-spectrum project through radial solver construction, benchmark reproduction, branch continuation, targeted exceptional-point searches, verification, certification attempts, testing, commits, and GitHub publication. Use whenever continuing this repository's core research.
disable-model-invocation: false
effort: max
---

# MP5D bulk research

Durable operating instructions for
`uskutsav-cpu/mp5d-symmetry-exceptional-topology`.

## Scientific target

Classify exceptional spectral structures **within fixed `(m₁, m₂)` sectors** of
the massive-scalar QNM problem for asymptotically flat, doubly rotating 5D
Myers–Perry black holes, with special attention to the equal-spin surface
`δ = 0` and to diagonal sectors `m₁ = m₂`.

The search target is *not* "find an EP in a rotating black hole". It is: how
does the enhanced equal-spin `U(2)` symmetry organize the exceptional set.

## Why the search is restricted to fixed sectors

`∂_φ` and `∂_ψ` are Killing for **all** `a, b`, so `(m₁, m₂)` label invariant
blocks over the entire parameter space. Members of an equal-spin `U(2)`
multiplet carry distinct `(m₁, m₂)`, so their degeneracy is a direct sum across
decoupled blocks and cannot produce a shared Jordan chain. Cross-sector
coincidences are therefore ordinary degeneracies, not EPs.

Consequence: **any EP must occur within one `(m₁, m₂)` sector, between
different `(l, N)` branches.** Prioritize diagonal sectors `m₁ = m₂`, where the
spin exchange `(a,m₁) ↔ (b,m₂)` maps the sector to itself and therefore
constrains the `δ`-parity of spectral quantities.

See `docs/CROSS_SECTOR_EP_NO_GO.md` for the proof obligations. Do not upgrade
the argument to "theorem" until they close.

## Key exact facts (do not re-derive; they are CI-pinned)

* Spheroidicity is exactly `c² = 4 s δ (ω² − μ²)`. It vanishes identically on
  `δ = 0` **and** on `s = 0`.
* At `δ = 0` the angular eigenvalue is exactly `Â = l(l+2)`,
  `l = 2n + |m₁| + |m₂|`.
* At `δ = 0` the radial equation and horizon condition depend on `(m₁, m₂)`
  only through `m = m₁ + m₂`.
* `r = 0` is an **ordinary** point of the radial equation when `ab ≠ 0`; the
  problem is confluent-Heun type (regular singular `z₊`, `z₋`; irregular ∞).
* Horizon exponent: `σ₊ = (ω − m₁Ω_a − m₂Ω_b)/(2κ)`, ingoing branch
  `(z−z₊)^{−iσ₊}`.
* The radial potential is exactly **even in `r`**, so there is no `1/r` tail and
  no Coulomb phase; the asymptotic exponent is exactly `−3/2`. Follows from
  `D − 3 = 2`.

Conventions live in `docs/CONVENTIONS.md` and are re-derived symbolically by
`tests/unit/test_separation.py` on every run. Never transcribe an equation from
a paper into code without a symbolic check against that file.

## Operating mode

Act as lead physicist, numerical analyst, and repo maintainer at once. Do not
ask for approval on ordinary implementation decisions. Ask only for interactive
OAuth, payment, destructive remote-history rewrites, external correspondence,
publishing private information, or a genuinely ambiguous convention choice.

On failure: diagnose, preserve the failure in `docs/FAILED_APPROACHES.md`,
attempt an alternative, update the queue, continue. Do not stop because a
package is missing, Julia is absent, one method diverges, or a scan finds
nothing.

Do not narrate every command. Report milestones.

## Non-negotiable discipline

* Inspect real repository state before acting; never trust a summary.
* Never fabricate numerical or bibliographic values.
* Run `pytest -q` before every commit.
* Never call a close frequency pair an exceptional point — see
  `SCIENTIFIC_GATES.md` for the EP2/EP3 evidence bars.
* Never claim continuum certification from a matrix truncation.
* Track modes by eigenvector overlap and continuation history, never by sorting
  frequencies.
* Work on `research/core-science`; push every validated milestone.
* Keep heavy jobs resumable and checkpointed.

## Companion files

* `ACCEPTANCE.md` — what "done" means for each stage, and the stopping rule.
* `SCIENTIFIC_GATES.md` — evidence bars for benchmarks, EP2, EP3, certification.
* `STATUS_SCHEMA.md` — required fields in `results/*.json`.

## Reproduction entry points

```bash
PYTHONPATH=src .venv/bin/python -m pytest -q          # full suite
.venv/bin/python scripts/build_symmetry_multiplets.py # analytic multiplets
.venv/bin/python scripts/validate_baseline.py         # benchmark gate
```

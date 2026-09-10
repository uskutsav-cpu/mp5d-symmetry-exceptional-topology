# mp5d-symmetry-exceptional-topology

Research code and results for the question: does the enhanced `U(2)` symmetry of
the equal-spin surface organize the exceptional-set topology of the
massive-scalar quasinormal-mode spectrum of the asymptotically flat, doubly
rotating five-dimensional Myers–Perry black hole?

Parameter space: `(s, δ, μ)` with `s = (a+b)/2`, `δ = (a−b)/2`, and `μ` the
scalar-field mass. `δ = 0` is the equal-spin surface. The Myers–Perry mass
parameter is `M`; it is never written `μ`.

`results/status.json` is the authoritative machine-readable status. Where this
README and that file disagree, the JSON is correct.

## Result

**Outcome C — a bounded negative classification, plus two exact structural
results.**

**Exact (independent of any numerics).** Because `∂_φ` and `∂_ψ` are Killing for
*arbitrary* rotation parameters, the quasinormal pencil is a direct sum over
azimuthal sectors `(m₁,m₂)`, and the Jordan chains of a direct sum are the union
of the summands' chains. A coincidence of frequencies from *different* sectors
is therefore always semisimple. The equal-spin `U(2)` multiplet degeneracies are
protected from splitting **and** forbidden from being defective by one and the
same fact — the most natural route to exceptional structure is closed exactly.
The surviving discrete symmetry is a Klein four-group acting on parameters *and*
sector labels; diagonal sectors `m₁=m₂` are exactly even in `δ`, which forces
any diagonal-sector exceptional line to occur in mirror pairs and terminate
quadratically on the equal-spin surface.

**Numerical.** Searching the one channel the no-go leaves open — collisions
between different `(ℓ,N)` branches inside a fixed sector — over 50,558 atlas
points in 7 sectors:

| quantity | bound over the searched domain |
| --- | --- |
| pairwise branch gap | **0.5116** |
| scale-free root separation `\|a₁/a₂\|` | **0.3107** |
| discriminant `\|D\| = \|ω₊−ω₋\|²` | **0.2618** |
| eigenvalue condition number at the strongest interaction | **8.73e3** (finite; ordinary control 18.3) |

An EP2 needs the first three to vanish together and the last to diverge. The ten
tightest candidates are **rejected** as avoided crossings by argument-principle,
monodromy and Puiseux tests. Domain: `s ≤ 0.42`, `|δ| ≤ 0.20`, `μ ≤ 1.8`,
`r₋ ≤ 0.2287`, extremality `≥ 0.294`; Solver A and the recurrence-free Solver C
agree to `2.06e-9` across it.

Details and the full limitation list: `docs/FINAL_NEGATIVE_RESULT.md`.

## What is *not* claimed

* **No exceptional point has been verified.** A close pair of frequencies is not
  an exceptional point, and this repository will not call one that.
* The numerical exclusion is **bounded**, over a **declared domain**, and is
  **not a theorem**. Its minimum sits on the boundary of the searched box. It
  must not be conflated with the cross-sector no-go, which *is* exact.
* **The quasiresonant regime is not searched.** There the recurrence-free
  solvers are inapplicable *by formulation* (the required complex-scaling angle
  → 90° as `Re Ω → 0`), so the evidentiary standard cannot be met. This is the
  one region of the studied sectors where an EP could still hide.
* **Novelty is not asserted.** Exceptional lines *are* established for the de
  Sitter members of this family; the `Λ=0` result here is a contrast between
  spacetimes, not a contradiction. See `docs/NOVELTY_AUDIT.md`.
* Certification is labelled by what it certifies. The angular eigenvalue of the
  truncated Jacobi matrix is enclosed in ball arithmetic; **the radial problem
  is not certified**, because the coefficients come from a DFT with no rigorous
  aliasing bound. No continuum certificate is claimed anywhere.
* Bibliography entries in `bibliography/` are **indexed**; only those actually
  relied upon have been metadata-verified, and `manuscript/references.bib`
  contains only those.

### Two defects worth knowing about

Both produced small computed gaps that looked physical, and both are preserved
in `docs/FAILED_APPROACHES.md` with their mechanism:

1. a **branch-label collapse** that reported a `4.3e-15` "gap" at 200 points —
   caught only because the search runs two diagnostics with different failure
   modes;
2. a **branch jump** admitted by a `|ω|`-relative tolerance, which contaminated
   the then-reported tightest candidate.

A prior bound based on `min |dF/dω|` has been **withdrawn**: that quantity is
not invariant under rescaling of the spectral condition and bounds nothing.

## Reading order

1. `docs/FINAL_NEGATIVE_RESULT.md` — the primary result and its limits.
2. `docs/CONVENTIONS.md` — metric, separation, signs. Pinned by
   `tests/unit/test_separation.py`.
3. `docs/CROSS_SECTOR_EP_NO_GO.md` — the no-go; all obligations discharged.
4. `docs/SYMMETRY_GROUP.md` — the Klein four-group and the termination law.
5. `docs/EP_VERIFICATION_GATE.md` — what a candidate must survive.
6. `docs/CLAIM_LEDGER.md` — every claim with its evidence level.
7. `docs/FAILED_APPROACHES.md` — preserved negative results.

## Install and run

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
PYTHONPATH=src .venv/bin/python -m pytest -q     # 252 tests
```

Reproduce the science:

```bash
./scripts/reproduce_core_science.sh --quick   # analytical + machinery, ~1 min
./scripts/reproduce_core_science.sh           # + atlas + collision search, ~1 h
./scripts/reproduce_full_science.sh --resume  # everything
```

The 24 MB branch atlas is **not** committed; it is regenerated by
`scripts/build_branch_atlas.py` and verified against the per-shard SHA-256 in
`results/full_branch_atlas.json`.

`uv` and Julia are not installed on the development machine; `requirements.lock`
is a `pip freeze` and plays the role of `uv.lock`. `python-flint` (Arb) is
available and used for ball arithmetic.

## Manuscript

`manuscript/` holds a complete draft (paper, supplement, two cover letters,
referee response template) and five figures generated from `results/*.json` —
no figure number is typed by hand. Build with:

```bash
./scripts/build_manuscript.sh
```

Requires `tectonic`. Figures additionally require `requirements-figures.lock`,
kept deliberately out of `requirements.lock` so that neither the results nor CI
depend on a plotting stack; the committed figure PDFs let the documents build
without it.

**Not submitted.** Authorship, affiliation and acknowledgments are intentionally
unassigned. PRL was considered and explicitly declined — the reasoning is
recorded in `manuscript/cover_letter_prl.tex` so the decision is deliberate
rather than an omission.

### Stability

APIs, file layouts, and JSON schemas change without notice. This is a research
repository, not a library.

## License

MIT, first-party code only. No third-party source is vendored. The bibliography
contains metadata only — no copyrighted article text or PDFs.

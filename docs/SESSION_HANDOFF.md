# Session handoff

## Session 8 (2026-08-02) — Outcome C reached: bounded multi-sector EP exclusion

### State

* Branch `research/absolute-final`, PR **#6** open against `main`, CI green.
* **252 tests passing.** Repository public, `main` protected.
* Atlas data (~24 MB) is **untracked**; regenerate with
  `scripts/build_branch_atlas.py` and verify against the sha256 per shard in
  `results/full_branch_atlas.json`.

### Achieved

**Analytical (exact, independent of the numerics)**

* Cross-sector no-go **closed unconditionally** (O6.5–O6.7 discharged). Jordan
  chains of a direct sum are the union of the summands' chains, so cross-sector
  coincidence cannot create defectiveness. O6.5 needs no `L²` — a strongly
  continuous torus action on any Banach space suffices. O6.6 was removed: the
  proof never used angular separation.
* Symmetry group is Klein four `{1,E,P,EP}` on (parameters × sector labels).
  Diagonal sectors even in `δ`, anti-diagonal even in `s`, `(0,0)` both.
  **C13 corrected** — false in general, true only in `(0,0)`.
* O7 discharged for all `l`: multiplet dimension `l+1`, level total `(l+1)²`.
* Termination law recorded **before** the search (C43, conditional).

**Numerical**

* Atlas: **50 558 points**, 7 sectors × 3 `l` × 4 `N`.
* Bounds: branch gap `≥ 0.5116`, scale-free separation `≥ 0.3107`,
  `min|D| ≥ 0.2618`, `κ ≤ 8.73e3` (ordinary control 18.3).
* Ten tightest candidates **rejected** — avoided crossings (G2/G3 pass,
  G1/G4/G6 fail).
* Long-lived branches tracked to damping `8.65e-10`.
* A-vs-C agreement `2.06e-9` up to `r₋ = 0.2287`.

**Corrections to the prior record**

* **C32 withdrawn**: `|dF/dω|` is not rescaling-invariant and bounds nothing.
* **Near-extremal breakdown is not extremality** — it is the quasiresonant
  limit (`Re Ω → 0` ⇒ required scaling angle → 90°).
* `r₂` is the inner horizon *radius*, not its square.

### Not done

* **Gate 16 certification is PARTIAL.** Krawczyk and Arb machinery exist and are
  validated both directions, and the angular eigenvalue is enclosed for the
  `N`-truncated Jacobi matrix. The **radial** problem is *not* certified: the
  polynomial coefficients come from a DFT with no rigorous aliasing bound, so
  `E_truncation` for the recurrence cannot currently be bounded. This is the
  exact obstruction, and it is why no "truncated-recurrence certified" label is
  used anywhere.
* **Solvers E (Wronskian) and F (hyperboloidal) not built** — deliberately
  deferred until gate 13 said whether a better solver could help. It now has:
  hyperboloidal is the principled fix for the quasiresonant regime.
* Gate 18 novelty audit done; gate 19 manuscript is a filled draft with
  authorship intentionally unassigned.
* `μ > 1.8`, `s > 0.45`, `N > 3` and higher `l` unsearched.

### Next, in order

1. **Solver F (hyperboloidal).** The only method that remains applicable where
   `Re Ω → 0`, i.e. the one region of the spectrum this search could not cover,
   and the only place an EP could still hide within the studied sectors.
2. Extend the atlas to `μ > 1.8` and `N > 3` once F exists.
3. Close gate 16 for the radial problem: either bound the DFT aliasing error
   rigorously, or replace the coefficient extraction with an exact symbolic
   route so Arb can enclose it end to end.

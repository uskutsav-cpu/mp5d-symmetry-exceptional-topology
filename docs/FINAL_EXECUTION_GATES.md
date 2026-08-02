# Final execution gates

A gate is COMPLETE only when its acceptance condition passes — never because
code exists. Status vocabulary: COMPLETE / PARTIAL / NOT STARTED / BLOCKED.

| # | Gate | Status | Evidence |
| --- | --- | --- | --- |
| 1 | Repository recovery and audit | **COMPLETE** | `docs/FINAL_RECOVERY_AUDIT.md`; PR #5 merged after CI-green check |
| 2 | Public repository, CI, protected `main` | **COMPLETE** | public; CI green; branch protection active |
| 3 | Persistent project configuration | **COMPLETE** | `.claude/skills/mp5d-complete-project/` |
| 4 | Solvers A–D with declared domains | **COMPLETE** | `docs/FINAL_SOLVER_DOMAINS.md` |
| 5 | Published-benchmark gate | **COMPLETE** (prior sessions) | three-solver agreement `<1e-5` on 10/10 |
| 6 | Cross-sector no-go, all obligations | **COMPLETE** | `docs/CROSS_SECTOR_EP_NO_GO.md`; O6.5–O6.7 discharged |
| 7 | Symmetry group and `δ`-parity | **COMPLETE** | `docs/SYMMETRY_GROUP.md`; C13 corrected; O7 discharged for all `l` |
| 8 | EP machinery validated on synthetic systems | **COMPLETE** | 21 tests incl. negative controls |
| 9 | Fixed-sector branch atlas | see `results/status.json` | 7 sectors × 3 `l` × 4 `N` |
| 10 | Adaptive collision search | see `results/status.json` | two-tier gap + scale-free separation |
| 11 | EP verification / rejection | see `results/status.json` | gate in `docs/EP_VERIFICATION_GATE.md` |
| 12 | Long-lived / quasi-bound branches | see `results/status.json` | adaptive-step `μ` continuation |
| 13 | Near-extremal spectral classification | see `results/status.json` | `θ`-independence test |
| 14 | Effective local model | see `results/status.json` | `D = (ω₊−ω₋)²`, out-of-sample scored |
| 15 | Targeted pseudospectra | see `results/status.json` | operator, norm and caveat stated |
| 16 | Certification | see `results/status.json` | scope-labelled per `STATUS_SCHEMA.md` |
| 17 | Reproducibility freeze | see `results/status.json` | `scripts/reproduce_*.sh` |
| 18 | Novelty audit | see `results/status.json` | — |
| 19 | Manuscript package | see `results/status.json` | — |

Gates 9–19 are updated from the machine-readable state rather than by hand, so
this table cannot drift from `results/status.json`.

## Deliberately out of scope, with reasons

* **Solver E (complex-contour Wronskian) and Solver F (hyperboloidal).** Two
  recurrence-free solvers (C and D) already exist and agree with each other in
  the ordinary regime. A third and fourth would add confirmation only in the
  near-extremal region, and the `θ`-independence test (gate 13) settles the
  *nature* of that region far more cheaply than building two more solvers. If
  gate 13 shows an isolated resonance that the existing solvers merely fail to
  resolve, Solver E becomes the right next step; if it shows continuum or
  accumulation, no solver would have helped and building them would have been
  wasted effort. The decision is therefore made **after** the diagnostic, not
  before it.
* **Higher raw-recurrence tail orders (`u₂`, `u₃`).** Blocked by the Gaussian
  reduction degradation past `n ≈ 500–800`, recorded in
  `docs/FAILED_APPROACHES.md`. The leading order already gives ~1.35× depth
  reduction and the benchmark agreements do not depend on it.

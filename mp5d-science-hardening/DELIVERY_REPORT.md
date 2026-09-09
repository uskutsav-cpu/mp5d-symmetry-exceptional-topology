# MP5D hardening — delivery and executed validation

**Delivery: implemented research-hardening code and development evidence.**
**Scientific result: `SCIENCE_FREEZE_BLOCKED`. All twelve research tasks are not complete.**

The baseline is `research/absolute-final` at
`72abf0bf68d41d07de0fbd5056c40b7a713a6cb5`, not the older default `main`.
GitHub rejected the write attempt. No remote branch, pull request, commit,
merge, workflow run or freeze tag was created.

## Executed checks

| Check | Actual outcome |
|---|---|
| New numerical/validation tests, including their slow-marked tests | **126 passed; 0 failed; 0 skipped** |
| Installer safety tests on temporary synthetic Git fixtures | **9 passed; 0 failed; 0 skipped** |
| Ten closest interactions | **120 recorded solver evaluations** |
| Saved N2 frequencies: three CF depths each | **30/30 failed cross-inversion validity** |
| Saved N2 frequencies: three corrected C resolutions each | **30/30 failed same-root independent confirmation** |
| Saved N3 frequencies | **60/60 passed the audit's pointwise checks** |
| Published static benchmarks | **90 evaluations; 8/15 frequencies pass all six settings; 7/15 unresolved** |
| Strongest interaction, seven-axis ladder attempt | **23 rungs fail physical-pair preflight; all seven axes inconclusive** |
| Adaptive physical boundary search | **Blocked at invalid initial pair; no minimum produced** |
| Exact transformed radial degree/horizon identities | **Passed; degrees (9, 8, 7)** |
| Exact-rational finite Krawczyk certificate | **Unique finite-polynomial root certified, not continuum QNM** |
| Hyperboloidal Pöschl–Teller testbed | **Testbed check passes; not a massive MP5D Solver F** |
| Lean kernel build | **Not run: Lake unavailable** |
| Result-envelope integrity/schema | **Passed without promoting scientific status** |
| Authoritative provenance / science freeze | **Correctly rejected / blocked** |
| Full upstream tests, upstream Ruff lint and GitHub CI | **Not executed** |

Full logs, XML results and JSON receipts are in `evidence/`. CSV files expose
actual recorded solver samples without producing publication figures.
`RUN_CHECKS.py` reproduces the standalone checks and returns a nonzero overall
exit while science is blocked. Interrupted candidate computations can be
resumed with `--resume-candidates` only when source, environment, settings,
fixture, order and integrity match; no result is promoted or relabeled.

## Important new finding

At the tightest saved location, the historical N2 frequency is approximately

```
omega = 2.5094126155849015 - 1.5480379310838372 i
```

The new arbitrary-precision continued-fraction calculation reproduces this
frequency to about fourteen decimal places, with its chosen inversion residual
around `1e-43` at the three tested depths. Nevertheless, other inversions do not
agree. Corrected recurrence-free Solver C converges to a *different* frequency,
near the saved N3 root, about `0.51162` away. A tiny selected residual plus depth
stability therefore cannot validate this historical branch on their own.

This pattern occurs for all ten saved N2 roots in the audit. It is evidence of
an inversion-specific numerical/root-identity problem, not a rigorous proof
that those continuum eigenfrequencies are impossible. The affected branches
must be independently reidentified. A substitute root must not be called N2
merely because it has a plausible damping rate or is nearby.

The old minimum gap `0.5116202072812408` is **not recertified**. This audit does
not find an EP, establish a new exclusion bound or validate the full old atlas.
The N3 checks use a `1e-5` same-root agreement tolerance; some Solver C errors
worsen as radial resolution increases. Their pointwise pass is **not** a claim
of arbitrarily high accuracy or a completed convergence ladder.

One implementation error underlying the new C solver work is exact algebra:
for the smallest left/right singular vectors of the same matrix M,
`u* M v = sigma_min`. This is a nonnegative singular value, not the previously
claimed holomorphic complex root function. The replacement uses fixed borders
and fixed row scaling for a locally analytic residual and separately tests the
radiation sheet and differential-equation residual.

## Checklist-by-checklist status

| Requested task | Implemented / executed | Still unresolved |
|---|---|---|
| 1. Correct EP evidence | Guarded EP2 gate, Jordan/nullity checks, Puiseux/monodromy controls, analytic-contour interfaces; unsafe legacy APIs retired by installer | No physical MP5D EP2 or continuum exclusion certified; legacy consumers must be migrated |
| 2. Ten permanent regressions | Historical-frequency fixtures, adversarial tracker tests, 120 physical audit evaluations | Historic N2 validity/identity fails; old rejection labels cannot be assumed correct |
| 3. Strong ladders | Seven-axis engine and physical adapters; all ten pairs tested at three A depths and three C resolutions; 23-rung strongest-point preflight attempted | No seven-axis convergence claim survives for the invalid pair |
| 4. Adaptive boundary | Tested cached, bounded, multi-route refinement engine and physical adapter | Actual search blocked before a valid minimum can be produced |
| 5. Quasiresonant Solver F | Explicit exclusion from proposed bounded scope; working exactly solvable hyperboloidal testbed | Full massive rotating MP5D hyperboloidal solver is not implemented or validated |
| 6. Higher modes | Mode-label generator and explicit N4/N5/admissible-ell probe configuration | Complete physical higher-mode probe and associated producer unfinished |
| 7. Radial certification | Direct symbolic coefficients; degree proof; optional Arb adapter; exact finite determinant and rigorous rational interval certificate | Arb adapter unexecuted; generic rotating finite certification and continuum tail bounds not established |
| 8. Submission CI | Fail-fast script and exact-commit manual workflow; no slow-test exclusion or swallowed scientific failures | Full upstream integration, original linter and Actions unrun; missing atlas/robustness/bibliography producers not magically supplied |
| 9. Provenance/claims | Hashed receipts; exact-commit, clean-source, environment/schema and claim checks; old status demoted without deleting history | Current receipts are development-only, not authoritative committed runs |
| 10. Literature inputs | Fifteen static values rechecked against original tables, with conversion and labels in a manifest | Full rotating/massive benchmark and bibliography audit incomplete; seven static cases need stronger convergence |
| 11. Lean | Pinned exact-layer proof drafts, no admissions inserted, failure-aware build receipt | Kernel compilation absent; full multiplet bijection and any continuum statement unproved here |
| 12. Freeze | Gate rejects current state; no publication figures generated from unfrozen data | **No science-freeze commit or tag exists** |

## Finite certificate: precise scope

The executed exact-rational interval calculation certifies a unique root of the
static massless `ell=0`, `K=12` raw finite radial determinant (degree 24) in a
rectangular complex box of radius `1e-30`. Its center is approximately

```
0.5344254629187461160842537200986021792118152145656172887
- 0.3832014857234140365562691232757638539024787873709776074 i
```

The exact coefficients, rational endpoints, strict Krawczyk inclusion and
contraction bound are in `evidence/benchmarks.json`. The certificate is not for
the infinite boundary-value problem, not generic rotating MP5D, and not the
continuum fundamental frequency. The implementation's exact-rational arithmetic
was executed and tested, but the interval program is not formally verified in
Lean.

## Source and execution provenance

The numerical environment was CPython 3.13.5, NumPy 2.3.5, SciPy 1.17.0,
SymPy 1.14.0, mpmath 1.3.0 and pytest 9.0.2 on Linux x86-64. python-flint and
Lean/Lake were unavailable. The delivered locks distinguish the actually
executed numerical environment from the unexecuted optional Arb dependency.
macOS and the full upstream installation have not been tested.

Every scientific receipt is explicitly `DEVELOPMENT_ONLY`, with no invented
source Git commit. The standalone overlay's source-file hash set is retained.
A candidate run was interrupted by the execution-tool timeout and safely
resumed from four completed records under the identical verified source and
environment. All ten records are now present. A complete standalone rerun does
not require access to the earlier interrupted run.

Published static values were checked visually against Matyjasek (2021),
arXiv:2107.04815, Tables I–III on printed pages 12–13. The CF/HD column uses
`omega_tilde=omega/T_H`, with `T_H=1/(2*pi)`, so the conversion is division by
`2*pi`. The manifest preserves table, column, mode labels and normalization.
No figures were reconstructed from screenshots and no OCR was used.

Macedo and Zenginoglu (2025), arXiv:2409.11478v2, sections 4.2–4.3, informs the
separation between a hyperboloidal coordinate testbed and a fully validated
radiation/regularity formulation. It does not supply a completed massive MP5D
Solver F for this project.

The Lean drafts pin mathlib commit
`e37d88a26f3791ed5a93daa1f949af1021b8d103` and its declared toolchain
`leanprover/lean4:v4.34.0-rc2`, retrieved from the official repository. The
original methods and defect locations are identified in the source comments
and `overlay/docs/SCIENCE_HARDENING.md`.

## Applying the bundle

Follow `README.md`: use a fresh review checkout of the exact audited base,
create a new branch, run `APPLY_UPDATE.py --dry-run`, inspect the proposed
changes and only then apply. The bundle is not a complete upstream clone and
should not be uploaded to GitHub as an opaque ZIP. Existing unsafe APIs are
intentionally disabled and may make old scripts/tests fail visibly. That is
an integration obligation, not permission to ignore failed gates.

**Do not merge to main or label this as submission-ready on the strength of
passing new software tests. Repair the physical evidence and complete the
missing scientific/integration work first.**

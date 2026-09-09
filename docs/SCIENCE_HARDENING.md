# Scientific hardening and revalidation

## Status: not a science freeze

This update targets `research/absolute-final` at
`72abf0bf68d41d07de0fbd5056c40b7a713a6cb5`. It is a research-hardening
implementation, not a claim that the original twelve-item checklist has been
fully discharged. The original `main` branch is not the audited baseline.

The ten closest historical interactions require revalidation. The new audit
performed 120 evaluations: ten interactions, two saved frequencies, three
continued-fraction depths and three independent collocation resolutions.
All ten saved `N2` frequencies failed the cross-inversion check at all three
depths. Corrected Solver C did not independently confirm those same frequencies
at any of the three tested resolutions. The saved `N3` frequencies passed the
checks in this run. These results challenge the identity/validity of the
historical roots; they do not establish an EP, a replacement overtone label,
or a continuum exclusion theorem.

The historically reported minimum gap `0.5116202072812408` is therefore not
recertified by this update. The installer wraps the previous `results/status.json`
as explicitly historical data and retains its original provenance without
inventing a generating commit. The new claim registry is
`config/science/claims.json`, with `science_ready: false`.

The numerical receipts distributed with the bundle have `DEVELOPMENT_ONLY`
authority. They were generated from the standalone overlay, not a clean full
checkout of the upstream repository. Their exact source-file and environment
hashes are recorded. They must not be renamed into an authoritative freeze.

## What changed

### Branch tracking

`mp5d_science.tracking` implements global one-to-one assignment, comparison with
the second-best assignment, phase-invariant eigenvector overlap, predictor
checks, uncertainty-aware duplicate rejection and adaptive step subdivision.
Accepted roots are paired with their actual accepted parameter values. Rejected
steps cannot advance the continuation history. Budget exhaustion and persistent
ambiguity produce explicit incomplete results.

The installer corrects the original predictor call from
`extrapolate(hist, ts[:-1] + [arc], arc)` to
`extrapolate(hist, ts[:-1], arc)`. The new tracker is a separate explicit API;
the original single-branch routine is not misrepresented as having acquired
all multi-branch safeguards automatically.

Synthetic regression tests exercise forward and reverse paths, independent
seeds, step subdivision, avoided crossings, eigenvector-assisted identity,
greedy-assignment failures, duplicate labels, and one- and two-loop EP2
transport. EP3 and semisimple controls are included. The ten physical fixtures
store the historical branch labels and actual saved frequencies, not unreliable
roots extracted by the old meromorphic contour moments. Their intended result
is revalidation, not unconditional acceptance of an old rejection verdict.

### EP evidence and contours

`evidence.py` separates a resolved pointwise gap, numerical support for EP2 and
an inconclusive calculation. Numerical failure is not evidence that an EP is
absent. An EP2 gate requires a combination of independent root location,
nullity-one/Jordan evidence, square-root splitting, monodromy and independent
radial solver lineages. A/B computations sharing the recurrence are not counted
as two independent radial formulations.

The withdrawn Taylor ratio and the discriminant-as-independent-evidence cannot
be used by the new gates. The legacy `root_separation`, unrestricted
`zero_count`, contour-moment and Jordan entry points are explicitly retired
rather than silently accepting their former assumptions. Calls raise an
informative error. The old EP synthetic test file is replaced with retirement
tests; the stronger positive and negative controls reside under
`tests/submission/`.

`spectral.py` requires an explicitly scoped analytic spectral-function object.
Its constructors distinguish polynomials and finite analytic determinants from
uncontrolled meromorphic continued fractions. Contours undergo refinement,
phase-increment, finite-value, dynamic-range and boundary-zero checks. Moment
reconstruction uses centered/scaled coordinates. These are numerical contour
checks, not interval contour certificates. Declaring a general callable
analytic still requires the caller to justify that assertion over its domain.

### Corrected Solver A and Solver C

The new Solver A uses algebraically assembled polynomial coefficients and true
mpmath working precision. A small residual in one inverted continued fraction
is not sufficient: multiple inversions are checked at the final frequency. The
legacy Solver A wrapper retains its returned frequency but refuses a successful
`converged` flag when this additional validity check fails.

The old Solver C scalar `u.conj() @ M @ v`, with `u` and `v` chosen from the SVD
of that same matrix at every iterate, equals the smallest singular value. It is
not the claimed holomorphic complex root function. The replacement uses a
locally analytic bordered-system residual with fixed border vectors and fixed
row scaling. Radiation-sheet, finite-contour suppression, singular-value and
off-grid differential-equation residual checks remain explicit. The original
Solver C public function is redirected to this implementation by the installer.

A corrected solver can converge to a *different* root. Candidate confirmation
therefore also requires proximity to the historical root and branch identity;
a small matrix residual alone is never reported as independent confirmation.
The two radial formulations remain coupled to the same angular model, which
limits their independence and is recorded as such.

### Exact radial algebra and finite certification

`radial_polynomial.py` constructs the transformed coefficients algebraically,
without FFT coefficient recovery or numerical degree trimming. The symbolic
check establishes polynomial degree bounds `(9, 8, 7)` in the reduced
representation, using the horizon and dispersion identities. The calculation
is an exact algebraic identity, not an infinite-recurrence convergence proof.

A common `(1-u)^4` factor is restored for the Gaussian-reduced continued-fraction
representation to match the recurrence setting being audited. Removing common
polynomial factors preserves an interior differential equation, but does not
by itself justify a new minimal-solution selection for an infinite recurrence.
The implementation does not conflate those two facts.

An optional `ball_coefficients` adapter uses python-flint. That adapter has not
been executed in the delivery environment, where python-flint was absent.
Do not describe the generic rotating coefficient path as Arb-certified.

`certification.py` also provides an independently executable exact-rational
interval Krawczyk calculation. Binary floats are rejected as exact inputs. A
strict inclusion and contraction check certifies a unique zero of an explicitly
specified finite polynomial. The executed example is the raw finite radial
Hill determinant for the static massless `ell=0`, `K=12` problem: degree 24,
with a rectangular root box of radius `1e-30`. Its root is approximately
`0.5344254629187461 - 0.3832014857234140 i`.

This is deliberately a certificate for a **finite determinant**, not the
continuum frequency and not generic rotating MP5D certification. The exact
polynomial, rational bounds, preconditioner-derived inclusion and contraction
bound are retained in the benchmark receipt. Python/SymPy exact arithmetic has
been tested; this is not a Lean kernel certificate of the interval program.

### Convergence and adaptive refinement

`convergence.py` verifies that actual solver settings change between ladder
rungs and that reported effective precision/resolution match the request.
Failed or missing rungs cannot be silently skipped. `research.py` connects the
physical solvers to depth, angular, precision, radial, contour-length,
complex-angle and continuation-step ladders. Continuation step changes trigger
actual path transport rather than relabeling a fixed-frequency solve.

The generic adaptive boundary engine is implemented and tested. It caches
bounded local stencils, halves their scale, tracks unsuccessful evaluations and
reports only local numerical stabilization, never a global lower bound.
The physical wrapper approaches points through distinct continuation routes
and checks independent Solver C agreement.

The default physical boundary/ladders preflight currently stops at the invalid
historical anchor. This is the correct blocked state, not a completed adaptive
minimum search. The executed ten-candidate audit covers three depth and three
radial rungs, not all seven axes. The configured higher-mode labels `N=4,5`
and subsequent admissible `ell` branches are a plan; the corresponding complete
physical probe and refreshed atlas have not been executed or packaged as a
finished end-to-end producer.

### Published benchmarks

`data/regressions/matyjasek_2021.json` contains 15 numbers independently checked
against Tables I–III on printed pages 12–13 of Jerzy Matyjasek,
*Accurate quasinormal modes of the five-dimensional Schwarzschild–Tangherlini
black holes*, arXiv:2107.04815 (2021). Values come from the CF/HD column,
not the neighboring WKB column. The table's convention is
`omega_tilde = omega/T_H`, `T_H=1/(2*pi)` for the normalized horizon;
conversion in the runner is therefore `omega=omega_tilde/(2*pi)`.

The audit makes 90 solver evaluations: 15 frequencies with three A and three C
settings each. At the declared `1e-5` comparison threshold, eight frequencies
pass every tested setting and seven do not. The discrepancies and failures are
retained. This rechecks that static benchmark subset; it is not a completed
independent audit of every rotating or massive benchmark and bibliography item
in the original repository.

### Hyperboloidal and Lean scope

`hyperboloidal.py` is a real Pöschl–Teller quadratic-pencil testbed. With
`y=tanh(x)` and `tau=t-log(cosh(x))`, its separated equation is

```
(1-y^2) u'' + 2(i omega - 1) y u' + (omega^2 + i omega - V0) u = 0.
```

It is tested against its known frequencies. It is **not** a validated massive
MP5D Solver F. The massive quasiresonant regime is explicitly outside the
current proposed submission scope. A full physical Solver F still needs the
correct radiation/regularity class, angular coupling, derivation and independent
benchmark validation. Macedo and Zenginoglu, *Hyperboloidal Approach to
Quasinormal Modes*, arXiv:2409.11478v2 (2025), especially sections 4.2–4.3,
distinguish geometric regularization from the analytic mode-selection problem.

`lean/` contains proof drafts for product/block linear-map structure, invariant
sectors, Klein-four actions, spectral-set parity and elementary identities.
Mathlib and its declared toolchain are pinned. No `sorry`, `admit` or new axiom
is inserted. **Compilation has not been checked:** Lean/Lake were unavailable.
The `Fin (ell+1)` cardinality result is not the complete bijection proof for all
allowed azimuthal pairs. No finite algebraic theorem is advertised as a
continuum-QNM theorem. See `lean/README.md` for exact limits.

### Provenance, installation and CI

All new scientific receipts record generating source commit when available,
source-file hashes, environment hash, solver/version, effective settings and
actual residuals. JSON serialization rejects non-finite values. Hashes detect
accidental alterations; they are not cryptographic signatures establishing who
ran a calculation. Authoritative runs require a clean committed checkout and
an externally specified exact source commit, not merely an ancestor.

Results are generated **after** committing source, in an external output
directory. This avoids the circular request to embed a Git commit's own hash
inside a file already contained in that same commit. Historical data is never
relabeled as generated by the new commit.

`reproduce_submission_science.sh` installs locked dependencies, runs the
original and new tests without excluding slow tests, runs the available audit
producers and requires every configured final artifact. Errors propagate. It
will currently fail, and that is intentional. Fresh atlas, complete robustness
and bibliography producers are not fully integrated, some required artifacts
remain absent, and numerical audits fail. The existence of a strict shell
wrapper is not completion of submission-science reproduction.

The original full test suite, original linter, full patched-repository
integration and GitHub Actions have **not** been executed in the delivery
environment. The new tests and bundle-installer tests were executed separately.
The installer tests use temporary synthetic Git repositories to test safety;
they are not a substitute for upstream integration testing. GitHub rejected the
write attempt, so no remote branch, PR, commit, merge, tag or workflow run was
created by this work.

## Running the delivered code

From the standalone bundle, using Python 3.13 and a suitable virtual environment:

```bash
python -m pip install -r overlay/requirements-science.lock
PYTHONPATH=overlay/src python -m pytest -c overlay/pytest-science.ini overlay/tests/submission -q
python -m pytest test_bundle_installer.py -q
python overlay/scripts/run_candidate_audit.py --development --dps 45 --output candidate_audit.json
python overlay/scripts/run_benchmark_audit.py --development --output benchmark_audit.json
python overlay/scripts/run_boundary_or_ladders.py boundary --development --output boundary.json
python overlay/scripts/check_lean.py --development --output lean.json
```

The last four commands can return exit code 2 with a saved failure/blocked
receipt. Do not convert those exit codes to successes. `--development` records
limited provenance; it does not change the mathematical acceptance criteria.

Use the top-level bundle README for the exact-base, new-branch installation
procedure. Do not upload this ZIP as a single opaque file into the repository.
Do not merge to `main`, publish new bound figures or create `science-freeze-v1`
on the strength of these development receipts.

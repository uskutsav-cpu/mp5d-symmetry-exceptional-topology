# Fresh MP5D numerical campaign — evidence schema 2

This adds an executable campaign around the existing `mp5d_science` physical
solvers. It does **not** report that the outstanding numerical study has passed.
It does not replace the physical MP5D equation with a test potential.

## Scope and entry point

Use `python -m mp5d_campaign.cli`, with `PYTHONPATH=src`. Version-2 receipts have
an explicit schema version and a new strict validator. Legacy version-1 scripts
are left unchanged; these files are not disguised as their version-1 envelopes.
The publication/release entry point for this campaign is the version-2 CLI.

The default atlas is a **declared finite grid in (m1,m2,ell)=(2,0,4)**, augmented
with the historical dangerous coordinates. It is not an enumeration of every
admissible angular sector or proof of coverage between samples. The sector list
and grid are explicit committed inputs in `config/campaign/campaign.json`.
Review this finite scope before making a manuscript claim about a larger atlas.
The formal overtones are 0–3. Robustness transports 0,4,5 **jointly** in each of
ell, ell+2, ell+4, so 4 and 5 cannot separately alias onto one accepted root.

Quasiresonances are excluded. The operational numerical scope also requires
`-Im(omega) >= 0.001` and `abs(sqrt(omega^2-mu^2)) >= 0.05`, with M=1. These
preregistered cutoffs are not a theorem identifying the onset of quasiresonance.
All physical contour-radiation checks must also pass. An excluded or failed point
is recorded, not silently removed from required task coverage. There is no
physical massive MP5D hyperboloidal Solver F in this update. Pöschl–Teller tests
are not a replacement. No continuum, global no-EP, or global-minimum theorem is
activated by this numerical campaign.

## Implemented producers

- `atlas.json`: published static anchors, three-resolution anchor preflight,
  continuous-separation-constant transport to integer ell, one-to-one adaptive
  joint continuation, and two different parameter-space approaches. Historical
  frequencies and branch names are never consumed. `VALIDATED` means all
  declared finite tasks passed, not global domain coverage.
- `candidates.json`: requires all ten historical coordinates to be accounted for
  using the fresh atlas. It does not reaffirm the old close-pair interpretation.
- `convergence.json`: all exact ladders from the merged research plan. The
  terminal-three *pairwise diameter*, not just neighboring-rung differences,
  measures empirical drift. Both A and C root drift are inspected. Precision
  sweeps genuinely vary A; C is explicitly still float64. Continuation-step
  sweeps traverse nontrivial local closed routes from each validated point, with
  different terminal approach axes. This tests local transport-step stability,
  not a fresh complete atlas reconstruction at each step. A repeated evaluation at a fixed
  point cannot masquerade as a continuation-step study. Every rung must also
  remain inside its fixed-label identity neighborhood; agreement between A and C
  alone cannot silently exchange two mode labels.
- `boundary.json`: starts at the validated fresh finite-atlas minimum. Complete
  bounded 3-dimensional stencils are evaluated from independent approaches.
  The radius does not shrink while significant descent continues. Two complete
  stable terminal polls and a new full seven-axis study at the refined candidate
  are required. The accepted state is local numerical stencil stability only.
- `robustness.json`: joint high-mode identity, two approaches, A/C agreement, and
  independent CF-depth and radial-resolution ladders at a separately registered
  high-overtone contour. These do not substitute for the core seven-axis study.
  A missing source-labeled N=5 anchor
  blocks acceptance. An inversion index alone cannot manufacture that label.
- `benchmarks.json`: all 15 configured static frequencies plus the six rotating/
  massive Table-II reference rows, each at three A depths and three C radial
  resolutions, with explicit M/rh/2pi conventions. The benchmark contour uses
  88.5 degrees and length 120, explicitly separate from the core seven-axis ladders. Foundation regression tests
  are also executed. This numerical subset is **not** a complete literature audit.
- `bibliography.json`: a separate complete-input inventory, metadata, conventions,
  transcription, applicability, primary-source evidence and review-note hashes.
  Every indexed bibliography row needs a review or a justified not-used
  disposition. Existing metadata flags and eight passing benchmarks do not make
  the remaining sources verified. The supplied review inventory is deliberately
  marked incomplete because that audit has not been supplied.
- `lean.json`: actually invokes Lake, records toolchain output, resolves pinned
  dependencies, runs `lake build`, and separately compiles every first-party Lean
  file to an `.olean`. Missing executables or unsuccessful compiler return codes
  never become `BUILD_PASS`. This only certifies compilation of those statements,
  not the desired future full operator/bijection theorem.
- `arb.json`: optional executed generic rotating ball-coefficient evaluation.
  Its scope is finite coefficient evaluation, not a rotating Krawczyk inclusion,
  not a continuum remainder estimate, and not continuum QNM certification.

## Exact seven-axis ladders

| Axis | Required values |
|---|---|
| CF depth | 160, 240, 320, 480 |
| Angular resolution | 24, 40, 64 |
| Decimal precision | 30, 50, 80 |
| Radial resolution | 160, 220, 280, 360 |
| Contour length | 40, 60, 90 |
| Scaling angle, degrees | 60, 65, 70 |
| Continuation step | 0.01, 0.005, 0.0025 |

Anchor preflight has separately recorded settings, including its steeper contour.
It is never counted as execution of these seven ladders. A rung that violates the
outgoing-contour criterion fails; it is not silently retried at a different angle
and reported as a pass at the requested angle.

## Installation and execution

Apply the additive package, review its source and scope, then commit it before
an authoritative run. No source registry readiness flag is changed by installation.
Use the repository-pinned Python dependency environment. A Lean/Lake installation
compatible with `lean/lean-toolchain` and the pinned mathlib revision is separately
required for actual proof compilation. The authoritative runner checks installed
versions against simple `name==version` pins in the committed requirements lock.
It also records the full installed-package inventory and machine/Python details.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD/src"
python -m mp5d_campaign.cli plan --root "$PWD"
python -m pytest -q tests/campaign

# Development-only diagnostic: cannot validate an atlas or release a paper.
python -m mp5d_campaign.cli diagnose --root "$PWD" \
  --output "$HOME/mp5d-evidence/inversion-diagnostic"

# Full post-commit run; choose a NEW directory outside the source checkout.
python -m mp5d_campaign.cli run --root "$PWD" \
  --output "$HOME/mp5d-evidence/campaign-run-001" --stages all
```

A full `run` returns exit code 2 unless every required scientific stage status
passes. It still preserves negative receipts and proceeds to independent audits,
such as bibliography and Lean, when numerical prerequisites fail. A failure is
not a software test success and not scientific completion. Per-task failure
records include the exception and a traceback where an exception occurred.

The staged commands are `benchmarks,atlas,candidates,convergence,boundary,
robustness,bibliography,lean,arb`. For example:

```bash
python -m mp5d_campaign.cli run --development --root "$PWD" \
  --output "$HOME/mp5d-evidence/build-and-audit-check" \
  --stages bibliography,lean,arb
```

## Resume and negative evidence

`--resume` requires identical source, configuration, environment, and authority.
Task cache keys include actual solver inputs, seeds, resolution, and run context.
Failed tasks are cached as failures; use `--retry-failed` to retry them under the
same context. Changed source or dependencies require a fresh committed run and
new evidence directory, not recycling the old run's accepted flags.

```bash
python -m mp5d_campaign.cli run --root "$PWD" \
  --output "$HOME/mp5d-evidence/campaign-run-001" --resume --stages all
```

The store has a single-writer lock. After an interrupted process, inspect the PID
in `.writer.lock` and confirm that no process is using the directory before
removing that lock. No automatic lock stealing or simultaneous mpmath threads
are used. `--point-limit` is allowed only with `--development`, remains explicitly
partial even when a selected point passes, and does not bound the anchor preflight.

## Missing scientific inputs must not be fabricated

The supplied static table contains N=0–4 for ell=0,1,2; **it contains no trusted
N=5 frequency**. Obtain and verify that label independently before adding an
anchor. A numerical guess obtained using inversion 5 is not equivalent to such
verification. New anchors require a source/normalization and a label audit;
additional numerical rows must also enter the submission bibliography inventory.

`rotating_references.json` covers Table II only. The remaining rotating/massive
literature and tables must be added to the committed reference inventory and
source reviews. `source_reviews.json` is an explicit incomplete audit record,
not a fake completed bibliography. The numerical producer can test configured
values, but cannot certify that unseen primary sources were read.

Existing finite exact Krawczyk results retain their original scope. Running the
optional Arb coefficient audit does not supply tail bounds, a continuum Fredholm
argument, or a generic rotating inclusion theorem. Those stronger claims remain
out of scope for this numerical-paper workflow.

## Claims and release without circular provenance

Readiness is not set by `run`. After all eight required artifacts pass and their
source/environment/task coverage is verified:

```bash
python -m mp5d_campaign.cli prepare-claims --root "$PWD" \
  --output "$HOME/mp5d-evidence/campaign-run-001" --resume
```

This writes **an external proposal**, `claims.proposed.json`; it does not edit
tracked source. Review it, copy the approved change to `config/science/claims.json`,
and commit it. Then run the complete campaign again into a new directory, using
that new clean source commit. Old receipts are invalid for the new commit.

Only the fresh successful second run can release:

```bash
python -m mp5d_campaign.cli release --root "$PWD" \
  --output "$HOME/mp5d-evidence/campaign-run-002" --resume
# Optional explicit local annotated tag (never pushed automatically):
# add --tag numerical-v0.1.0
```

The release step rechecks the actual physical evidence structure, exact required
coverage, seven-axis rungs, refined-minimum polls, higher-mode label coverage,
Lean return codes, source inventory, and provenance before producing publication
PNG/SVG figures, plotted data, frozen claims, and a release manifest. The tag
receipt is separate so creating a tag does not alter the release file it hashes.

The installed GitHub workflow always offers software tests. The expensive manual
science job requires a separately provisioned self-hosted `mp5d-science` runner
with Lean/Lake; adding this workflow does not create or pay for such a runner.
A software-green CI badge does not imply a science-green campaign.

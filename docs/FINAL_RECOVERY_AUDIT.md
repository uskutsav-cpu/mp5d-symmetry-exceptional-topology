# Final recovery audit

Verified state at the start of the absolute-final session, established by
running commands rather than by reading the previous handoff.

## Repository

| item | verified value |
| --- | --- |
| remote | `uskutsav-cpu/mp5d-symmetry-exceptional-topology` |
| visibility | **PUBLIC** (confirmed via `gh repo view`) |
| default branch | `main` |
| branch at session start | `research/final-science` @ `63e20ad`, clean tree |
| PRs | #1–#4 merged; **#5 open, CI green, `MERGEABLE`/`CLEAN`** |
| tags | `milestone/stage-A`, `milestone/radial-structure`, `milestone/first-qnm` |
| test suite at session start | **147 passed in 464 s** |

PR #5 satisfied the merge conditions (CI passing, no candidate described as
verified, Solver D limitations explicit, clean tree) and was merged. Work
continues on `research/absolute-final`.

## Environment

| item | value |
| --- | --- |
| Python | 3.13 venv at `.venv` (symlinked to miniconda base interpreter) |
| numpy / scipy / sympy | 2.5.1 / 1.18.0 / 1.14.0 |
| mpmath | 1.3.0 |
| python-flint | **0.9.0 present** (Arb ball arithmetic available) |
| Julia | **not installed** |
| `uv` | **not installed**; `requirements.lock` is a `pip freeze` |
| cores | 8 physical |

## Solvers, as verified in code and tests

| solver | formulation | recurrence-free | validated domain |
| --- | --- | --- | --- |
| A | Leaver continued fraction | no | `z₋ ≲ 0.20`, sub-extremal |
| B | raw recurrence / Hill / Wynn | no | as A |
| C | exterior complex-scaled collocation | **yes** (enforced by a grep test) | ordinary modes incl. `z₋ = 0.19` |
| D | multidomain near-horizon spectral | **yes** | ordinary modes only; fails at `z₋ = 0.44` |

`z₋` is the **inner** horizon radius squared, the near-extremality difficulty
parameter — not `r₊`. Three-solver agreement `< 1e-5` on 10/10 benchmark points.

## Claims found to need correction

**C32 was not a valid bound.** The recorded exclusion read "min `|dF/dω|` =
`1.02e-2` … therefore no EP2 in this box". The spectral condition `F` is defined
only up to an arbitrary nonvanishing analytic prefactor — Solver A's continued
fraction can be rescaled or inverted at a different depth without changing its
zero set — so `|dF/dω|` at a root carries no invariant meaning and can be made
arbitrarily small by rescaling. A bound built from it does not bound the
physics.

This is demonstrated, not asserted:
`tests/unit/test_ep_machinery_synthetic.py::test_root_separation_is_invariant_under_rescaling`
exhibits two spectral conditions with identical zero sets whose `|dF/dω|` differ
by more than `1e5`, while the scale-free diagnostic agrees to 2 %.

The exclusion has been **re-derived** using the invariant root separation
`|a₁/a₂|` and the pairwise branch gap. The conclusion (no EP2 in that box)
survives; the *reasoning* had to be replaced.

## Documentation found stale

* `CLAUDE.md` said "work on `research/core-science`" — superseded.
* `docs/EXECUTION_GATES.md` gate table predated Solver D and PR #5.

## Verified-not-done at session start

No branch atlas, no collision search beyond a single sector/`l`/overtone, no
Jordan-chain / Puiseux / monodromy machinery, no certification, no effective
model, no pseudospectral or time-domain analysis, no manuscript. Gates 5 and
7–15 were untouched.

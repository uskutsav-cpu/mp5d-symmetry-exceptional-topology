# Execution gates

Status vocabulary: COMPLETE / PARTIAL / NOT STARTED. A gate is COMPLETE only
when its acceptance conditions pass — never because code exists.

| # | Gate | Status | Evidence |
| --- | --- | --- | --- |
| 1 | Repository recovery and consistency | **COMPLETE** | state verified at `c2ab227f`; stale "9-term" claims corrected (see below) |
| 2 | Raw-recurrence asymptotic tail | **PARTIAL** | leading order only (`u₁ = −√(−2c)`), ~1.35× depth reduction; `u₂`+ not derived; raw-recurrence tail path not built |
| 3 | Independent Solver C | **COMPLETE** | `docs/SOLVER_C.md`; complex-scaled collocation, recurrence-free (enforced by test) |
| 4 | Three-solver validation | **COMPLETE** | 10/10 points agree < `1e-5`, incl. r₂ = 0.140 and 0.190 |
| 5 | Long-lived massive branches | **NOT STARTED** | — |
| 6 | Cross-sector no-go closure | **PARTIAL** | O6.1–O6.4 discharged; O6.5–O6.7 open |
| 7 | Fixed-sector branch atlas | **NOT STARTED** | — |
| 8 | Candidate collision localization | **NOT STARTED** | — |
| 9 | EP verification | **NOT STARTED** | — |
| 10 | Exceptional-locus continuation | **NOT STARTED** | — |
| 11 | Higher-order and junction search | **NOT STARTED** | — |
| 12 | Effective local theory | **NOT STARTED** | — |
| 13 | Pseudospectral and time-domain | **NOT STARTED** | — |
| 14 | Numerical certification | **NOT STARTED** | — |
| 15 | Final reproducibility freeze | **NOT STARTED** | — |

## Documentation corrections made in this session

* The raw recurrence width was stated as **9** in `docs/FIRST_QNM_VALIDATION.md`,
  `docs/CLAIM_LEDGER.md` (C20) and `src/mp5d/radial/qnm.py`. That is the value
  for the **static** case `a = b = 0` only; for generic two-spin parameters the
  width is **13**. Corrected in all four places.
* Solver C's absence is no longer stated anywhere as current.

## Gate 3 acceptance detail

Solver C independence is enforced structurally: a test greps the module for
`reduce_to_three_term`, `continued_fraction`, `hill_determinant`,
`wynn_epsilon`, `recurrence_row` and fails if any appears.

## What blocks gate 5 onward

Nothing technical — gates 5–15 were simply not reached in this session. The
critical-path blocker (no independent solver) is now removed, so the branch
atlas and everything downstream are unblocked for the first time.

# Session handoff

## Session 2 (2026-08-01) — bulk research takeover

### Completed

**Administrative**
* Recovery audit (`docs/RECOVERY_AUDIT.md`): entry state independently verified,
  56/56 tests passing, handoff accurate.
* History-aware publication audit over every reachable commit — zero
  credentials, zero PDFs, zero large blobs, zero personal paths. Repository is
  now **PUBLIC** with history intact (no rewrite, nothing to revoke).
* CI **green**: `tests and lint` + `secret scan`, minimal `contents: read`
  permissions, no `pull_request_target`, slow tests excluded from the per-commit
  path, plus a manual `research-validation.yml` for the full suite.
* `main` protected: PR required, both checks required, force pushes and
  deletions blocked, linear history, conversation resolution.
* Durable project skill at `.claude/skills/mp5d-bulk-research/`.

**Scientific**
* `docs/CROSS_SECTOR_EP_NO_GO.md` written properly: the no-go now rests on the
  direct-sum structure of the operator pencil **and its resolvent**, giving a
  simple (order-1) resolvent pole at a cross-sector coincidence and hence
  semisimplicity. O6.1–O6.4 discharged; **O6.5–O6.7 open**.
* Radial singular structure and exponents (session 1) unchanged and still CI-pinned.
* 67 tests passing.

### NOT completed — and this is the blocking item

**No QNM frequency has been computed. Zero benchmarks. Zero EP searches.**

Two radial formulations were built; both are documented failures with
diagnoses in `docs/FAILED_APPROACHES.md` (entries 3–5):

* **Collocation (Solver A) — structurally wrong, abandoned.** After peeling
  `e^{iΩr}`, the ingoing solution `e^{−2iΩr}` *decays* for `Im ω < 0`, so
  demanding boundedness does not exclude it. Matrix singular everywhere in the
  lower half plane. Not fixable by tuning.
* **Leaver/Hill (Solver B) — partially working, the real path forward.** The
  Leaver variable, prefactor, and polynomial-coefficient recovery all work and
  are self-validating (degrees 12/11/10, stable). The hard Hill truncation does
  **not** isolate roots: the last columns have reduced norm and produce spurious
  near-null vectors at the truncation edge.

### Exact next commands

The blocker is the truncation, not the construction. Fix it in
`src/mp5d/radial/leaver.py`:

1. Reduce the multi-term recurrence (bandwidth set by `deg A = 12`) to a
   three-term recurrence by Leaver's Gaussian elimination, row by row.
2. Replace `hill_matrix` / `smallest_singular_value` with the continued-fraction
   condition on the three-term recurrence, using the standard large-`n` tail
   estimate `a_{n+1}/a_n → 1 − √(c)/√n + …` for a rank-1 irregular point, so the
   tail is closed analytically instead of truncated.
3. Validate in this order, and do not skip: 5D Schwarzschild–Tangherlini
   (`a = b = 0`, `l = 0`) first, then convergence in `depth`, then the
   `(a, m₁) ↔ (b, m₂)` exchange identity, which is exact and is the most
   sensitive diagnostic available.

```bash
cd ~/mp5d-symmetry-exceptional-topology
PYTHONPATH=src .venv/bin/python -m pytest -q            # expect 67 passed
PYTHONPATH=src .venv/bin/python -c "
from mp5d.radial.leaver import LeaverProblem
from mp5d.geometry import MPGeometry
p = LeaverProblem(MPGeometry(a=0.0,b=0.0,M=1.0), 0.0, 0,0,0, depth=50)
print(p.smallest_singular_value(0.9-0.9j))"   # currently ~1e-8 everywhere: the bug
```

Only after a QNM is trustworthy do the benchmark gate
(`.claude/skills/mp5d-bulk-research/SCIENTIFIC_GATES.md`), then the branch atlas
in diagonal sectors `m₁ = m₂` first, then the targeted EP search.

### Do not repeat

* Fully symbolic sympy derivations with all parameters symbolic (fails, twice).
* Any method that selects QNMs by "peel the outgoing factor, demand bounded".
* Leaving the clearing factor's `u^k` in place (zeroes the leading recurrence rows).

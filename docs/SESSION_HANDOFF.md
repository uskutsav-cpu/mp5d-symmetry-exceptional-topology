# Session handoff

## Session 3 (2026-08-01) — radial solver rescue: SUCCESS

### Achieved

**First trustworthy QNM.** 5D Schwarzschild–Tangherlini, massless scalar,
`l = 0`, `n = 0`, `r₊ = 1`:

```
omega = 0.533835574268 - 0.383375368513 i
```

matching Matyjasek (arXiv:2107.04815, Table I) to all 12 published digits
(`4.98e-13`). Four benchmarks total (`l = 0,1` × `n = 0,1`), all cross-solver
verified. Full detail in `docs/FIRST_QNM_VALIDATION.md`.

**Two independent solvers.**
* Solver A — Leaver CF: width-9 Frobenius recurrence, general Gaussian
  reduction to three terms, `n`-th CF inversion.
* Solver B — Hill/Wynn: raw recurrence, `a_N` condition, Wynn acceleration.
  No reduction, no continued fraction; independent root condition.
They agree to `2.15e-12`.

**U(2) degeneracy confirmed numerically.** The session-1 analytic prediction
(claims C6/C7) now has numerical backing: `(1,1)`, `(2,0)`, `(0,2)` at `l=2`,
`m=2`, `a=b=0.25` share one frequency to `4.67e-15`; `δ ≠ 0` splits them by
`5.38e-02`. Exchange symmetry holds to `4.15e-15`.

93 tests passing.

### Not done, and needed before the science resumes

1. **Rotating benchmarks against published values.** The rotating results are
   validated by exact internal symmetries and A/B agreement only. The
   Huang–Huang two-spin comparison (arXiv:2502.11764) is still open, as are the
   singly-rotating literature points. This is the next gate.
2. **Precision escalation.** Everything is double precision; `A, B, C` are
   extracted by FFT, flooring residuals near `1e-12`. Gate requirement 3
   ("stable with arithmetic precision") is therefore **not discharged**. Move
   the extraction to exact rational or `mpmath` arithmetic.
3. **Nollert asymptotic tail.** Not implemented; the CF uses plain truncation.
   Fine for these modes, needed for high overtones and long-lived massive modes.

### Then, and only then

Branch atlas in diagonal sectors `m₁ = m₂` first (they have exact `δ`-parity),
then the targeted within-sector EP search. No exceptional-point work has been
done and none should start before gate 1 above.

### Reproduction

```bash
PYTHONPATH=src .venv/bin/python -m pytest -q            # 93 passed
PYTHONPATH=src .venv/bin/python scripts/first_qnm.py    # regenerates results/first_qnm.json
```

### Do not repeat

* Fully symbolic Leaver derivation with symbolic parameters (fails, three times).
* Peel-and-demand-boundedness selection for `Im ω < 0` (FAILED_APPROACHES #3).
* SVD of the equilibrated Hill matrix (#5/#5b — dynamic range, not edge modes).
* Widening the Muller initial triangle (#6 — made things worse).

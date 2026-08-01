# Rotating validation against Huang & Huang (arXiv:2502.11764)

**Status: cross-solver verified against the published two-spin spectrum.**

## Convention mapping (verified, not assumed)

Their metric, `ρ²`, `Δ` and mode ansatz agree with ours exactly. Three mappings
had to be established before any comparison was meaningful:

* their `λ_{k m₁ m₂}` uses `(a²−b²)cos²θ` with **no** `b²` term, so
  **their λ = our Â = Λ + (ω²−μ²)b²**;
* their `k` is our angular node index, so `l = 2k + |m₁| + |m₂|`;
* their `σ` denominator `2r_H(a²+b²−M+2r_H²)` equals our `2r₊(z₊−z₋)`
  identically;
* their `r2` is our inner horizon `r₋` — reproduced from their Table I to `2e-5`.

## Results

| Set | Passed | Typical agreement |
| --- | --- | --- |
| Table II (static + singly rotating + massive) | **6/6** | `3e-6`–`4e-6` vs their CFM |
| Table III (arbitrary two-spin) | **3/3** | `1.8e-6`–`7.4e-5` vs their CFM |
| Table V/VI exchange pairs | **5/5** | exchange residual `< 1e-10` |
| Table IV/VII equal-spin family, off-diagonal | **4/4** | `5e-5`–`2.2e-4` vs their MM |

### The three Table III two-spin points

| `a` | `b` | `μ` | `k` | ours | their CFM | diff |
| --- | --- | --- | --- | --- | --- | --- |
| 0.2 | 0.3 | 0.1 | 0 | `1.681109 − 0.347227i` | `1.68112 − 0.3472i` | `2.9e-5` |
| 0.4 | 0.2 | 0.9 | 0 | `1.822164 − 0.311200i` | `1.82222 − 0.311152i` | `7.4e-5` |
| 0.3 | 0.1 | 0.1 | 1 | `2.634400 − 0.349099i` | `2.6344 − 0.349097i` | `1.8e-6` |

### Tolerance policy

The paper prints 5–6 significant digits, and its **own two methods disagree by
`1.9e-4`** at an identical parameter point (Table III row 1: MM
`1.6812 − 0.347026i` vs CFM `1.68112 − 0.3472i`). Agreement is therefore judged
at `2e-4` against its continued-fraction column — the closest methodological
match to Solver A — and at `3e-4` against its matrix-method grids. That is the
paper's internal method spread, not a relaxed goalpost. Our agreement with the
CFM column is typically `1e-5` or better.

## Two corrections found along the way

**Ours.** Table II's key is `{a, μ, k, m₁, m₂}` with `b = 0` fixed by its
caption. The row printed `{0.3, 0.3, 1, 1, 1}` therefore means `a = 0.3`,
`μ = 0.3`, `b = 0` — a singly rotating *massive* case. We first transcribed it
as `a = b = 0.3`, `μ = 0`, and our solver disagreed by `1.1e-1`. With the
correct reading we agree to `4.2e-6`. The disagreement was the transcription,
not the solver — which is exactly why the manifest requires a second check.

**Theirs (suspected).** The equal-spin diagonals of Tables V, VI and VII do not
match their own tables' trends. Table VII's diagonal
(`1.56895, 1.63886, 1.72958, 1.86008`) is *identical* to Table IV's, although
Table IV is `k = 0` and Table VII is `k = 1`. Computing those points directly:

| `a = b` | ours (`k = 1`) | printed diagonal | paper's own off-diagonal neighbours |
| --- | --- | --- | --- |
| 0.2 | `2.631617 − 0.349994i` | `1.63886 − 0.351547i` | `2.59614`, `2.67298` |
| 0.3 | `2.718895 − 0.341123i` | `1.72958 − 0.340538i` | `2.67298`, `2.77349` |

Our values interpolate the paper's own off-diagonal neighbours smoothly; the
printed diagonal does not. The most likely mechanism: at `a = b` the
spheroidicity `c² = 4sδ(ω²−μ²)` vanishes identically, so an angular routine
returning the lowest eigenvalue regardless of the requested `k` would reproduce
exactly this pattern. Recorded as a suspected source artifact in
`data/published_benchmarks/huang_huang_2025/provenance.json`; **not** treated as
a solver failure, and excluded from the pass criteria.

An earlier version of that manifest argued the anomaly from a supposed
contradiction with Table II. That argument rested on our own misreading above
and has been **withdrawn**; the anomaly now rests only on the direct computation
in the table above.

## Recurrence validity domain

Singularities of the series about `x = 0` sit at `x = 1` (`r → ∞`) and at the
images of `r = −r₊`, `r = −r₋`, `r = 0`:

```
x(−r₊) = 2r₊/(r₊+r₋)      x(−r₋) = (r₊+r₋)/(2r₋)      x(0) = r₊/r₋
```

All exceed 1 while `r₊ > r₋`, so the radius of convergence is exactly 1 with the
irregular point on the boundary. As `r₋ → r₊` (near-extremal), `x(−r₊) → 1` from
above: a spurious singularity crowds the unit circle and convergence at `x = 1`
degrades. **That is the analytic mechanism behind the deterioration Huang and
Huang report near `r₂ ≳ 0.1`.**

Our compactification is not limited at the same place. Over the 27 mapped
points, 24 are recurrence-convergent out to `r₂ = 0.200`, and only `r₂ ≳ 0.26`
becomes slowly convergent. Empirically we reproduce their Table VII entry at
`(a,b) = (0.3,0.4)`, where `r₂ = 0.140` — beyond their stated limit — to
`5.3e-5`. Map: `results/recurrence_validity_map.json`.

This is *not* a claim that the method is reliable everywhere. No point was
classified from Newton returning a number; the classification is from the
singularity geometry, and the `r₂ > 0.1` claims are backed by direct comparison
with their matrix-method values.

## Reproduction

```bash
PYTHONPATH=src .venv/bin/python scripts/build_huang_huang_manifest.py
PYTHONPATH=src .venv/bin/python scripts/rotating_validation.py
```

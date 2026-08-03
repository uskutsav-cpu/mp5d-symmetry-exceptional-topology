# Novelty audit

Performed 2026-08-02, before any novelty language was used anywhere in the
repository. The rule applied: **no claim is described as new, first, or
universal**; each is recorded with its nearest prior result and the exact
distinction, and "priority uncertain" is used wherever an exhaustive check was
not possible.

## Established prior work — must not be presented as new

| result | source | status |
| --- | --- | --- |
| Massive-scalar QNMs of doubly rotating 5D Myers–Perry, two arbitrary spins; continued-fraction and matrix methods; long-lived modes at large scalar mass | Huang & Huang, [arXiv:2502.11764](https://arxiv.org/abs/2502.11764) (Phys. Lett. B) | **established**; used here as the baseline |
| Massive-scalar QNMs / bound states of 5D Myers–Perry | Li et al., [arXiv:2309.08203](https://arxiv.org/pdf/2309.08203) | **established** |
| Isometry enhancement at equal angular momenta, `R × U(1) × U(1) → R × SU(2) × U(1)` in 5D; cohomogeneity-1 metric | Murata & Soda, [arXiv:0803.1371](https://arxiv.org/pdf/0803.1371) | **established**. The enhancement itself is *not* claimed here. |
| **Exceptional lines** in massive-scalar Kerr–de Sitter and Myers–Perry–de Sitter, in the (scalar mass, spin) plane and the Nariai limit | Nakamoto & Oshita, [arXiv:2601.00704](https://arxiv.org/abs/2601.00704) | **established** |
| Exceptional points / eigenvalue repulsion in Kerr, Kerr–Newman | Cavalcante et al. 2407.20850; Dias et al. 2109.13949 | **established** |
| Jordan structure at exceptional points; EP framework for ringdown | Panosso Macedo et al., [arXiv:2512.02110](https://arxiv.org/html/2512.02110) | **established** |
| Complex-scaling approach to black-hole QNMs | [arXiv:2604.20442](https://arxiv.org/html/2604.20442) | **established**; Solver C is an instance of this general method |

## The critical contrast

Exceptional lines **do** exist in the de Sitter members of this family
(Kerr–dS, Myers–Perry–dS), including in the scalar-mass/spin plane. This
project studies the **asymptotically flat** (`Λ = 0`) doubly rotating case and
finds none in the searched domain.

That is a contrast between different spacetimes, **not** a contradiction, and it
must be stated that way. The de Sitter exceptional lines are tied to the
near-Nariai structure — the cosmological horizon and the Pöschl–Teller limit —
which has no analogue at `Λ = 0`. Nothing here challenges those results, and
nothing here is evidence about them.

## Claims made in this repository, with nearest prior result

| claim | nearest prior work | distinction | novelty status |
| --- | --- | --- | --- |
| C42 — cross-sector coincidences cannot be defective; Jordan chains of a direct sum are the union of the summands' | The direct-sum lemma is standard operator theory (Kato, *Perturbation Theory for Linear Operators*) | The **lemma is not new**. Its use to close the equal-spin `U(2)` multiplet channel for this system is what is claimed. | **priority uncertain** — elementary enough that it may be folklore; no source found stating it for this system |
| C39/C40 — Klein four-group `{1,E,P,EP}`; diagonal sectors even in `δ`, anti-diagonal even in `s`; `(0,0)` both | Exchange symmetry `(a,m₁)↔(b,m₂)` is standard | The **anti-diagonal even-in-`s`** statement and the stabilizer-based classification were not found in the literature | **priority uncertain**, not claimed as new |
| C43 — diagonal-sector exceptional lines must come in mirror pairs and terminate quadratically on `δ = 0` | Symmetry-protected exceptional chains in non-Hermitian crystals (Zhang et al. 2204.08052) — analogous mechanism in a different setting | Conditional and **untested**: no exceptional line was found to exhibit it | **conjectured (conditional)**; explicitly not a discovery |
| C37 — bounded multi-sector EP2 exclusion for asymptotically flat MP5D | No prior systematic EP search for this spacetime was found | A negative result over a declared domain | **priority uncertain**; stated as a bounded exclusion, never as "there are no EPs" |
| C44 — the near-extremal breakdown is the quasiresonant limit (`Re Ω → 0` ⇒ required scaling angle → 90°) | Complex-scaling QNM literature states the admissibility condition on `θ` | The condition is standard; identifying it as **the** cause of the failure previously attributed to extremality is a correction internal to this project | **internal correction**, not a literature claim |
| C41 — multiplet dimension `l+1`, level total `(l+1)²` | `(l+1)²` is the standard SO(4) scalar-harmonic level dimension | Only the `l+1`-at-fixed-`m` count and its use here are stated | **not new**; standard representation theory |
| Scale-free EP diagnostic `\|a₁/a₂\|` | Contour/moment root finding (Kravanja & Van Barel); standard Puiseux theory | The observation that `\|dF/dω\|` cannot bound anything, and the rescaling-invariant replacement | **priority uncertain**; presented as methodology, not discovery |

## Words not used, and why

`first`, `new`, `universal`, `topological`, `symmetry-protected`, `certified`,
`theorem`, `PRL-level` do not appear as descriptions of any result in this
repository, except:

* "theorem" is used for the **direct-sum Jordan lemma**, which is proved in two
  lines from standard facts and is not claimed as original;
* "symmetry derived" is used as an *evidence level* in the claim ledger, with
  the meaning defined in `STATUS_SCHEMA.md`, not as a novelty claim.

## Limits of this audit

Searches were run against the open literature in August 2026 for the specific
statements above. They are not exhaustive: in particular the direct-sum no-go is
elementary enough that it may appear in the black-hole perturbation literature
in a form these searches did not surface. Any priority claim should be
re-checked before submission, and none is asserted here.

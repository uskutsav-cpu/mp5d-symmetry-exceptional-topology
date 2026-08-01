# results/*.json schema

Every scientific record must carry:

| field | meaning |
| --- | --- |
| `mode_labels` | `{m1, m2, l, n, overtone}` |
| `parameters` | `{a, b, s, delta, mu, M}` |
| `solver` | solver name + version |
| `precision` | working precision (bits or dps) |
| `resolution` | radial/angular truncation |
| `residual` | root residual of the nonlinear system |
| `boundary_residual` | horizon/infinity boundary residuals |
| `convergence` | history of the root iteration |
| `evidence_level` | see the claim-ledger vocabulary |
| `independent_check` | which second solver agreed, and to what tolerance |
| `certification_level` | none / matrix / recurrence / algebraic / continuum |
| `commit` | git commit hash at generation time |
| `env_hash` | hash of `requirements.lock` |

Claim-ledger evidence vocabulary: literature-established, symbolically derived,
minimally reproduced, numerically observed, convergence verified, cross-solver
verified, symmetry derived, matrix certified, recurrence certified, continuum
proved, conjectured, falsified, unresolved.

Required files:

```
results/status.json
results/recovery_state.json
results/baseline_validation.json
results/equal_spin_branch_catalog.json
results/degeneracy_candidates.json
results/exceptional_candidates.json
results/verified_exceptional_points.json
results/negative_search_regions.json
results/certified_candidates.json
results/artifact_manifest.json
results/compute_manifest.json
```

Provenance rule: `commit` must be a real hash. If an artifact is generated
before its commit exists, regenerate it afterwards rather than leaving
`"uncommitted"`.

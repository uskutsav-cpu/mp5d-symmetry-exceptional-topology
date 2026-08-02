# Evidence vocabulary

Ordered roughly by strength. Never use a stronger label than the evidence.

| status | means |
| --- | --- |
| `literature-established` | published elsewhere; cite, do not claim |
| `symbolically derived` | exact symbolic derivation, test-enforced |
| `symmetry derived` | follows from an exact symmetry of the system |
| `minimally reproduced` | published number reproduced to stated tolerance |
| `numerically observed` | computed once, no convergence study |
| `convergence verified` | stable under resolution/depth/precision refinement |
| `cross-solver verified` | two **independent** formulations agree |
| `matrix certified` | interval/ball certificate for the finite matrix |
| `truncated-recurrence certified` | ditto for the truncated recurrence |
| `algebraic-system certified` | ditto for the nonlinear algebraic system |
| `continuum proved` | statement about the differential operator itself |
| `conjectured` / `falsified` / `unresolved` | self-explanatory |

Forbidden without the matching evidence: *first*, *new*, *universal*,
*topological*, *symmetry-protected*, *certified*, *theorem*, *PRL-level*.

Error separation is mandatory for any certification claim:

    E_total <= E_discretization + E_truncation + E_arithmetic + E_root

A certificate for a truncated object is never a continuum theorem.

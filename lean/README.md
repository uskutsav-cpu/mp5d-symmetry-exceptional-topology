# Exact algebra layer — build NOT verified in the delivery environment

Pinned to mathlib commit `e37d88a26f3791ed5a93daa1f949af1021b8d103`
(2026-09-09) and its declared Lean `v4.34.0-rc2` toolchain. No floating-point
calculation, QNM table, or continuum spectral claim is formalized here.

```
cd lean
lake update
lake exe cache get
lake build
```

There are no admitted propositions in these source files. That is not a claim
that compilation has succeeded: Lean/Lake were unavailable in the execution
environment, and these proof drafts must be kernel checked. The theorem about
`Fin (ell+1)` counts an index set; the full bijection to allowed azimuthal pairs
is deliberately NOT represented as already proved. The operator results are
algebraic linear-map/product results, not continuum-QNM domain theorems.

Keep this entire layer outside active manuscript claims until the build has
passed and the precise theorem scope has been reviewed. Do not generate a
`lean.json` success receipt by hand.

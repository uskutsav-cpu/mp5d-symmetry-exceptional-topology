# Scientific evidence gates

## Benchmark gate (must pass before any wide EP scan)

1. one 5D Schwarzschild-Tangherlini scalar QNM;
2. two singly rotating 5D points;
3. three general two-spin points;
4. one equal-spin point;
5. one massive long-lived / quasiresonant point;
6. one exact spin-exchange pair `(a,m1) <-> (b,m2)`;
7. at least one nontrivial point computed by **both** radial solvers.

Engineering target: relative error `1e-6` where the source prints enough
digits. Do not declare failure only because a paper prints few digits; record
the published precision alongside the comparison.

## EP2 gate — a candidate is NOT an EP2 until all hold

* converged frequency coalescence under resolution AND precision increase;
* converged angular-eigenvalue behaviour;
* radial eigenfunction coalescence;
* angular eigenfunction coalescence where applicable;
* algebraic multiplicity 2, geometric multiplicity 1;
* an explicit generalized eigenvector / Jordan chain with a residual bound;
* local square-root (Puiseux) behaviour of the branches;
* 2-cycle mode permutation around a closed parameter loop (monodromy);
* agreement of two genuinely independent solvers;
* stability under resolution and precision refinement.

Two nearby frequencies are never sufficient. Near-defectiveness is not
defectiveness.

## EP3 gate

Everything in the EP2 gate, plus: three converged branches, algebraic
multiplicity 3 with geometric multiplicity < 3, Jordan-chain evidence of length
3, cube-root (or symmetry-permitted) local scaling, 3-cycle monodromy, and
independent numerical confirmation. A three-mode cluster on a plot is not an
EP3.

## Certification gate

Distinguish, always and explicitly, certification of:

1. a finite-dimensional collocation matrix;
2. a truncated recurrence;
3. an algebraic determinant system;
4. the continuum differential operator.

Levels 1-3 never imply level 4. Attempt the error separation

```
E_total <= E_discretization + E_truncation + E_arithmetic + E_root
```

and report which terms are bounded rigorously vs. estimated. A failed
certification must be recorded with its exact obstruction.

## Independence rule

Two wrappers around the same discretized matrix are NOT independent. Acceptable
pairs: collocation vs. continued fraction; Chebyshev vs. Bernstein; collocation
vs. Hill determinant; direct nonlinear eigenproblem vs. determinant root;
distinct language implementations.

## Forbidden vocabulary without matching evidence

"first", "new", "universal", "topological", "symmetry-protected", "certified",
"theorem", "PRL-level".

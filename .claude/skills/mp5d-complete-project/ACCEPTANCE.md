# Acceptance criteria and stopping rule

## Administrative

* repository safely public (history audited);
* CI active, or the exact external blocker documented;
* `main` protected where the plan supports it;
* clean working tree, current branch pushed.

## Solver

* Radial Solver A production-ready, with the coupled dependence
  `Lambda = Lambda(omega, a, b, mu, m1, m2, l)` treated self-consistently
  during root iteration -- NOT as a frozen constant;
* first trustworthy QNM computed (residual + convergence recorded);
* published baseline set reproduced (see SCIENTIFIC_GATES.md);
* Solver B independently verifies nontrivial two-spin cases.

## Mathematical

* cross-sector no-go proof obligations resolved or explicitly left open;
* diagonal-sector even-`delta` claim proved, corrected, or falsified;
* branch-identity framework implemented (overlap-based, not sort-based).

## Research

* bounded equal-spin branch atlas completed;
* candidate interactions ranked by real diagnostics;
* targeted symmetry-breaking search executed;
* EP candidates verified or rejected against the EP2 gate;
* negative search regions recorded;
* strongest candidate subjected to a certification attempt.

## Stopping rule

Stop only on one of:

* **A** verified exceptional structure, independently confirmed;
* **B** substantial negative classification over the bounded target region,
  with candidates rejected by validated methods;
* **C** a demonstrated technical blocker (must be shown, not asserted, after
  multiple justified approaches);
* **D** a resource ceiling, with all cheaper work finished, checkpoints pushed,
  and the exact next commands written down.

Finding no EP is an acceptable outcome. Claiming one without the gate is not.

# CLAUDE.md — mp5d-symmetry-exceptional-topology

Repository-level rules. Full instructions live in the project skill
`.claude/skills/mp5d-complete-project/`.

* **Inspect real repository state before acting.** Never trust a summary,
  including a handoff document. Run the tests first.
* **Continue through the task queue.** Do not stop after one subtask, one
  computed value, or repository setup.
* **Never fabricate numerical or bibliographic results.** An unverified
  bibliographic entry is "indexed", not "verified".
* **Run `pytest -q` before every commit.**
* **Record failed methods** in `docs/FAILED_APPROACHES.md` so they are not
  rediscovered.
* **Preserve checkpoints.** Expensive work must be resumable.
* **Distinguish evidence levels** using the claim-ledger vocabulary in
  `.claude/skills/mp5d-complete-project/STATUS_SCHEMA.md`.
* **Never call a close frequency pair an exceptional point.** Every applicable
  condition in `docs/EP_VERIFICATION_GATE.md` must pass.
* **Never use `|dF/domega|` as a bound** — it is not invariant under rescaling
  of the spectral condition. Use `root_separation` or the branch gap.
* **Never label a branch by sorting frequencies.** Identity comes from
  continuation history plus predictor agreement.
* **Never claim continuum certification from a matrix truncation.**
* **Stay inside the validated domain** (`z_minus <= 0.20`, extremality
  `>= 0.05`); `z_minus` is the *inner* horizon radius squared.
* **Work on a research branch**, not `main`.
* **Push every validated milestone.**

## Conventions

`docs/CONVENTIONS.md` is authoritative and is re-derived symbolically on every
test run. Code must match it, not a paper.

Scalar mass is `mu`. The Myers-Perry mass parameter is `M`. Never conflate them.
`Phi ~ e^{-i omega t}`, so damped modes have `Im omega < 0`.

## Environment

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
PYTHONPATH=src .venv/bin/python -m pytest -q
```

`uv` and Julia are not installed on the development machine; `requirements.lock`
is a `pip freeze` and substitutes for `uv.lock`.

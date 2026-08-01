# CLAUDE.md — mp5d-symmetry-exceptional-topology

Repository-level rules. Full instructions live in the project skill
`.claude/skills/mp5d-bulk-research/`.

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
  `.claude/skills/mp5d-bulk-research/STATUS_SCHEMA.md`.
* **Never call a close frequency pair an exceptional point.** See the EP2 gate.
* **Never claim continuum certification from a matrix truncation.**
* **Work on `research/core-science`**, not `main`.
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

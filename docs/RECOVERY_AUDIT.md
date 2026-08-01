# Recovery audit

Performed at the start of the bulk-research takeover, before any modification.

## Verified state

| Item | Finding |
| --- | --- |
| Entry commit | `bf67e35`, matching the handoff |
| Working tree | clean |
| Branches | `main`, `research/core-science` (both pushed) |
| Tags | `milestone/stage-A`, `milestone/radial-structure` |
| Test suite | **56 passed, 0 failed** — re-run before any edit |
| Python | 3.13.13 (project venv) |
| Dependencies | numpy 2.5.1, scipy 1.18.0, sympy 1.14.0, mpmath 1.3.0, python-flint 0.9.0 |
| Missing tooling | `uv`, `julia`, `gitleaks`, `trufflehog`, SLEPc/PETSc |

The handoff summary was accurate in every checked particular. One thing had
**changed** since it was written: the GitHub token now carries the `workflow`
scope, so CI could be activated in this session.

## Claimed artifacts confirmed to exist in code

* symbolic separation derivation and its 5 regression tests;
* geometry module (horizons, `Ω_a`, `Ω_b`, `κ`, `T_H`, `s`/`δ`);
* exact angular Jacobi solver;
* independent angular finite-difference solver;
* radial *structure* module (characteristic exponents only).

Confirmed **absent**, as claimed: any radial eigenvalue solver, any computed QNM
frequency, any numerical exceptional-point result.

## Publication audit

`gitleaks` and `trufflehog` are not installed. Substituted `detect-secrets`
(27 plugins) run over the extracted file tree of **every reachable commit**,
plus direct pattern and object-level scans of full history.

| Check | Result |
| --- | --- |
| Credential patterns across all commits | none |
| `.env` / `.pem` / `.key` / token files ever committed | none |
| PDFs in history | none |
| Blobs > 1 MB | none |
| Files > 25 MB on disk | none |
| Personal absolute paths (`/Users/...`) in history | none |
| detect-secrets findings | 1, a false positive |

The single finding is `lock_sha256_16` in `results/symmetry_multiplets.json`,
flagged as a high-entropy hex string. It is the SHA-256 provenance hash of
`requirements.lock`. Not a secret.

**Verdict: safe to publish.** No credential requires revocation or rotation.

## License

MIT, committed at repository creation, applying to first-party code only. No
third-party source has been vendored, so there is no compatibility conflict. No
copyrighted article text or PDF is committed; the bibliography holds metadata
only.

## Corrections made as a result of this audit

* `results/symmetry_multiplets.json` had provenance `"commit": "uncommitted"`
  because it was generated before its commit existed. Regenerated so provenance
  is a real hash, and `STATUS_SCHEMA.md` now forbids the placeholder.

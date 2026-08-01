# Pending CI workflow

`ci.yml.pending` is the intended `.github/workflows/ci.yml`.

It is parked here because the GitHub token in use has scopes
`gist, read:org, repo` but **not** `workflow`, so pushing a file under
`.github/workflows/` is rejected by the remote.

To activate:

```bash
gh auth refresh -h github.com -s workflow
git mv ci/ci.yml.pending .github/workflows/ci.yml
git commit -m "Enable CI" && git push
```

Until then CI does not run on the remote; tests are run locally with
`PYTHONPATH=src .venv/bin/python -m pytest -q`.

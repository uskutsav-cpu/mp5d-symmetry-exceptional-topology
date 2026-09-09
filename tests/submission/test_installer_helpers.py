"""The installer is available in the delivery bundle, not inside an applied repo.
Its structural helper tests are executed through the bundle's own test script.
This file tests the safety-related source properties that remain after applying.
"""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_strict_workflow_has_no_error_swallowing():
    text = (ROOT / "scripts/reproduce_submission_science.sh").read_text()
    assert "set -euo pipefail" in text
    assert "|| echo" not in text and "|| true" not in text
    assert '-m "not slow"' not in text
    assert "--strict-markers" in text


def test_all_new_python_sources_parse():
    for path in (ROOT / "src/mp5d_science").glob("*.py"):
        ast.parse(path.read_text(), filename=str(path))


def test_no_unproved_lean_placeholders():
    import re

    for path in (ROOT / "lean").rglob("*.lean"):
        assert not re.search(r"\b(sorry|admit|axiom)\b", path.read_text())


def test_lean_toolchain_is_pinned_not_master():
    assert "e37d88a26f3791ed5a93daa1f949af1021b8d103" in (ROOT / "lean/lakefile.toml").read_text()

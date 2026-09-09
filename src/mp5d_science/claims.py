"""Claim registry, withdrawn-diagnostic checks and exact-commit freeze gates."""

from __future__ import annotations

from pathlib import Path

from .evidence import WITHDRAWN
from .provenance import check_result, file_hash, read_json, source_snapshot


def withdrawn_paths(value, path="$", *, historical=False):
    """Only explicitly historical/WITHDRAWN records are exempt, never active gates."""
    problems = []
    if isinstance(value, dict):
        local_historical = historical or value.get("authority") in {"HISTORICAL_ONLY", "WITHDRAWN"}
        if value.get("active") is True:
            local_historical = False
        for key, item in value.items():
            here = f"{path}.{key}"
            if key in WITHDRAWN and not local_historical:
                problems.append(here)
            if key == "evidence_ids" and isinstance(item, list) and not local_historical:
                problems.extend(f"{here}[{i}]" for i, x in enumerate(item) if x in WITHDRAWN)
            problems.extend(withdrawn_paths(item, here, historical=local_historical))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            problems.extend(withdrawn_paths(item, f"{path}[{i}]", historical=historical))
    return problems


def check_claim_registry(registry):
    errors = []
    if registry.get("schema_version") != 1:
        errors.append("claim schema must be version 1")
    claims = registry.get("claims")
    if not isinstance(claims, list) or not claims:
        return errors + ["nonempty claim registry required"]
    ids = [c.get("id") for c in claims]
    if len(set(ids)) != len(ids) or any(not x for x in ids):
        errors.append("claim IDs must be unique and nonempty")
    for claim in claims:
        if claim.get("active"):
            if claim.get("status") not in {"SUPPORTED", "NUMERICALLY_SUPPORTED"}:
                errors.append(claim["id"] + ": active claim lacks supported status")
            if not claim.get("evidence_files"):
                errors.append(claim["id"] + ": evidence files missing")
            if claim.get("scope") not in {
                "FINITE_ALGEBRA",
                "FINITE_DISCRETIZATION",
                "BOUNDED_NUMERICAL_SCAN",
            }:
                errors.append(claim["id"] + ": unsupported or unspecified scope")
            for path in withdrawn_paths(claim):
                errors.append(claim["id"] + ": withdrawn diagnostic at " + path)
    return errors


def validate_freeze(root: Path, artifact_dir: Path, expected_commit: str, registry: dict) -> dict:
    """This gate cannot promote a missing, development, historical, or failed result."""
    errors = check_claim_registry(registry)
    if registry.get("science_ready") is not True:
        errors.append("claim registry explicitly marks science as not ready")
    if not any(c.get("active") for c in registry.get("claims", [])):
        errors.append("there are no active supported submission claims")
    expected_digest = source_snapshot(root)["sha256"]
    artifacts = {}
    required = registry.get("required_artifacts", {})
    if not required:
        errors.append("required artifact set is empty")
    for claim in registry.get("claims", []):
        if claim.get("active"):
            for name in claim.get("evidence_files", []):
                if name not in required:
                    errors.append(
                        f"{claim['id']}: evidence file is not a required checked artifact: {name}"
                    )
    base = Path(artifact_dir).resolve()
    for name, contract in required.items():
        path = (base / name).resolve()
        if not path.is_relative_to(base):
            errors.append(f"{name}: unsafe artifact path")
            continue
        if not path.is_file():
            errors.append(f"{name}: required artifact missing")
            continue
        try:
            result = read_json(path)
            errors.extend(
                f"{name}: {e}"
                for e in check_result(
                    result,
                    expected_commit=expected_commit,
                    expected_source_digest=expected_digest,
                    authoritative=True,
                )
            )
            if result.get("kind") != contract["kind"]:
                errors.append(f"{name}: wrong result kind")
            if result.get("payload", {}).get("status") not in contract["accepted_statuses"]:
                errors.append(f"{name}: scientific status not accepted")
            errors.extend(
                f"{name}: active withdrawn diagnostic at {p}"
                for p in withdrawn_paths(result.get("payload"))
            )
            artifacts[name] = file_hash(path)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            errors.append(f"{name}: {exc}")
    if registry.get("quasiresonant_scope") not in {"EXCLUDED", "VALIDATED"}:
        errors.append("quasiresonant scope must be explicitly excluded or independently validated")
    return {
        "status": "SCIENCE_FREEZE_PASS" if not errors else "BLOCKED",
        "errors": errors,
        "source_commit": expected_commit,
        "source_tree_sha256": expected_digest,
        "artifact_hashes": artifacts,
        "continuum_qnm_theorem": False,
    }

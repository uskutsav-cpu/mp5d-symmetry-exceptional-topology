"""Installer safety tests on temporary, synthetic Git repositories.

These do NOT claim an integration run on the full upstream scientific repo.
The audited baseline SHA is replaced only within these isolated test fixtures.
"""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('bundle_installer', HERE / 'APPLY_UPDATE.py')
installer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(installer)


def command(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args], text=True).strip()


@pytest.fixture
def repository(tmp_path, monkeypatch):
    root = tmp_path / 'repo'
    root.mkdir()
    files = {
        'src/mp5d/continuation/branch.py':
            'def sample(hist, ts, arc):\n    return extrapolate(hist, ts[:-1] + [arc], arc)\n',
        'src/mp5d/exceptional/diagnostics.py':
            'def root_separation(*args, **kwargs):\n    return 1\n\ndef zero_count(*args, **kwargs):\n    return 1\n',
        'src/mp5d/exceptional/verification.py':
            '\n\n'.join('def '+name+'(*args, **kwargs):\n    return 1' for name in
                         ('roots_in_disc', 'track_roots_on_loop', 'jordan_chain_matrix'))+'\n',
        'src/mp5d/radial/qnm.py': '''def solve_qnm_cf(a, b, field_mass, m1, m2, ell, M=1,
                 depth_schedule=(100,200), angular_N=24, overtone=0):
    omega = 1j
    ok = True
    return QNMSolution(
        converged=ok,
    )
''',
        'src/mp5d/radial/solver_c.py': '''def solve_qnm_c(a, b, mu, m1, m2, ell, initial_frequency,
                resolution=160, angular_N=24, L=40, theta=1, M=1, tol=1e-10, maxiter=60):
    return None
''',
        'tests/unit/test_ep_machinery_synthetic.py': 'def test_old():\n    assert True\n',
        'pyproject.toml': '[tool.hatch.build.targets.wheel]\npackages = ["src/mp5d"]\n',
        'README.md': '# Historical fixture\n',
        'results/status.json': '{"status":"historical-complete","original_commit":"preserved"}\n',
        '.github/workflows/research-validation.yml': 'name: historical fixture\n',
    }
    for name, content in files.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    command(root, 'init', '-q')
    command(root, 'add', '.')
    command(root, '-c', 'user.name=Bundle fixture', '-c', 'user.email=fixture@example.test',
            'commit', '-qm', 'isolated installer fixture')
    command(root, 'switch', '-qc', 'review-test')
    command(root, 'remote', 'add', 'origin', 'https://github.com/'+installer.REPOSITORY+'.git')
    monkeypatch.setattr(installer, 'BASE', command(root, 'rev-parse', 'HEAD'))
    return root, files


def test_preflight_is_nonmutating_and_generates_guards(repository):
    root, original = repository
    changes = installer.build_changes(root)
    assert all((root / name).read_text() == content for name, content in original.items())
    assert command(root, 'status', '--porcelain') == ''
    assert b'ts[:-1], arc' in changes['src/mp5d/continuation/branch.py']
    assert b'_root_guard["valid"]' in changes['src/mp5d/radial/qnm.py']
    assert b'from mp5d_science.physical import solve_c' in changes['src/mp5d/radial/solver_c.py']
    status = json.loads(changes['results/status.json'])
    assert status['status'] == 'SCIENTIFIC_REVALIDATION_REQUIRED'
    assert status['historical_status']['original_commit'] == 'preserved'


def test_dirty_checkout_is_refused(repository):
    root, _ = repository
    (root / 'README.md').write_text('uncommitted work\n')
    with pytest.raises(ValueError, match='not clean'):
        installer.build_changes(root)


def test_protected_branch_is_refused(repository):
    root, _ = repository
    command(root, 'branch', '-M', 'main')
    with pytest.raises(ValueError, match='new working branch'):
        installer.build_changes(root)


def test_wrong_repository_is_refused(repository):
    root, _ = repository
    command(root, 'remote', 'set-url', 'origin', 'https://github.com/example/unrelated.git')
    with pytest.raises(ValueError, match='origin is not'):
        installer.build_changes(root)


def test_wrong_revision_is_refused(repository, monkeypatch):
    root, _ = repository
    monkeypatch.setattr(installer, 'BASE', '0'*40)
    with pytest.raises(ValueError, match='exact audited'):
        installer.build_changes(root)


def test_function_replacement_preserves_signature_and_neighbors():
    source = 'def f(a,\n      b=2):\n    """old"""\n    return a+b\n\ndef g():\n    return 3\n'
    result = installer.replace_function_body(source, 'f', 'raise RuntimeError("retired")')
    assert result.startswith('def f(a,\n      b=2):')
    assert 'def g():\n    return 3' in result
    with pytest.raises(ValueError):
        installer.replace_once('abcabc', 'abc', 'def')


def test_payload_tamper_is_refused(tmp_path, monkeypatch):
    root = tmp_path / 'bundle'
    (root / 'overlay').mkdir(parents=True)
    (root / 'overlay/x.py').write_text('changed=1\n')
    (root / 'PAYLOAD_SHA256.json').write_text('{"overlay_files":{"x.py":"bad"}}')
    monkeypatch.setattr(installer, 'HERE', root)
    with pytest.raises(ValueError, match='hash mismatch'):
        installer.overlay_files()


def test_transaction_rolls_back_content_on_write_error(repository, monkeypatch):
    root, original = repository
    write = installer.atomic_write
    attempts = 0
    def failing_once(path, data):
        nonlocal attempts
        attempts += 1
        if attempts == 3:
            raise OSError('injected write failure')
        return write(path, data)
    monkeypatch.setattr(installer, 'atomic_write', failing_once)
    monkeypatch.setattr(sys, 'argv', ['APPLY_UPDATE.py', str(root)])
    assert installer.main() == 2
    assert all((root / name).read_text() == content for name, content in original.items())
    assert command(root, 'status', '--porcelain') == ''


def test_new_shell_script_is_executable(tmp_path):
    target = tmp_path / 'script.sh'
    installer.atomic_write(target, b'#!/bin/sh\nexit 0\n')
    assert target.stat().st_mode & 0o111 == 0o111

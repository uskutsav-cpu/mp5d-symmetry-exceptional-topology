from copy import deepcopy
from pathlib import Path
import json
import pytest
from mp5d_science.provenance import envelope,check_result,digest,atomic_json,read_json
from mp5d_science.claims import withdrawn_paths,check_claim_registry,validate_freeze


def development(tmp_path):
    (tmp_path/'src').mkdir(exist_ok=True)
    (tmp_path/'src/source.py').write_text('x=1\n')
    return envelope('test',{'status':'PASS'},tmp_path,allow_development=True)

def test_development_is_valid_but_not_authoritative(tmp_path):
    record=development(tmp_path)
    assert not check_result(record)
    assert check_result(record,expected_commit='a'*40,authoritative=True)

def test_altered_artifact_is_detected(tmp_path):
    record=development(tmp_path);record['payload']['status']='CERTIFIED'
    assert 'integrity hash mismatch' in check_result(record)

def test_ancestor_is_not_freeze_commit(tmp_path):
    record=development(tmp_path)
    record['provenance'].update(source_commit='a'*40,dirty_checkout=False,authority='COMMITTED_SOURCE_RUN')
    record['integrity_sha256']=digest({k:v for k,v in record.items() if k!='integrity_sha256'})
    assert not check_result(record,expected_commit='a'*40,authoritative=True)
    assert any('exact' in x for x in check_result(record,expected_commit='b'*40,authoritative=True))

def test_atomic_write_rejects_nan_without_destroying_old_file(tmp_path):
    path=tmp_path/'test.json';atomic_json(path,{'ok':1})
    with pytest.raises(ValueError):atomic_json(path,{'bad':float('nan')})
    assert read_json(path)=={'ok':1}

@pytest.mark.parametrize('key',['root_separation','min_root_separation','a1_over_a2','min_abs_discriminant'])
def test_withdrawn_active_diagnostic_fails(key):
    assert withdrawn_paths({'active':True,key:.2})
    assert not withdrawn_paths({'authority':'WITHDRAWN',key:.2})
    assert withdrawn_paths({'authority':'WITHDRAWN','child':{'active':True,key:.2}})

def test_evidence_ids_checked_too():
    assert withdrawn_paths({'evidence_ids':['discriminant_independent']})

def test_claim_registry_does_not_upgrade_unknown_scope():
    registry={'schema_version':1,'claims':[{'id':'x','active':True,'status':'SUPPORTED','scope':'CONTINUUM_QNM','evidence_files':['a']}]}
    assert check_claim_registry(registry)

def test_freeze_requires_every_artifact_and_active_supported_claim(tmp_path):
    registry={'schema_version':1,'science_ready':True,'quasiresonant_scope':'EXCLUDED','claims':[
        {'id':'x','active':True,'status':'SUPPORTED','scope':'FINITE_ALGEBRA','evidence_files':['a.json']}],
        'required_artifacts':{'a.json':{'kind':'test','accepted_statuses':['PASS']}}}
    result=validate_freeze(tmp_path,tmp_path/'missing','a'*40,registry)
    assert result['status']=='BLOCKED'
    assert any('missing' in e for e in result['errors'])


def test_active_claim_cannot_reference_unchecked_evidence(tmp_path):
    registry = {
        'schema_version': 1, 'science_ready': True, 'quasiresonant_scope': 'EXCLUDED',
        'claims': [{'id': 'unsupported-reference', 'active': True, 'status': 'SUPPORTED',
                    'scope': 'FINITE_ALGEBRA', 'evidence_files': ['unchecked.json']}],
        'required_artifacts': {'checked.json': {'kind': 'test', 'accepted_statuses': ['PASS']}},
    }
    result = validate_freeze(tmp_path, tmp_path / 'artifacts', 'a' * 40, registry)
    assert result['status'] == 'BLOCKED'
    assert any('not a required checked artifact' in message for message in result['errors'])

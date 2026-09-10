
import pytest

from mp5d_campaign.core import Blocked, Store, read_json
from mp5d_campaign.release import _global_claims_safe, accepted_evidence, check_physical_row
from mp5d_campaign.audits import lean_build

from copy import deepcopy
from mp5d_campaign.core import AXES

def plan():
    return {'axes': deepcopy(AXES), 'higher_modes': {'overtones': [4,5], 'ell_offsets':[0,2,4]},
            'formal_domain_proposal': {'overtones':[0,1,2,3]},
            'adaptive_boundary': {'independent_continuation_directions':['forward','reverse']}}

def source(tmp_path):
    root=tmp_path/'src';root.mkdir();(root/'hello.py').write_text('x = 1\n');return root



def test_global_claim_stays_inactive():
    assert _global_claims_safe({'claims':[{'id':'global-no-ep','active':False,'status':'NOT_ESTABLISHED'}]})
    assert not _global_claims_safe({'claims':[{'id':'global-no-ep','active':True,'status':'NOT_ESTABLISHED'}]})
    assert not _global_claims_safe({'claims':[{'id':'global-no-ep','active':False,'status':'PROVED'}]})


def test_development_never_releases(tmp_path):
    with Store(source(tmp_path), tmp_path/'run', {}, plan(), development=True) as store:
        with pytest.raises(Blocked, match='Development evidence'):
            accepted_evidence(store, {}, plan())


def test_status_string_without_physics_is_rejected():
    with pytest.raises(Blocked):
        check_physical_row({'status':'PASS'}, 1e-6)


def test_missing_lake_is_not_build_pass(tmp_path, monkeypatch):
    root=source(tmp_path)
    (root/'lean').mkdir()
    (root/'lean/Test.lean').write_text('theorem trivial_test : True := True.intro\n')
    with Store(root,tmp_path/'run',{},plan(),development=True) as store:
        monkeypatch.setattr('mp5d_campaign.audits.shutil.which',lambda name:None)
        result=lean_build(store,{'lean_timeout_seconds':1})
        assert result['status']=='NOT_BUILT'
        assert result['compilation_checked'] is False
        assert read_json(store.output/'lean.json')['payload']['status']!='BUILD_PASS'

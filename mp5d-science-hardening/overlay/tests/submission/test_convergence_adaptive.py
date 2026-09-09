from dataclasses import replace
import pytest
from mp5d_science.convergence import Resolution,Measurement,run_ladder,axis_ladders,robustness_labels
from mp5d_science.adaptive import GapObservation,refine_minimum

BASE=Resolution()

def measured(r):
    return Measurement((1+1/r.cf_depth**4-1j,2-1j),(1e-12,1e-12),'test','test',r,True)

def test_convergence_all_rungs_required():
    levels=[replace(BASE,cf_depth=n) for n in (80,160,320)]
    assert run_ladder(measured,levels)['status']=='PASS'

def test_fake_precision_is_rejected():
    def fake(r):return replace(measured(r),effective_resolution=BASE)
    out=run_ladder(fake,[replace(BASE,precision_dps=n) for n in (30,50,80)])
    assert out['status']=='INCONCLUSIVE'
    assert len(out['failures'])==3

def test_one_failed_rung_cannot_be_hidden():
    def sometimes(r):
        if r.cf_depth==160:raise ArithmeticError('unstable pivot')
        return measured(r)
    out=run_ladder(sometimes,[replace(BASE,cf_depth=n) for n in (80,160,320,640)])
    assert out['status']=='INCONCLUSIVE'

def test_unresolved_pair_is_not_separated():
    def collapsed(r):return replace(measured(r),frequencies=(1-1j,1+1e-11-1j))
    assert run_ladder(collapsed,[replace(BASE,cf_depth=n) for n in (80,160,320)])['status']=='INCONCLUSIVE'

def test_unknown_axis_fails():
    with pytest.raises(ValueError):axis_ladders(BASE,{'fake_axis':[1,2,3]})

@pytest.mark.parametrize('axis',list(BASE.__dict__))
def test_every_axis_is_representable(axis):
    assert axis in axis_ladders(BASE,{axis:[getattr(BASE,axis)]})

def test_robustness_not_silently_added_to_formal_box():
    rows=robustness_labels(2,0,4)
    assert {(r['ell'],r['overtone']) for r in rows}=={(l,n) for l in (4,6,8) for n in (4,5)}
    assert all(r['scope']=='OUT_OF_DOMAIN_ROBUSTNESS_PROBE' for r in rows)

def good(x):return GapObservation(.5+x[0]**2+x[1]**2,1e-12,('a','b'),True,('forward','reverse'))

def test_adaptive_interior_minimum():
    result=refine_minimum(good,[.25,.25],[(-1,1)]*2,initial_radius=[.25,.25],max_levels=18)
    assert result['status']=='LOCALLY_STABLE_NUMERICAL_MINIMUM'
    assert result['gap']==pytest.approx(.5)
    assert not result['global_exclusion_proved']

def test_adaptive_boundary_minimum():
    result=refine_minimum(lambda x:good((x[0]-2,x[1])),[.5,0],[(0,1),(-1,1)],initial_radius=[.5,.5],max_levels=18)
    assert result['point'][0]==1
    assert result['boundary_active'][0]

def test_adaptive_missing_independent_directions_never_passes():
    result=refine_minimum(lambda x:replace(good(x),continuation_directions=('forward',)),[0,0],
                          [(-1,1)]*2,initial_radius=[.1,.1])
    assert result['status']=='INCONCLUSIVE'
    assert result['failures']

def test_adaptive_budget_exhaustion_is_failure():
    result=refine_minimum(good,[.1,.1],[(-1,1)]*2,initial_radius=[.1,.1],max_evaluations=2)
    assert result['status']=='INCONCLUSIVE'

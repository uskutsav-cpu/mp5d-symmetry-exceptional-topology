import numpy as np
import pytest
from mp5d_science.evidence import *


def witnesses(omega=0j, lineage='contour'):
    return [SolverEvidence('A','recurrence',omega,1e-14,1e-12,True,True),
            SolverEvidence('C',lineage,omega+1e-13,1e-14,1e-12,True,True)]


def test_jordan_ep2():
    T=np.array([[0,1],[0,0]],complex)
    r=jordan_test(T,-np.eye(2),np.zeros((2,2)))
    assert r['nullity']==1 and r['order_two_supported']
    assert not r['continuum_certificate']


@pytest.mark.parametrize('scale',[1e-4,1.,1e4])
def test_jordan_scale_invariance(scale):
    r=jordan_test(scale*np.array([[0,1],[0,0]]),-scale*np.eye(2),np.zeros((2,2)))
    assert r['order_two_supported']


def test_zero_matrix_is_not_defective_ep2():
    r=jordan_test(np.zeros((2,2)),-np.eye(2),np.zeros((2,2)))
    assert r['nullity']==2 and not r['order_two_supported']


def test_invertible_matrix_cannot_be_ep2():
    r=jordan_test(np.eye(2),np.array([[0,1],[0,0]]),np.zeros((2,2)))
    assert r['nullity']==0 and not r['chain_supported']


def test_semisimple_repeated_root_is_not_ep2():
    r=jordan_test(np.diag([0,0,2]),-np.eye(3),np.zeros((3,3)))
    assert not r['geometric_multiplicity_one']
    assert not r['order_two_supported']


def test_ep3_not_misreported_as_ep2():
    T=np.diag([1.,1.],1)
    r=jordan_test(T,-np.eye(3),np.zeros((3,3)))
    assert r['chain_supported']
    assert not r['order_two_supported']


def test_nonlinear_pencil_derivative_needed_for_exact_order():
    r=jordan_test([[0,1],[0,0]],-np.eye(2))
    assert r['chain_supported'] and not r['order_two_supported']


@pytest.mark.parametrize('kind',['constant','linear','sqrt','avoided'])
def test_puiseux_distinguishes_scaling(kind):
    t=np.geomspace(1e-8,1e-2,16)
    values={'constant':np.ones_like(t),'linear':3*t,'sqrt':2*np.sqrt(t),
            'avoided':np.sqrt(.1+t)}[kind]
    assert puiseux_test(t,values)['supported']==(kind=='sqrt')


def test_unresolved_puiseux_does_not_pass():
    t=np.geomspace(1e-8,1e-2,16)
    assert not puiseux_test(t,np.sqrt(t),uncertainties=np.sqrt(t))['supported']


@pytest.mark.parametrize('bad',[np.nan,np.inf,-1,0])
def test_invalid_puiseux_inputs_fail(bad):
    t=np.geomspace(1e-8,1e-2,16); y=np.sqrt(t); y[0]=bad
    with pytest.raises(ValueError):puiseux_test(t,y)


def test_missing_evidence_inconclusive_not_exclusion():
    r=evaluate_ep2(EP2Evidence([]))
    assert r['verdict']=='INCONCLUSIVE'


def test_shared_recurrence_a_b_not_independent():
    e=EP2Evidence(witnesses(lineage='recurrence'))
    r=evaluate_ep2(e)
    assert any('lineages' in x for x in r['reasons'])


def test_resolved_pair_is_only_pointwise():
    e=EP2Evidence(witnesses(1),partner_frequencies=witnesses(2))
    r=evaluate_ep2(e)
    assert r['verdict']=='DISTINCT_AT_TESTED_POINT'
    assert r['scope']=='TESTED_POINT_ONLY'


@pytest.mark.parametrize('name',sorted(WITHDRAWN))
def test_withdrawn_statistic_forbidden(name):
    with pytest.raises(ValueError):evaluate_ep2(EP2Evidence([],active_diagnostics=[name]))


def test_complete_synthetic_ep2_witness():
    t=np.geomspace(1e-8,1e-2,16)
    e=EP2Evidence(witnesses(),matrix=jordan_test([[0,1],[0,0]],-np.eye(2),np.zeros((2,2))),
                  puiseux=puiseux_test(t,2*np.sqrt(t)),first_loop=[1,0],second_loop=[0,1],
                  loops_converged=True,eigenvector_overlap=1.,eigenvector_trend_verified=True,
                  coalescence_located=True)
    r=evaluate_ep2(e)
    assert r['verdict']=='EP2_NUMERICALLY_SUPPORTED'
    assert not r['continuum_certificate']

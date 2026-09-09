import numpy as np
import pytest
from mp5d_science.research import PairAnchor,PairEvaluator,coordinate_route,parameters
from mp5d_science.convergence import Resolution


def anchor():return PairAnchor((.42,-.2,1.8),2,0,4,('N2','N3'),(2,3),(2.5-1.5j,2.6-2.05j))

def test_anchor_uses_decimal_parameter_arithmetic():
    p=parameters(anchor(),(.42,-.2,1.8))
    assert p.a=='.22' or p.a=='0.22'
    assert p.b=='0.62'

@pytest.mark.parametrize('step',[.1,.05,.025])
def test_coordinate_routes_use_actual_step(step):
    route=coordinate_route((0,0,0),(.2,.1,.3),step=step)
    assert route[0]==(0,0,0) and route[-1]==pytest.approx((.2,.1,.3))
    assert np.max(np.linalg.norm(np.diff(route,axis=0),axis=1))<=step+1e-14

def test_routes_change_geometry_not_just_names():
    a=coordinate_route((0,0,0),(.2,.1,.3),(0,1,2))
    b=coordinate_route((0,0,0),(.2,.1,.3),(2,1,0))
    assert a!=b

def test_unvalidated_anchor_cannot_start_continuation():
    p=PairEvaluator(anchor(),Resolution(precision_dps=30))
    with pytest.raises(ValueError,match='preflight'):p.follow([anchor().point,(.4,-.2,1.8)])

def test_high_precision_claim_requires_high_precision_configuration():
    with pytest.raises(ValueError,match='precision'):PairEvaluator(anchor(),Resolution())

def test_bad_route_order_is_not_silently_accepted():
    with pytest.raises(ValueError):coordinate_route((0,0,0),(0,0,1),(0,0,1))

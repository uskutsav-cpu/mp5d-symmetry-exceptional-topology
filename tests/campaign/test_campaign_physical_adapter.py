"""A real physical adapter smoke test; not a seven-axis science result."""
import math

from mp5d_science.convergence import Resolution
from mp5d_campaign.core import Store
from mp5d_campaign.numerics import Physics, z

from copy import deepcopy
from mp5d_campaign.core import AXES

def plan():
    return {'axes': deepcopy(AXES), 'higher_modes': {'overtones': [4,5], 'ell_offsets':[0,2,4]},
            'formal_domain_proposal': {'overtones':[0,1,2,3]},
            'adaptive_boundary': {'independent_continuation_directions':['forward','reverse']}}

def source(tmp_path):
    root=tmp_path/'src';root.mkdir();(root/'hello.py').write_text('x = 1\n');return root



def test_static_lambda_homotopy_matches_physical_static_root(tmp_path):
    root=source(tmp_path)
    config={'residual_tolerance':1e-8,'frequency_tolerance':1e-6}
    resolution=Resolution(cf_depth=160,angular_n=24,precision_dps=30,
                          radial_n=160,contour_length=60,scaling_angle_deg=70,continuation_step=.01)
    seed=complex(9.49117521848,-2.24647282923)/(2*math.pi)
    with Store(root,tmp_path/'run',config,plan(),development=True) as store:
        physics=Physics(store,config)
        transported=physics.static_lambda_root('8',seed,0,resolution)
        assert transported['status']=='PASS', transported
        physical=physics.solve('A',[0,0,0],{'m1':2,'m2':0,'ell':2},0,seed,resolution)
        assert physical['status']=='PASS', physical
        assert abs(z(transported['omega'])-z(physical['omega']))<1e-8

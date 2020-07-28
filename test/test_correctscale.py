import numpy as np
from deerlab import dipolarkernel, correctscale
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp

def test_correction():
# ======================================================================
    "Check that scale correction works properly"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = dipolarkernel(t,r,0.25,B)@P
    
    Vc,_ = correctscale(V*scale,t)

    assert max(abs(Vc-V)) < 1e-1
# ======================================================================

def test_model_deer():
# ======================================================================
    "Check scale correction using the deer model works"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = scale*dipolarkernel(t,r,0.25,B)@P
    
    _,V0 = correctscale(V,t,model='deer')

    assert abs(1 - V0/scale) < 0.1
# ======================================================================

def test_model_gauss():
# ======================================================================
    "Check scale correction using the gauss model works"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = scale*dipolarkernel(t,r,0.25,B)@P
    
    _,V0 = correctscale(V,t,model='gauss')

    assert abs(1 - V0/scale) < 0.15
# ======================================================================


def test_tmax():
# ======================================================================
    "Check scale correction using a manually specified tmax"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = dipolarkernel(t,r,0.25,B)@P
    
    Vc,_ = correctscale(V*scale,t,tmax=3.5)

    assert max(abs(Vc-V)) < 1e-1
# ======================================================================

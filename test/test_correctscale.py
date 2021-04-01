import numpy as np
from deerlab import dipolarkernel, correctscale
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp

def test_correction():
# ======================================================================
    "Check that scale correction using DEER model works"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = dipolarkernel(t,r,mod=0.25,bg=B)@P
    
    Vc = correctscale(V*scale,t,model='deer')

    assert max(abs(Vc-V)) < 1.5e-1
# ======================================================================

def test_model_gauss():
# ======================================================================
    "Check scale correction using the Gauss model works"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = dipolarkernel(t,r,mod=0.25,bg=B)@P
    
    Vc = correctscale(scale*V,t,model='gauss')

    assert max(abs(Vc-V)) < 1.5e-1
# ======================================================================


def test_tmax():
# ======================================================================
    "Check scale correction using a manually specified tmax"

    t = np.linspace(-1,5,300)
    r = np.linspace(3,7,200)
    P = dd_gauss(r,[5,0.2])
    B = bg_exp(t,0.3)
    scale = 1e8
    V = dipolarkernel(t,r,mod=0.25,bg=B)@P
    
    Vc = correctscale(V*scale,t,tmax=3.5)

    assert max(abs(Vc-V)) < 1.5e-1
# ======================================================================

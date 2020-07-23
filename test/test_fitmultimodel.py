
import numpy as np
from numpy import pi, inf, NaN
from deerlab import dipolarkernel, regoperator, regparamrange, selregparam, whitegaussnoise, fitmultimodel
from deerlab.dd_models import dd_gauss,dd_gauss3
from deerlab.utils import ovl

def test_multigauss():
#=======================================================================
    "Check that the fit of a multi-Gauss model works"
        
    r = np.linspace(2,6,300)
    t = np.linspace(-1,3,500)
    K = dipolarkernel(t,r)
    parin = [4, 0.2, 0.4, 4, 1, 0.4, 3, 0.4, 0.2]
    P = dd_gauss3(r,parin)
    V = K@P
    Pfit,_ = fitmultimodel(V,K,r,dd_gauss,5,'aicc',normP=False, uqanalysis=False)

    assert ovl(P,Pfit) > 0.99 # more than 90% overlap
#=======================================================================


import numpy as np
from numpy import pi, inf, NaN
from deerlab import dipolarkernel, regoperator, regparamrange, selregparam, whitegaussnoise, fitmultimodel
from deerlab.dd_models import dd_gengauss, dd_gauss, dd_rice, dd_rice3, dd_gauss3
from deerlab.utils import ovl

def test_multigauss():
#=======================================================================
    "Check that the fit of a multi-Gauss model works"
        
    r = np.linspace(2,6,300)
    t = np.linspace(0,6,500)
    K = dipolarkernel(t,r)
    parin = [4, 0.2, 0.4, 4, 1, 0.4, 3, 0.4, 0.2]
    P = dd_gauss3(r,parin)
    V = K@P

    Pfit,_ = fitmultimodel(V,K,r,dd_gauss,4,'aicc', uqanalysis=False)

    assert ovl(P,Pfit) > 0.99 # more than 90% overlap
#=======================================================================

def test_multirice():
#=======================================================================
    "Check that the fit of a multi-Rician model works"
        
    r = np.linspace(2,6,300)
    t = np.linspace(0,6,500)
    K = dipolarkernel(t,r)
    parin = [4, 0.2, 0.4, 4, 1, 0.4, 3, 0.4, 0.4]
    P = dd_rice3(r,parin)
    V = K@P

    Pfit,_ = fitmultimodel(V,K,r,dd_rice,4,'aicc', uqanalysis=False)

    assert ovl(P,Pfit) > 0.99 # more than 90% overlap
#=======================================================================


def test_multigengauss():
#=======================================================================
    "Check that the fit of a multi-generalized-Gaussian model works"
        
    r = np.linspace(2,6,300)
    t = np.linspace(0,6,500)
    K = dipolarkernel(t,r)
    P = dd_gengauss(r,[2.5, 0.5, 5]) + 0.8*dd_gengauss(r,[3, 0.7, 2])
    P /= np.trapz(P,r)
    V = K@P

    Pfit,_ = fitmultimodel(V,K,r,dd_gengauss,4,'aicc', uqanalysis=False)

    assert ovl(P,Pfit) > 0.99 # more than 90% overlap
#=======================================================================
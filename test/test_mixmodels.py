import numpy as np
from deerlab import mixmodels, dipolarkernel, fitparamodel
from deerlab.dd_models import dd_gauss, dd_rice

def test_gaussgauss():
# ======================================================================
    "Check the construction of a mixed model of Gaussian-Gaussian"

    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    parIn1 = [3, 0.5]
    P1 = dd_gauss(r,parIn1)
    parIn2 = [4, 0.5]
    P2 = dd_gauss(r,parIn2)
    P = 0.7*P2 + 0.3*P1

    mixedModel = mixmodels(dd_gauss,dd_gauss)
    parInMix = [3, 0.5, 0.3, 4, 0.5, 0.7]
    Pmix = mixedModel(r,parInMix)

    assert max(abs(Pmix - P)) < 1e-8
# ======================================================================


def test_gaussrice():
# ======================================================================
    "Check the construction of a mixed model of Gaussian-Gaussian"

    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    parIn1 = [3, 0.5]
    P1 = dd_gauss(r,parIn1)
    parIn2 = [4, 0.5]
    P2 = dd_rice(r,parIn2)
    P = 0.7*P2 + 0.3*P1

    mixedModel = mixmodels(dd_gauss,dd_rice)
    parInMix = [3, 0.5, 0.3, 4, 0.5, 0.7]
    Pmix = mixedModel(r,parInMix)

    assert max(abs(Pmix - P)) < 1e-8
# ======================================================================

def test_fit():
# ======================================================================
    "Check the fitting of a mixed model"
    t = np.linspace(-0.5,5,200)
    r = np.linspace(2,6,100)
    mixedModel = mixmodels(dd_gauss,dd_gauss)
    parInMix = [3, 0.5, 0.3, 4, 0.5, 0.7]
    Pmix = mixedModel(r,parInMix)

    K = dipolarkernel(t,r)
    V = K@Pmix
    Vmodel = lambda par: K@mixedModel(r,par)
    info = mixedModel()
    par0 = info['Start']
    lb = info['Lower'] 
    ub = info['Upper'] 
    parfit,_,_ = fitparamodel(V,Vmodel,par0,lb,ub)

    Pfit = mixedModel(r,parfit)

    assert max(abs(Pmix - Pfit)) < 1e-4
# ======================================================================

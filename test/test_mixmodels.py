import numpy as np
from deerlab import mixmodels
from deerlab.dd_models import dd_gauss, dd_rice

def test_gaussgauss():
# ======================================================================
    "Check the construction of a mixed model of Gaussian-Gaussian"

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


import numpy as np
from deerlab import dipolarkernel, whitegaussnoise, bootstrap_analysis, snlls
from deerlab.dd_models import dd_gauss
from deerlab.utils import assert_docstring

x = np.linspace(0,5,200)
sigma = 0.01
phaseshift = np.exp(1j*np.pi/2)
model = lambda par: par[0] + par[1]*x
model_global = lambda par: [model([par[0],par[1]]), model([par[2],par[3]])]
model_complex = lambda par: model(par)*phaseshift
yexp = model([1,3]) + whitegaussnoise(x,sigma)
yglobal1 = model([1,3]) + whitegaussnoise(x,sigma)
yglobal2 = model([2,4]) + whitegaussnoise(x,sigma)
yexp_complex = model_complex([1,3]) + whitegaussnoise(x,sigma)

def fitfcn(yexp):
    fit = snlls(yexp,model,[1,3],uq=False)
    return fit.nonlin*fit.lin

def fitfcn_global(ys):
    fit = snlls(ys,model_global,[1,3,2,4],uq=False)
    return fit.nonlin*fit.lin

def fitfcn_multiout(yexp):
    fit = snlls(yexp,model,[1,3],uq=False)
    return fit.nonlin*fit.lin, fit.model 

def fitfcn_complex(yexp):
    fit = snlls(yexp,model,[1,3],uq=False)
    return fit.nonlin*fit.lin


def test_basics():
# ======================================================================
    "Check the basic functionality of the bootstrapping"

    parfit = fitfcn(yexp)
    paruq = bootstrap_analysis(fitfcn,yexp,model(parfit),10)

    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

def test_resampling():
# ======================================================================
    "Check that both bootstrap resampling method yield similar results"

    parfit = fitfcn(yexp)
    paruq1 = bootstrap_analysis(fitfcn,yexp,model(parfit),10,resampling='residual')
    paruq2 = bootstrap_analysis(fitfcn,yexp,model(parfit),10,resampling='gaussian')

    assert all(abs(paruq1.mean - paruq2.mean) < 1.6e-2)
# ======================================================================


def test_multiple_ouputs():
# ======================================================================
    "Check that both bootstrap handles the correct number outputs"

    parfit,yfit = fitfcn_multiout(yexp)
    paruq = bootstrap_analysis(fitfcn_multiout,yexp,model(parfit),10)

    assert len(paruq)==2 and all(abs(paruq[0].mean - parfit)) and all(abs(paruq[1].mean - yfit))
# ======================================================================

def test_multiple_datasets():
# ======================================================================
    "Check bootstrapping when using multiple input datasets"

    parfit = fitfcn_global([yglobal1,yglobal2])
    paruq = bootstrap_analysis(fitfcn_global,[yglobal1,yglobal2],[model([1,3]),model([2,4])],5)

    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

def test_parallelization():
# ======================================================================
    "Check that bootstrap_analysis can run with multiple cores"

    parfit = fitfcn(yexp)
    paruq = bootstrap_analysis(fitfcn,yexp,model(parfit),10,cores=-1)

    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(bootstrap_analysis)
# ======================================================================

def test_complex_values():
# ======================================================================
    "Check the functionality of the bootstrapping with complex-valued outputs"

    parfit = fitfcn_complex(yexp_complex)
    paruq = bootstrap_analysis(fitfcn_complex,yexp_complex,model_complex(parfit),3)

    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

import pytest
# ======================================================================
def test_memory_limit():
    "Check that the memory limit works for too large analyses"
    # Request one million samples, approx. 80GB
    yfit = np.linspace(0,1,int(1e4))    
    yexp = yfit + whitegaussnoise(yfit,0.01)
    def bootfcn(y):
        return y
    with pytest.raises(MemoryError):
        paruq = bootstrap_analysis(bootfcn,yexp,yfit,1e6)
# ======================================================================

def test_noiselvl():
# ======================================================================
    "Check that noise levels can be specified"

    parfit = fitfcn(yexp)
    paruq = bootstrap_analysis(fitfcn,yexp,model(parfit),10, noiselvl=sigma)

    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

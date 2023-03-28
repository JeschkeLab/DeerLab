
import numpy as np
from deerlab import whitegaussnoise, bootstrap_analysis, snlls
from deerlab.utils import assert_docstring
import pytest 

# Fixtures
# --------------------------------------------------------------------------
x = np.linspace(0,5,200)
noiselvl = 0.01
phaseshift = np.exp(1j*np.pi/2)
paramset1 = [1,3]
paramset2 = [2,4]
 
def model(param):
    return param[0] + param[1]*x

def model_complex(param):
    return model(param)*phaseshift

def model_global(param):
    return [model([param[0],param[1]]), model([param[2],param[3]])]

@pytest.fixture(scope='module')
def dummy_data():
    return model(paramset1) + whitegaussnoise(x,noiselvl)

@pytest.fixture(scope='module')
def dummy_complex_data():
    return model_complex(paramset1) + whitegaussnoise(x,noiselvl)

@pytest.fixture(scope='module')
def dummy_global_data():
    data1 = model(paramset1) + whitegaussnoise(x,noiselvl)
    data2 = model(paramset2) + whitegaussnoise(x,noiselvl)
    return [data1,data2]

def fitfcn(yexp):
    fit = snlls(yexp,model,[1,3],uq=False)
    return fit.nonlin*fit.lin

def fitfcn_global(ys):
    fit = snlls(ys,model_global,[1,3,2,4],uq=False)
    return fit.nonlin*fit.lin

def fitfcn_multiout(yexp):
    fit = snlls(yexp,model,[1,3],uq=False)
    return fit.nonlin*fit.lin, fit.model, fit.nonlin[0]

def fitfcn_complex(yexp):
    fit = snlls(yexp,model,[1,3],uq=False)
    return fit.nonlin*fit.lin
# --------------------------------------------------------------------------


def test_basics(dummy_data):
# ======================================================================
    "Check the basic functionality of the bootstrapping"
    parfit = fitfcn(dummy_data)
    paruq = bootstrap_analysis(fitfcn,dummy_data,model(parfit),10)
    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

def test_resampling(dummy_data):
# ======================================================================
    "Check that both bootstrap resampling method yield similar results"
    parfit = fitfcn(dummy_data)
    paruq1 = bootstrap_analysis(fitfcn,dummy_data,model(parfit),10,resampling='residual')
    paruq2 = bootstrap_analysis(fitfcn,dummy_data,model(parfit),10,resampling='gaussian')
    assert all(abs(paruq1.mean - paruq2.mean) < 1.6e-2)
# ======================================================================

def test_multiple_ouputs(dummy_data):
# ======================================================================
    "Check that both bootstrap handles the correct number outputs"
    parfit,yfit,_ = fitfcn_multiout(dummy_data)
    paruq = bootstrap_analysis(fitfcn_multiout,dummy_data,model(parfit),10)
    assert len(paruq)==3 and all(abs(paruq[0].mean - parfit)) and all(abs(paruq[1].mean - yfit))
# ======================================================================

def test_multiple_datasets(dummy_global_data):
# ======================================================================
    "Check bootstrapping when using multiple input datasets"
    parfit = fitfcn_global(dummy_global_data)
    paruq = bootstrap_analysis(fitfcn_global,dummy_global_data,[model([1,3]),model([2,4])],5)
    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

def test_parallelization(dummy_data):
# ======================================================================
    "Check that bootstrap_analysis can run with multiple cores"
    parfit = fitfcn(dummy_data)
    paruq = bootstrap_analysis(fitfcn,dummy_data,model(parfit),10,cores=-1)
    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

def test_complex_values(dummy_complex_data):
# ======================================================================
    "Check the functionality of the bootstrapping with complex-valued outputs"
    parfit = fitfcn_complex(dummy_complex_data)
    paruq = bootstrap_analysis(fitfcn_complex,dummy_complex_data,model_complex(parfit),3)
    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

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

def test_noiselvl(dummy_data):
# ======================================================================
    "Check that noise levels can be specified"
    parfit = fitfcn(dummy_data)
    paruq = bootstrap_analysis(fitfcn,dummy_data,model(parfit),10, noiselvl=noiselvl)
    assert all(abs(paruq.mean - parfit) < 1.5e-2)
# ======================================================================

# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(bootstrap_analysis)
# ======================================================================

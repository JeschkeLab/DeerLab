import numpy as np
from deerlab import noiselevel, whitegaussnoise, dipolarkernel
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp
from deerlab.utils import assert_docstring
import pytest 

# Fixtures 
# ------------------------------------------------------------
t = np.linspace(0,3,200)
true_noiselvl = 0.03 
@pytest.fixture(scope='module')
def noise():
    return whitegaussnoise(t,true_noiselvl,seed=1,rescale=True)

@pytest.fixture(scope='module')
def mock_dataset(noise):
    r = np.linspace(2,6,100)
    P = dd_gauss(r,3, 0.5)
    lam = 0.25
    B = bg_exp(t,1.5)
    return dipolarkernel(t,r,mod=lam,bg=B)@P + noise

@pytest.fixture(scope='module')
def mock_multiple_datasets(noise, mock_dataset):
    N = 500
    dataset = mock_dataset - noise 
    datasets = np.zeros((len(t),N))
    for i in range(N):
        datasets[:,i] = dataset + whitegaussnoise(t,true_noiselvl,seed=i,rescale=True)
    return datasets

@pytest.fixture(scope='module')
def mock_dataset_complex(noise, mock_dataset):
    dataset = mock_dataset - noise 
    dataset = dataset*np.exp(-1j*2/5*np.pi);
    return dataset + noise + 1j*noise
# ------------------------------------------------------------


def test_filtered_movmean(mock_dataset):
#============================================================
    "Check estimation of noiselevel using a moving-mean filter"
    approxlevel = noiselevel(mock_dataset,'movmean',3)
    assert abs(approxlevel - true_noiselvl) < 1e-2
#============================================================

def test_der(mock_dataset):
#============================================================
    "Check estimation of noiselevel using a moving-mean filter"
    approxlevel = noiselevel(mock_dataset,'der')
    assert abs(approxlevel - true_noiselvl) < 1e-2
#============================================================

def test_reference(noise, mock_dataset):
#============================================================
    "Check estimation of noiselevel using a reference signal"
    approxlevel = noiselevel(mock_dataset,'reference',mock_dataset - noise)
    assert abs(approxlevel - true_noiselvl) < 1e-2
#============================================================

def test_filtered_savgol(mock_dataset):
#============================================================
    "Check estimation of noiselevel using a Savitzky-Golay filter"
    approxlevel = noiselevel(mock_dataset,'savgol',11,3)
    assert abs(approxlevel - true_noiselvl) < 1e-2
#============================================================

def test_multiscan(mock_multiple_datasets):
#============================================================
    "Check estimation of noiselevel using multiple scans of a signal"
    approxlevel = noiselevel(mock_multiple_datasets,'scans')
    assert abs(approxlevel - true_noiselvl) < 1e-2
#============================================================

def test_complex(mock_dataset_complex):
#============================================================
    "Check estimation of noiselevel using a complex signal"
    approxlevel = noiselevel(mock_dataset_complex,'complex')
    assert abs(approxlevel - true_noiselvl) < 1e-2
#============================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(noiselevel)
# ======================================================================

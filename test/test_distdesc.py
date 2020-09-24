
import numpy as np
from deerlab import dipolarkernel, distdesc, fitregmodel
from deerlab.dd_models import dd_gauss

# Gaussian distribution parameters
Pmean = 4 # nm
Psigma = 0.4 # nm

def assert_descriptor(key,truth):
# ----------------------------------------------------------------------
    # High-resolution distance axis
    r = np.linspace(2,6,2000)
    dr = np.mean(np.diff(r))
    # Gaussian distribution
    P = dd_gauss(r,[Pmean, Psigma])
    
    descriptor = distdesc(r,P)[0][key]
    
    assert abs(descriptor - truth) < dr
# ----------------------------------------------------------------------

def test_mean():
# ======================================================================
    "Check that the mean is correctly computed"

    truth = Pmean
    assert_descriptor('mean',truth)
# ======================================================================

def test_median():
# ======================================================================
    "Check that the median is correctly computed"

    truth = Pmean
    assert_descriptor('median',truth)
# ======================================================================

def test_mode():
# ======================================================================
    "Check that the mode is correctly computed"

    truth = Pmean
    assert_descriptor('mode',truth)
# ======================================================================

def test_iqm():
# ======================================================================
    "Check that the interquartile mean is correctly computed"

    truth = Pmean
    assert_descriptor('iqm',truth)
# ======================================================================

def test_iqr():
# ======================================================================
    "Check that the interquartile range is correctly computed"

    truth = 1.349*Psigma
    assert_descriptor('iqr',truth)
# ======================================================================

def test_std():
# ======================================================================
    "Check that the standard deviation is correctly computed"

    truth = Psigma
    assert_descriptor('std',truth)
# ======================================================================

def test_variance():
# ======================================================================
    "Check that the variance is correctly computed"

    truth = Psigma**2
    assert_descriptor('var',truth)
# ======================================================================

def test_mad():
# ======================================================================
    "Check that the mean absolute deviation is correctly computed"

    truth = np.sqrt(2/np.pi)*Psigma
    assert_descriptor('mad',truth)
# ======================================================================

def test_modality():
# ======================================================================
    "Check that the modality is correctly computed"

    truth = 1
    assert_descriptor('modality',truth)
# ======================================================================

def test_skewness():
# ======================================================================
    "Check that the skewness is correctly computed"

    truth = 0
    assert_descriptor('skewness',truth)
# ======================================================================

def test_kurtosis():
# ======================================================================
    "Check that the excess kurtosis is correctly computed"

    truth = 0
    assert_descriptor('kurtosis',truth)
# ======================================================================


t = np.linspace(0,5,200)
r = np.linspace(2,6,100)
P = dd_gauss(r,[Pmean, Psigma])
K = dipolarkernel(t,r)
fit = fitregmodel(K@P,K,r)


def assert_uncertainty(key):
# ----------------------------------------------------------------------
    desc,uq = distdesc(r,fit.P,fit.uncertainty)
    desc = desc[key]
    ci = uq[key].ci(95)
    assert (desc >= ci[0]) & (desc <= ci[1])
# ----------------------------------------------------------------------

def test_mean_uncertainty():
# ======================================================================
    "Check that the mean confidence intervals are correct."
    assert_uncertainty('mean')
# ======================================================================

def test_median_uncertainty():
# ======================================================================
    "Check that the median confidence intervals are correct."
    assert_uncertainty('median')
# ======================================================================

def test_iqm_uncertainty():
# ======================================================================
    "Check that the interquartile mean confidence intervals are correct."
    assert_uncertainty('iqm')
# ======================================================================

def test_std_uncertainty():
# ======================================================================
    "Check that the standard deviation confidence intervals are correct."
    assert_uncertainty('std')
# ======================================================================

def test_iqr_uncertainty():
# ======================================================================
    "Check that the interquartile range confidence intervals are correct."
    assert_uncertainty('iqr')
# ======================================================================

def test_mad_uncertainty():
# ======================================================================
    "Check that the mean absolute deviation confidence intervals are correct."
    assert_uncertainty('mad')
# ======================================================================

def test_var_uncertainty():
# ======================================================================
    "Check that the variance confidence intervals are correct."
    assert_uncertainty('var')
# ======================================================================

def test_skewness_uncertainty():
# ======================================================================
    "Check that the skewness confidence intervals are correct."
    assert_uncertainty('skewness')
# ======================================================================

def test_kurtosis_uncertainty():
# ======================================================================
    "Check that the kurtosis confidence intervals are correct."
    assert_uncertainty('kurtosis')
# ======================================================================
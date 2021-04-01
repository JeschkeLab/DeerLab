
import numpy as np
from deerlab import dipolarkernel, diststats, fitregmodel
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
    
    descriptor = diststats(r,P)[0][key]
    
    assert abs(descriptor - truth) < dr
# ----------------------------------------------------------------------

def test_rmin():
# ======================================================================
    "Check that the minimal distribution distance is correctly computed"

    truth = min(r)
    assert_descriptor('rmin',truth)
# ======================================================================

def test_rmax():
# ======================================================================
    "Check that the maximal distribution distance is correctly computed"

    truth = max(r)
    assert_descriptor('rmax',truth)
# ======================================================================

def test_int():
# ======================================================================
    "Check that the distribution integral is correctly computed"

    truth = 1
    assert_descriptor('int',truth)
# ======================================================================

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

def test_modes():
# ======================================================================
    "Check that the modes are correctly computed"

    truth = Pmean
    assert_descriptor('modes',truth)
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

def test_entropy():
# ======================================================================
    "Check that the entropy is correctly computed"

    truth = 1/2*np.log(2*np.pi*np.e*Psigma**2)
    assert_descriptor('entropy',truth)
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

def test_moment1():
# ======================================================================
    "Check that the 1st moment is correctly computed"

    truth = Pmean
    assert_descriptor('moment1',truth)
# ======================================================================

def test_moment2():
# ======================================================================
    "Check that the 2nd moment is correctly computed"

    truth = Psigma**2
    assert_descriptor('moment2',truth)
# ======================================================================

def test_moment3():
# ======================================================================
    "Check that the 3rd moment is correctly computed"

    truth = 0
    assert_descriptor('moment3',truth)
# ======================================================================

def test_moment4():
# ======================================================================
    "Check that the 4th moment is correctly computed"

    truth = 3
    assert_descriptor('moment4',truth)
# ======================================================================

t = np.linspace(0,5,200)
r = np.linspace(2,6,100)
P = dd_gauss(r,[Pmean, Psigma])
K = dipolarkernel(t,r)
fit = fitregmodel(K@P,K,r)


def assert_uncertainty(key):
# ----------------------------------------------------------------------
    desc,uq = diststats(r,fit.P,fit.Puncert)
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

def test_entropy_uncertainty():
# ======================================================================
    "Check that the entropy confidence intervals are correct."
    assert_uncertainty('entropy')
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

def test_moment1_uncertainty():
# ======================================================================
    "Check that the 1st moment confidence intervals are correct."
    assert_uncertainty('moment1')
# ======================================================================

def test_moment2_uncertainty():
# ======================================================================
    "Check that the 2nd moment confidence intervals are correct."
    assert_uncertainty('moment2')
# ======================================================================

def test_moment3_uncertainty():
# ======================================================================
    "Check that the 3rd moment confidence intervals are correct."
    assert_uncertainty('moment3')
# ======================================================================

def test_moment4_uncertainty():
# ======================================================================
    "Check that the 4th moment confidence intervals are correct."
    assert_uncertainty('moment4')
# ======================================================================


def assert_descriptor_nonuniform_r(key,truth):
# ----------------------------------------------------------------------
    # High-resolution non-uniform distance axis
    r = np.sqrt(np.linspace(2**2,6**2,2000))
    dr = np.min(np.diff(r))
    # Gaussian distribution
    P = dd_gauss(r,[Pmean, Psigma])
    
    descriptor = diststats(r,P)[0][key]
    assert abs(descriptor - truth) < dr
# ----------------------------------------------------------------------

def test_nonuniform_r_rmin():
# ======================================================================
    "Check that the minimal distribution distance is correctly computed"

    truth = min(r)
    assert_descriptor_nonuniform_r('rmin',truth)
# ======================================================================

def test_nonuniform_r_rmax():
# ======================================================================
    "Check that the maximal distribution distance is correctly computed"

    truth = max(r)
    assert_descriptor_nonuniform_r('rmax',truth)
# ======================================================================

def test_nonuniform_r_int():
# ======================================================================
    "Check that the distribution integral is correctly computed"

    truth = 1
    assert_descriptor_nonuniform_r('int',truth)
# ======================================================================

def test_nonuniform_r_mean():
# ======================================================================
    "Check that the mean is correctly computed"

    truth = Pmean
    assert_descriptor_nonuniform_r('mean',truth)
# ======================================================================

def test_nonuniform_r_median():
# ======================================================================
    "Check that the median is correctly computed"

    truth = Pmean
    assert_descriptor_nonuniform_r('median',truth)
# ======================================================================

def test_nonuniform_r_mode():
# ======================================================================
    "Check that the mode is correctly computed"

    truth = Pmean
    assert_descriptor_nonuniform_r('mode',truth)
# ======================================================================

def test_nonuniform_r_modes():
# ======================================================================
    "Check that the modes are correctly computed"

    truth = Pmean
    assert_descriptor_nonuniform_r('modes',truth)
# ======================================================================

def test_nonuniform_r_iqm():
# ======================================================================
    "Check that the interquartile mean is correctly computed"

    truth = Pmean
    assert_descriptor_nonuniform_r('iqm',truth)
# ======================================================================

def test_nonuniform_r_iqr():
# ======================================================================
    "Check that the interquartile range is correctly computed"

    truth = 1.349*Psigma
    assert_descriptor_nonuniform_r('iqr',truth)
# ======================================================================

def test_nonuniform_r_std():
# ======================================================================
    "Check that the standard deviation is correctly computed"

    truth = Psigma
    assert_descriptor_nonuniform_r('std',truth)
# ======================================================================

def test_nonuniform_r_variance():
# ======================================================================
    "Check that the variance is correctly computed"

    truth = Psigma**2
    assert_descriptor_nonuniform_r('var',truth)
# ======================================================================

def test_nonuniform_r_entropy():
# ======================================================================
    "Check that the entropy is correctly computed"

    truth = 1/2*np.log(2*np.pi*np.e*Psigma**2)
    assert_descriptor_nonuniform_r('entropy',truth)
# ======================================================================

def test_nonuniform_r_mad():
# ======================================================================
    "Check that the mean absolute deviation is correctly computed"

    truth = np.sqrt(2/np.pi)*Psigma
    assert_descriptor_nonuniform_r('mad',truth)
# ======================================================================

def test_nonuniform_r_modality():
# ======================================================================
    "Check that the modality is correctly computed"

    truth = 1
    assert_descriptor_nonuniform_r('modality',truth)
# ======================================================================

def test_nonuniform_r_skewness():
# ======================================================================
    "Check that the skewness is correctly computed"

    truth = 0
    assert_descriptor_nonuniform_r('skewness',truth)
# ======================================================================

def test_nonuniform_r_kurtosis():
# ======================================================================
    "Check that the excess kurtosis is correctly computed"

    truth = 0
    assert_descriptor_nonuniform_r('kurtosis',truth)
# ======================================================================

def test_nonuniform_r_moment1():
# ======================================================================
    "Check that the 1st moment is correctly computed"

    truth = Pmean
    assert_descriptor_nonuniform_r('moment1',truth)
# ======================================================================

def test_nonuniform_r_moment2():
# ======================================================================
    "Check that the 2nd moment is correctly computed"

    truth = Psigma**2
    assert_descriptor_nonuniform_r('moment2',truth)
# ======================================================================

def test_nonuniform_r_moment3():
# ======================================================================
    "Check that the 3rd moment is correctly computed"

    truth = 0
    assert_descriptor_nonuniform_r('moment3',truth)
# ======================================================================

def test_nonuniform_r_moment4():
# ======================================================================
    "Check that the 4th moment is correctly computed"

    truth = 3
    assert_descriptor_nonuniform_r('moment4',truth)
# ======================================================================
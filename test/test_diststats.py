
from deerlab.utils import assert_docstring
import numpy as np
from deerlab import dipolarkernel, diststats, snlls
from deerlab.dd_models import dd_gauss
from deerlab.utils import assert_docstring
from pytest import fixture,mark

# Fixtures
# -----------------------------------------------------------------------
Pmean = 4 # nm
Psigma = 0.4 # nm
integral = 1e5
rmin,rmax,dr = 2,6,0.002

@fixture(scope='module')
def statistical_descriptors():
    r = np.arange(rmin,rmax+dr,dr)
    P = integral*dd_gauss(r,Pmean, Psigma)    
    return diststats(r,P)[0]    

@fixture(scope='module')
def statistical_descriptors_nonuniform_r():
    r = np.sqrt(np.arange(rmin**2,(rmax)**2+dr,dr))
    P = integral*dd_gauss(r,Pmean,Psigma)    
    return diststats(r,P)[0]

@fixture(scope='module')
def fitted_statistical_descriptors():
    t = np.linspace(0,5,200)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,Pmean,Psigma)
    K = dipolarkernel(t,r)
    fit = snlls(K@P,K,lbl=np.zeros_like(r))
    return diststats(r,fit.param,fit.paramUncert)

truths = [
    ('rmin', rmin),
    ('rmax', rmax),
    ('int', integral),
    ('mean', Pmean),
    ('median', Pmean),
    ('mode', Pmean),
    ('modes',Pmean),
    ('iqm', Pmean),
    ('iqr', 1.349*Psigma),
    ('std', Psigma),
    ('var', Psigma**2),
    ('entropy', 1/2*np.log(2*np.pi*np.e*Psigma**2)),
    ('mad', np.sqrt(2/np.pi)*Psigma),
    ('modality', 1),
    ('skewness', 0),
    ('kurtosis', 0),
    ('moment1', Pmean),
    ('moment2', Psigma**2),
    ('moment3', 0),
    ('moment4', 3)
]
# ----------------------------------------------------------------------

# ======================================================================
@mark.parametrize('key, truth',truths)
def test_statistical_descriptors_values(statistical_descriptors,key,truth):
    "Check that the statistical descriptors are correctly computes"
    assert abs(statistical_descriptors[key] - truth) < dr
# ======================================================================

# ======================================================================
@mark.parametrize('key',['mean','median','iqm','iqr','std','var','entropy','mad',
                         'skewness','kurtosis','moment1','moment2','moment3','moment4',])
def test_statistical_descriptors_uncertainties(fitted_statistical_descriptors,key):
    values,uq = fitted_statistical_descriptors
    value = values[key]
    ci = uq[key].ci(95)
    assert (value >= ci[0]) & (value <= ci[1])
# ======================================================================

# ======================================================================
@mark.parametrize('key, truth',truths)
def test_statistical_descriptors_values_nonuniform_r(statistical_descriptors_nonuniform_r,key,truth):
    "Check that the statistical descriptors are correctly computes"
    assert abs(statistical_descriptors_nonuniform_r[key] - truth) < dr
# ======================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(diststats)
# ======================================================================

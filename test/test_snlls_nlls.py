
import numpy as np
from deerlab import dipolarkernel, whitegaussnoise, snlls, regoperator
from deerlab.dd_models import dd_gauss
from deerlab.utils import skip_on
from scipy.linalg import block_diag
import pytest 

# Fixtures
# -----------------------------------------------------------------------
x = np.linspace(0,3,300)
y = np.linspace(2,5,200)
par = [3,0.2]
par0 = [5, 0.5]
lb = [1, 0.1]
ub = [20, 1]
noiselvl = 0.03 
 
@pytest.fixture(scope='module')
def design_matrix():
    return dipolarkernel(x,y)

@pytest.fixture(scope='module')
def mock_model(design_matrix):
    return lambda p: design_matrix@dd_gauss(y,*p)

@pytest.fixture(scope='module')
def mock_data(mock_model):
    return mock_model(par)
# -----------------------------------------------------------------------


def test_NLLS_fit(mock_model, mock_data):
# ======================================================================
    "Check the NLLS fit"
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    assert all(abs(par - fit.nonlin) < 1e-2)
# ======================================================================

def test_NLLS_fit_unconstrained(mock_model, mock_data):
# ======================================================================
    "Check the fit of a model with unconstrained parameters"
    fit = snlls(mock_data,mock_model,par0)
    assert all(abs(par - fit.nonlin) < 1e-2)
# =======================================================================

def test_NLLS_fit_confidence_intervals(mock_model, mock_data):
# ======================================================================
    "Check that the confidence intervals of the fitted parameters are correct"
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    paruq, pfit = fit.nonlinUncert, fit.nonlin
    p95lb,p95ub = paruq.ci(95).T
    p50lb,p50ub = paruq.ci(50).T
    errors = []
    if not np.all(p95lb <= pfit) and not np.all(p50lb <= pfit):
        errors.append("Some fitted values are below the lower bound of the confidence intervals.")
    if not np.all(p95ub >= pfit) and not np.all(p50lb >= pfit):
        errors.append("Some fitted values are over the upper bound of the confidence intervals.")
    if not np.all(p95lb <= p50lb):
        errors.append("The 50%-CI has lower values than the 95%-CI")
    if not np.all(p95ub >= p50ub):
        errors.append("The 50%-CI has larger values than the 95%-CI")
    if not np.all(np.minimum(lb,p95lb)==lb):
        errors.append("The lower bounds are not satisfied by the confidence intervals.")
    if not np.all(np.maximum(ub,p95ub)==ub):
        errors.append("The upper bounds are not satisfied by the confidence intervals.")
    assert not errors, f"Errors occured:\n{chr(10).join(errors)}"
# ======================================================================

# ======================================================================
@skip_on('_tkinter.TclError', reason="A problem with the Tk backend occured")
def test_NLLS_fit_plot(mock_model, mock_data):
    "Check that the plot method works"
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    fig = fit.plot()
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ======================================================================

# ======================================================================
def test_NLLS_fit_cost_value(mock_model, mock_data):
    "Check that the cost value is properly returned"
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    assert isinstance(fit.cost,float) and np.round(fit.cost/(np.sum(fit.residuals**2)/len(fit.residuals)),5)==1
# ======================================================================

# ======================================================================
def test_NLLS_goodness_of_fit(mock_model, mock_data):
    "Check the goodness-of-fit statistics are correct even with arbitrary scaling"
    noisy_data =  mock_data+ whitegaussnoise(x,noiselvl,seed=1,rescale=True)
    fit = snlls(noisy_data, mock_model, par0, lb=[1, 0.1], ub=[20, 1], noiselvl=noiselvl)
    assert abs(fit.stats['chi2red'] - 1) < 0.05
# ======================================================================

# ======================================================================
def test_NLLS_goodness_of_fit_when_scaled(mock_model, mock_data):
    "Check the goodness-of-fit statistics are correct even with arbitrary scaling"
    scale = 1e5
    noisy_data =  scale*(mock_data+ whitegaussnoise(x,noiselvl,seed=1,rescale=True))
    fit = snlls(noisy_data, mock_model, par0, lb=[1, 0.1], ub=[20, 1], noiselvl=scale*noiselvl)
    assert abs(fit.stats['chi2red'] - 1) < 0.05
# ======================================================================

# ======================================================================
def test_NLLS_global_fit(mock_model, mock_data):
    "Check that global fitting yields correct results"
    fit = snlls([mock_data]*2, lambda p: [mock_model(p)]*2 ,par0,lb,ub)
    assert all(abs(par - fit.nonlin) < 1e-2)
# ======================================================================

# ======================================================================
def test_NLLS_global_fit_with_weights(mock_model, mock_data):
    "Check that the global weights properly work when specified"
    fit1 = snlls([mock_data]*2, lambda p: [mock_model(p)]*2 ,par0,lb,ub, weights=[1,1e-10])
    fit2 = snlls([mock_data]*2, lambda p: [mock_model(p)]*2 ,par0,lb,ub, weights=[1e-10,1])
    assert all(abs(fit1.nonlin/par-1) < 0.03) and all(abs(fit2.nonlin/par-1) < 0.03)
# ======================================================================

# ======================================================================
def test_NLLS_fit_with_extrapenalty(mock_model, mock_data):
    "Check that custom penalties can be passed and act on the solution"
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    tikhonov = lambda p,_: 1e-4*regoperator(y,2)@dd_gauss(y,*p)
    fit_tikh = snlls(mock_data,mock_model,par0,lb,ub, extrapenalty=tikhonov)
    assert  np.linalg.norm(fit.nonlin - par) < np.linalg.norm(fit_tikh.nonlin - par)
# ======================================================================

# ======================================================================
def test_NLLS_fit_with_frozen_param(mock_model, mock_data):
    "Check that parameters can be frozen during the optimization"
    frozenpars = [3,None]
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    fit_frozen = snlls(mock_data,mock_model,par0,lb,ub, nonlin_frozen=frozenpars)
    assert np.allclose(fit_frozen.model,fit.model)
# ======================================================================

# ======================================================================
def test_NLLS_fit_Nparam_with_frozen_param(mock_model, mock_data):
    "Check that the correct number of parameters are returned even with frozen parameters"
    frozenpars = [3,None]
    fit = snlls(mock_data,mock_model,par0,lb,ub)
    fit_frozen = snlls(mock_data,mock_model,par0,lb,ub, nonlin_frozen=frozenpars)
    assert len(fit.nonlin)==2 and len(fit_frozen.nonlin)==2
# ======================================================================
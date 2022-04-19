
import numpy as np
from deerlab import dipolarkernel, whitegaussnoise, snlls, regoperator
from deerlab.dd_models import dd_gauss, dd_gauss2, dd_rice
from deerlab.utils import skip_on, ovl, assert_docstring
from scipy.linalg import block_diag

def test_gaussian():
# ======================================================================
    "Check the fit of a dipolar evolution function originating from a Gaussian"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,*par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V,model,par0,lb,ub)

    assert all(abs(par - fit.nonlin) < 1e-2)
# ======================================================================

def test_rice():
# ======================================================================
    "Check the fit of a dipolar evolution function originating from a random coil"
        
    t = np.linspace(0,3,400)
    r = np.linspace(1,6,150)
    par = np.array([4, 0.17])
    P = dd_rice(r,*par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [3,0.2]
    lb = [1,0.05]
    ub = [20,5]
    model = lambda p: K@dd_rice(r,*p)
    fit = snlls(V,model,par0,lb,ub)

    assert all(abs(par - fit.nonlin) < 1e-2)
# ======================================================================


def test_unconstrained():
# ======================================================================
    "Check the fit of a model with unconstrained parameters"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,*par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V,model,par0)

    assert all(abs(par - fit.nonlin) < 1e-2)
# =======================================================================

 
def assert_confidence_intervals(pci50,pci95,pfit,lb,ub):
# ----------------------------------------------------------------------

    p95lb = pci95[:,0]
    p95ub = pci95[:,1]
    p50lb = pci50[:,0]
    p50ub = pci50[:,1]
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
#----------------------------------------------------------------------


def test_confinter_Pfit():
# ======================================================================
    "Check that the confidence intervals of the fitted parameters are correct"
    
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    P = dd_gauss(r,3,0.2)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V,model,par0,lb,ub)
    paruq = fit.nonlinUncert
    parfit = fit.nonlin

    assert_confidence_intervals(paruq.ci(50),paruq.ci(95),parfit,lb,ub)
# ======================================================================

def test_globalfit():
# ======================================================================
    "Check that global fitting yields correct results"
        
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,*par)

    t1 = np.linspace(0,3,300)
    K1 = dipolarkernel(t1,r)
    V1 = K1@P

    t2 = np.linspace(-0.5,4,200)
    K2 = dipolarkernel(t2,r)
    V2 = K2@P

    def Vmodel(par):
        Pfit = dd_gauss(r,*par)
        V1 = K1@Pfit
        V2 = K2@Pfit
        return [V1,V2]
    
    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    fit = snlls([V1,V2],Vmodel,par0,lb,ub)

    assert all(abs(par - fit.nonlin) < 1e-2)
# ======================================================================

def test_globalfit_scales():
#============================================================
    "Check that the global fit with arbitrary amplitudes works."
    t1 = np.linspace(0,5,300)
    t2 = np.linspace(0,2,300)
    r = np.linspace(3,5,100)
    P = dd_gauss(r,4,0.25)
    K1 = dipolarkernel(t1,r) 
    K2 = dipolarkernel(t2,r)
    scales = [1e3, 1e9]
    V1 = scales[0]*K1@P
    V2 = scales[1]*K2@P

    def Vmodel(par):
        Pfit = dd_gauss(r,*par)
        V1 = K1@Pfit
        V2 = K2@Pfit
        return block_diag(V1,V2).T

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    fit = snlls(np.concatenate([V1,V2]),Vmodel,par0,lb,ub)

    assert max(abs(np.asarray(scales)/fit.lin - 1)) < 1e-2 
#============================================================

@skip_on('_tkinter.TclError', reason="A problem with the Tk backend occured")
def test_plot():
# ======================================================================
    "Check that the plot method works"
  
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    P = dd_gauss(r,3,0.2)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V,model,par0,lb,ub)
    
    fig = fit.plot()
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ======================================================================

def test_cost_value():
# ======================================================================
    "Check that the cost value is properly returned"

    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    P = dd_gauss(r,3,0.2)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V,model,par0,lb,ub)
    assert isinstance(fit.cost,float) and np.round(fit.cost/np.sum(fit.residuals**2),5)==1
# ======================================================================

def test_global_weights():
# ======================================================================
    "Check that the global weights properly work when specified"

    t = np.linspace(0,5,300)
    r = np.linspace(2,8,600)
    K = dipolarkernel(t,r)

    param1 = [3,0.2]
    param2 = [5,0.2]
    P1 = dd_gauss(r,*param1)
    P2 = dd_gauss(r,*param2)
    V1 = K@P1 + whitegaussnoise(t,0.01,seed=1)
    V2 = K@P2 + whitegaussnoise(t,0.01,seed=1)

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: [K@dd_gauss(r,*p)]*2
    fit1 = snlls([V1,V2],model,par0,lb,ub,weights=[1,1e-10])
    fit2 = snlls([V1,V2],model,par0,lb,ub,weights=[1e-10,1])

    assert all(abs(fit1.nonlin/param1-1) < 0.03) and all(abs(fit2.nonlin/param2-1) < 0.03)
# ======================================================================

def test_global_weights_default():
# ======================================================================
    "Check the correct fit of two signals when one is of very low quality"

    t = np.linspace(0,5,300)
    r = np.linspace(2,6,90)
    param = [4.5, 0.25]
    P = dd_gauss(r,*param)

    K = dipolarkernel(t,r)
    scales = [1e3,1e9]
    V1 = scales[0]*K@P + whitegaussnoise(t,0.001,seed=1)
    V2 = scales[1]*K@P + whitegaussnoise(t,0.1,seed=1)
    
    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: [K@dd_gauss(r,*p)]*2
    fit = snlls([V1,V2],model,par0,lb,ub,weights=[1,0])

    assert all(abs(fit.nonlin/param-1) < 0.03)
# ======================================================================

def test_goodness_of_fit():
# ======================================================================
    "Check the goodness-of-fit statistics are correct even with arbitrary scaling"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    P = dd_gauss(r,3,0.2)
    K = dipolarkernel(t,r)
    sigma = 0.03
    V = K@P + whitegaussnoise(t,sigma,seed=1,rescale=True)

    par0 = [5, 0.5]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V, model, par0, lb=[1, 0.1], ub=[20, 1], noiselvl=sigma)
    
    assert abs(fit.stats['chi2red'] - 1) < 0.05
# ======================================================================

def test_goodness_of_fit_scaled():
# ======================================================================
    "Check the goodness-of-fit statistics are correct even with arbitrary scaling"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    P = dd_gauss(r,3,0.2)
    K = dipolarkernel(t,r)
    V0 = 1e5
    sigma = V0*0.03
    V = V0*K@P + whitegaussnoise(t,sigma,seed=1,rescale=True)

    par0 = [5, 0.5]
    model = lambda p: K@dd_gauss(r,*p)
    fit = snlls(V, model, par0, lb=[1, 0.1], ub=[20, 1], noiselvl=sigma)
    
    assert abs(fit.stats['chi2red'] - 1) < 0.05
# ======================================================================

def test_extrapenalty():
# ======================================================================
    "Check that custom penalties can be passed and act on the solution"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    P = dd_gauss2(r,3.5,0.5,0.5,4,0.1,0.5)
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.15,seed=1)

    par0 = [2.5, 0.01, 0.1, 4.5, 0.01, 0.6]
    lb = [1, 0.01, 0, 1, 0.01, 0]
    ub = [20, 1, 1, 20, 1, 1]
    # Fit case it fails, stuck at "spicky" Gaussians
    model = lambda p: K@dd_gauss2(r,*p)
    fit = snlls(V,model,par0,lb,ub)

    # Fit with Tikhonov penalty on the Gaussians model
    L = regoperator(r,2)
    alpha = 1e-4
    tikhonov = lambda p,_: alpha*L@dd_gauss2(r,*p)
    fit_tikh = snlls(V,model,par0,lb,ub, extrapenalty=tikhonov)

    Pfit = dd_gauss2(r,*fit.nonlin)
    Pfit_tikh = dd_gauss2(r,*fit_tikh.nonlin)

    assert  ovl(P,Pfit)<ovl(P,Pfit_tikh)
# ======================================================================

def test_frozen_param():
# ======================================================================
    "Check that parameters can be frozen during the optimization"
    r = np.linspace(0,6,300)
    def model(param):
        mean1,mean2,std1,std2,amp1,amp2 = param
        return amp1*dd_gauss(r,mean1,std1) + amp2*dd_gauss(r,mean2,std2)

    y = model([3,4,0.2,0.3,0.5,0.6])
    par0 = [3.1,4.2,0.2,0.3,0.5,0.6]
    lb = [0,0,0.01,0.01,0,0]
    ub = [10,10,5,5,1,1]

    frozenpars = [None,4,None,0.3,None,None]
    fit = snlls(y,model,par0,lb,ub)
    fit_frozen = snlls(y,model,par0,lb,ub, nonlin_frozen=frozenpars)

    assert np.allclose(fit_frozen.model,fit.model) and np.allclose(fit_frozen.model,y)
# ======================================================================

def test_frozen_Nparam():
# ======================================================================
    "Check that the correct number of parameters are returned even with frozen parameters"
    r = np.linspace(0,6,300)
    def model(param):
        mean1,mean2,std1,std2,amp1,amp2 = param
        return amp1*dd_gauss(r,mean1,std1) + amp2*dd_gauss(r,mean2,std2)
    y = model([3,4,0.2,0.3,0.5,0.6])
    par0 = [2,2,0.5,0.5,0.5,0.5]
    lb = [0,0,0.01,0.01,0,0]
    ub = [10,10,5,5,1,1]

    frozenpars = [None,4,None,0.3,None,None]
    fit = snlls(y,model,par0,lb,ub)
    fit_frozen = snlls(y,model,par0,lb,ub, nonlin_frozen=frozenpars)

    assert len(fit.nonlin)==6 and len(fit_frozen.nonlin)==6
# ======================================================================
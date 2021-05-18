
import numpy as np
from deerlab import dipolarkernel, whitegaussnoise, fitparamodel
from deerlab.dd_models import dd_gauss, dd_rice


def test_gaussian():
# ======================================================================
    "Check the fit of a dipolar evolution function originating from a Gaussian"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,p)
    fit = fitparamodel(V,model,par0,lb,ub)

    assert all(abs(par - fit.param) < 1e-2)
# ======================================================================


def test_rice():
# ======================================================================
    "Check the fit of a dipolar evolution function originating from a random coil"
        
    t = np.linspace(0,3,400)
    r = np.linspace(1,6,150)
    par = np.array([4, 0.17])
    P = dd_rice(r,par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = dd_rice.start
    lb = dd_rice.lower
    ub = dd_rice.upper
    model = lambda p: K@dd_rice(r,p)
    fit = fitparamodel(V,model,par0,lb,ub)

    assert all(abs(par - fit.param) < 1e-2)
# ======================================================================


def test_unconstrained():
# ======================================================================
    "Check the fit of a model with unconstrained parameters"
        
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    model = lambda p: K@dd_gauss(r,p)
    fit = fitparamodel(V,model,par0)

    assert all(abs(par - fit.param) < 1e-2)
# =======================================================================


def test_rescaling():
# =======================================================================
    "Check that rescaling does not alter the results"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[4, 0.2])
    K = dipolarkernel(t,r)

    scale = 1e9
    V = K@P
    par0 = [3, 0.3]
    lb = [1, 0.1]
    ub = [20, 5]
    mymodel = lambda param: dipolarkernel(t,r)@dd_gauss(r,param)

    fit1 = fitparamodel(V*scale,mymodel,par0,lb,ub,fitscale=True)
    fit2 = fitparamodel(V,mymodel,par0,lb,ub,fitscale=False)

    assert all(abs(fit1.param - fit2.param) < 1e-2)
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
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,p)
    fit = fitparamodel(V,model,par0,lb,ub)
    paruq = fit.paramUncert
    parfit = fit.param

    assert_confidence_intervals(paruq.ci(50),paruq.ci(95),parfit,lb,ub)
# ======================================================================

def test_manual_covmatrix():
# ======================================================================
    "Check that covariance matrix can be manually specified"
    
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    K = dipolarkernel(t,r)
    sig = 0.01
    V = K@P + whitegaussnoise(t,sig,seed=1)
    covmat = sig**2*np.eye(len(t))

    par0 = [5, 0.50]
    lb   = [2, 0.05]
    ub   = [5, 0.70]
    model = lambda p: K@dd_gauss(r,p)
    fitmanual = fitparamodel(V,model,par0,lb,ub, covmatrix = covmat)
    fitauto = fitparamodel(V,model,par0,lb,ub)

    paruq_manual = fitmanual.paramUncert
    paruq_auto = fitauto.paramUncert

    ci_manual = paruq_manual.ci(95)
    ci_auto = paruq_auto.ci(95)

    assert np.all(abs(ci_manual - ci_auto) < 1e-2)
# ======================================================================

def test_globalfit():
# ======================================================================
    "Check that global fitting yields correct results"
        
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    t1 = np.linspace(0,3,300)
    K1 = dipolarkernel(t1,r)
    V1 = K1@P

    t2 = np.linspace(-0.5,4,200)
    K2 = dipolarkernel(t2,r)
    V2 = K2@P

    def Vmodel(par):
        Pfit = dd_gauss(r,par)
        V1 = K1@Pfit
        V2 = K2@Pfit
        return [V1,V2]
    
    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    fit = fitparamodel([V1,V2],Vmodel,par0,lb,ub)

    assert all(abs(par - fit.param) < 1e-2)
# ======================================================================

def test_globalfit_scales():
#============================================================
    "Check that the global fit with arbitrary amplitudes works."
    t1 = np.linspace(0,5,300)
    t2 = np.linspace(0,2,300)
    r = np.linspace(3,5,100)
    P = dd_gauss(r,[4,0.25])
    K1 = dipolarkernel(t1,r) 
    K2 = dipolarkernel(t2,r)
    scales = [1e3, 1e9]
    V1 = scales[0]*K1@P
    V2 = scales[1]*K2@P

    def Vmodel(par):
        Pfit = dd_gauss(r,par)
        V1 = K1@Pfit
        V2 = K2@Pfit
        return [V1,V2]

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    fit = fitparamodel([V1,V2],Vmodel,par0,lb,ub,weights=[1,1])

    assert max(abs(np.asarray(scales)/np.asarray(fit.scale) - 1)) < 1e-2 
#============================================================

def test_plot():
# ======================================================================
    "Check that the plot method works"
  
    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,p)
    fit = fitparamodel(V,model,par0,lb,ub)
    
    fig = fit.plot(show=False)
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ======================================================================

def test_cost_value():
# ======================================================================
    "Check that the cost value is properly returned"

    t = np.linspace(0,3,300)
    r = np.linspace(2,5,200)
    par = np.array([3,0.2])
    P = dd_gauss(r,par)

    K = dipolarkernel(t,r)
    V = K@P

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: K@dd_gauss(r,p)
    fit = fitparamodel(V,model,par0,lb,ub)
    assert isinstance(fit.cost,float) and np.round(fit.cost/np.sum(fit.residuals**2),5)==1
# ======================================================================


def test_confinter_values():
# ======================================================================
    "Check that the values of the confidence intervals are correct"

    x = np.linspace(0.3, 10, 100)
    np.random.seed(0)
    p = [0.1,2]
    y = p[0]*x + p[1] + 0.5*np.random.randn(x.size)

    # Reference confidence intervals (calculated using lmfit package, see #116)
    cov = [99.73,95.45,68.27]
    a_ci_ref = [[0.02947257081676860, 0.13987679478796447],
                [0.04834630918638558, 0.12100143673822517], 
                [0.06664985961745347, 0.10269723504693312]]
    b_ci_ref = [[1.7844876967028012, 2.433179427358454], 
                [1.8953902386010164, 2.322276885466784], 
                [2.0029218541761944, 2.214745269892137]]
    
    fit = fitparamodel(y,lambda p: p[0]*x + p[1],[0.1,1],fitscale=False)
    a_ci = [fit.paramUncert.ci(cov[i])[0,:] for i in range(3)]
    b_ci = [fit.paramUncert.ci(cov[i])[1,:] for i in range(3)]

    ci_match = lambda ci,ci_ref,truth:np.max(abs(np.array(ci) - np.array(ci_ref)))/truth < 0.01
    assert ci_match(a_ci,a_ci_ref,p[0]) & ci_match(b_ci,b_ci_ref,p[1])
# ======================================================================

def test_global_weights():
# ======================================================================
    "Check that the global weights properly work when specified"

    t = np.linspace(0,5,300)
    r = np.linspace(2,8,600)
    K = dipolarkernel(t,r)

    param1 = [3,0.2]
    param2 = [5,0.2]
    P1 = dd_gauss(r,param1)
    P2 = dd_gauss(r,param2)
    V1 = K@P1 + whitegaussnoise(t,0.01,seed=1)
    V2 = K@P2 + whitegaussnoise(t,0.01,seed=1)

    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: [K@dd_gauss(r,p)]*2
    fit1 = fitparamodel([V1,V2],model,par0,lb,ub,weights=[1,0])
    fit2 = fitparamodel([V1,V2],model,par0,lb,ub,weights=[0,1])

    assert all(abs(fit1.param/param1-1) < 0.03) and all(abs(fit2.param/param2-1) < 0.03)
# ======================================================================

def test_global_weights_default():
# ======================================================================
    "Check the correct fit of two signals when one is of very low quality"

    t = np.linspace(0,5,300)
    r = np.linspace(2,6,90)
    param = [4.5, 0.25]
    P = dd_gauss(r,param)

    K = dipolarkernel(t,r)
    scales = [1e3,1e9]
    V1 = scales[0]*K@P + whitegaussnoise(t,0.001,seed=1)
    V2 = scales[1]*K@P + whitegaussnoise(t,0.1,seed=1)
    
    par0 = [5, 0.5]
    lb = [1, 0.1]
    ub = [20, 1]
    model = lambda p: [K@dd_gauss(r,p)]*2
    fit = fitparamodel([V1,V2],model,par0,lb,ub,weights=[1,0])

    assert all(abs(fit.param/param-1) < 0.03)
# ======================================================================
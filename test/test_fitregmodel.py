
import numpy as np
from numpy import pi, inf, NaN
from deerlab import fitregmodel,dipolarkernel, regoperator, regparamrange, selregparam, whitegaussnoise
from deerlab.dd_models import dd_gauss,dd_gauss2

def assert_solver(solver):
#============================================================

    np.random.seed(1)
    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.5])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)
    Pfit,_ = fitregmodel(V,K,r,'tikhonov','aic',solver=solver)

    assert np.max(abs(Pfit - P)) < 5.5e-2
#============================================================

def test_solver_fnnls():
#============================================================
    "Check that the 'fnnls' solver works"

    assert_solver('fnnls')
#============================================================

def test_solver_nnlsbpp():
#============================================================
    "Check that the 'nnlsbpp' solver works"

    assert_solver('nnlsbpp')
#============================================================

def test_solver_cvx():
#============================================================
    "Check that the 'cvx' solver works"

    assert_solver('cvx')
#============================================================

def test_unconstrained():
#============================================================
    "Check that unconstrained distributions are correctly fitted"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.5])
    K = dipolarkernel(t,r)
    L = regoperator(r,2)
    V = K@P

    Pfit,_ = fitregmodel(V,K,r,'tikhonov',alpha=0.05,nonnegativity=False,renormalize=False)
    Pfit_ref = np.linalg.solve(K.T@K + 0.05**2*L.T@L, K.T@V)

    assert np.max(abs(Pfit - Pfit_ref)) < 1e-10
#============================================================

def assert_confidence_intervals(Puq,Pfit):
    
    Pci50 = Puq.ci(50) 
    Pci95 = Puq.ci(95) 
    P95lb = Pci95[:,0]
    P95ub = Pci95[:,1]
    P50lb = Pci50[:,0]
    P50ub = Pci50[:,1]
    errors = []
    if not np.all(P95lb <= Pfit) and not np.all(P50lb <= Pfit):
        errors.append("Some fitted values are below the lower bound of the confidence intervals.")
    if not np.all(P95ub >= Pfit) and not np.all(P50lb >= Pfit):
        errors.append("Some fitted values are over the upper bound of the confidence intervals.")
    if not np.all(P95lb <= P50lb):
        errors.append("The 50%-CI has lower values than the 95%-CI")
    if not np.all(P95ub >= P50ub):
        errors.append("The 50%-CI has larger values than the 95%-CI")

    assert not errors, "Errors occured:\n{}".format("\n".join(errors))

def test_confidence_intervals():
#============================================================
    "Check that the confidence intervals are correctly computed"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.5])
    K = dipolarkernel(t,r)
    V = K@P

    Pfit,Puq = fitregmodel(V,K,r,'tikhonov','aic')

    assert_confidence_intervals(Puq,Pfit)
#============================================================

def test_renormalize():
#============================================================
    "Check that renormalization of the distribution can be correctly enabled/disabled"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.5])
    K = dipolarkernel(t,r)
    V = K@P

    Praw,_ = fitregmodel(V,K,r,'tikhonov','aic',renormalize = False)
    Prenorm,_ = fitregmodel(V,K,r,'tikhonov','aic',renormalize = True)

    assert max(abs(Praw - Prenorm*np.trapz(Praw,r))) < 1e-10
#============================================================

def test_scale_agnostic():
#============================================================
    "Check agnosticity of the w.r.t. an arbitrary scale of dipolar signal"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.5])
    K = dipolarkernel(t,r)
    scale = 1e8
    V = scale*K@P

    Pfit,_ = fitregmodel(V,K,r,'tikhonov','aic',renormalize = False)

    assert max(abs(P - Pfit/scale)) < 2e-6
#============================================================

def test_nonuniform_r():
#============================================================
    "Check that fit works correctly with non-uniform distance vectors"

    t = np.linspace(0,4,300)
    r = np.sqrt(np.linspace(1,7**2,200))
    P = dd_gauss(r,[3,0.5])
    K = dipolarkernel(t,r)
    V = K@P

    Pfit,_ = fitregmodel(V,K,r,'tikhonov','aic')
    Vfit = K@Pfit
    assert Vfit[0]==1
#============================================================
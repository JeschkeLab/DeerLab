
import numpy as np
from deerlab import fitregmodel,dipolarkernel, regoperator, whitegaussnoise
from deerlab.dd_models import dd_gauss,dd_gauss2
from deerlab.utils import ovl

def assert_solver(regtype,solver):
#============================================================

    np.random.seed(1)
    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)
    fit = fitregmodel(V,K,r,regtype,'aic',solver=solver)

    assert ovl(P,fit.P) > 0.95 # more than 95% overlap
#============================================================

def test_tikh_solver_fnnls():
#============================================================
    "Check that 'fnnls' solves a Tikhonov regularization problem"

    assert_solver('tikhonov','fnnls')
#============================================================

def test_tikh_solver_nnlsbpp():
#============================================================
    "Check that 'nnlsbpp' solves a Tikhonov regularization problem"

    assert_solver('tikhonov','nnlsbpp')
#============================================================

def test_tikh_solver_cvx():
#============================================================
    "Check that 'cvx' solves a Tikhonov regularization problem"

    assert_solver('tikhonov','cvx')
#============================================================

def test_tv_solver_fnnls():
#============================================================
    "Check that 'fnnls' solves a TV regularization problem"

    assert_solver('tv','fnnls')
#============================================================

def test_tv_solver_nnlsbpp():
#============================================================
    "Check that 'nnlsbpp' solves a TV regularization problem"

    assert_solver('tv','nnlsbpp')
#============================================================

def test_tv_solver_cvx():
#============================================================
    "Check that 'cvx' solves a TV regularization problem"

    assert_solver('tv','cvx')
#============================================================


def test_huber_solver_fnnls():
#============================================================
    "Check that 'fnnls' solves a Huber regularization problem"

    assert_solver('huber','fnnls')
#============================================================

def test_huber_solver_nnlsbpp():
#============================================================
    "Check that 'nnlsbpp' solves a Huber regularization problem"

    assert_solver('huber','nnlsbpp')
#============================================================

def test_huber_solver_cvx():
#============================================================
    "Check that 'cvx' solves a Huber regularization problem"

    assert_solver('huber','cvx')
#============================================================

def test_tikh_with_noise():
#============================================================
    "Check the TV regularization method"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(1,5,100)
    P = dd_gauss(r,[3,0.08])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)

    fit = fitregmodel(V,K,r,'tikhonov','aic')
    
    assert ovl(P,fit.P) > 0.95 # more than 95% overlap
#============================================================

def test_tv_with_noise():
#============================================================
    "Check the TV regularization method"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(1,5,100)
    P = dd_gauss(r,[3,0.08])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)

    fit = fitregmodel(V,K,r,'tv','aic')
    
    assert ovl(P,fit.P) > 0.95 # more than 95% overlap
#============================================================

def test_huber_with_noise():
#============================================================
    "Check the Huber regularization method"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(1,5,100)
    P = dd_gauss(r,[3,0.08])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)

    fit = fitregmodel(V,K,r,'huber','aic')
    
    assert ovl(P,fit.P) > 0.95 # more than 95% overlap
#============================================================


def test_unconstrained():
#============================================================
    "Check that unconstrained distributions are correctly fitted"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    L = regoperator(r,2)
    V = K@P

    fit = fitregmodel(V,K,r,'tikhonov',regparam=0.05,nonnegativity=False,renormalize=False)
    Pfit_ref = np.linalg.solve(K.T@K + 0.05**2*L.T@L, K.T@V)

    assert ovl(Pfit_ref,fit.P) > 0.99 # more than 99% overlap
#============================================================

def generate_global_dataset():
    
    t1 = np.linspace(0,3,100)
    t2 = np.linspace(0,4,200)
    r = np.linspace(2.5,5,80)
    P = dd_gauss2(r,[3.7,0.2,0.5,4.3,0.1,0.5])
    K1 = dipolarkernel(t1,r)
    K2 = dipolarkernel(t2,r)
    np.random.seed(1)
    V1 = K1@P + whitegaussnoise(t1,0.05)
    np.random.seed(2)
    V2 = K2@P + whitegaussnoise(t2,0.05)

    return r,P,V1,V2,K1,K2

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
    if not np.all(np.minimum(0,P95lb)==0):
        errors.append("The non-negativity constraint is not working.")
    assert not errors, "Errors occured:\n{}".format("\n".join(errors))

def test_confinter_tikh():
#============================================================
    "Check that the confidence intervals are correctly computed using Tikhonov regularization"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P

    fit = fitregmodel(V,K,r,'tikhonov','aic')

    assert_confidence_intervals(fit.uncertainty,fit.P)
#============================================================

def test_confinter_tv():
#============================================================
    "Check that the confidence intervals are correctly computed using TV regularization"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P

    fit = fitregmodel(V,K,r,'tv','aic')

    assert_confidence_intervals(fit.uncertainty,fit.P)
#============================================================


def test_confinter_huber():
#============================================================
    "Check that the confidence intervals are correctly computed using Huber regularization"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P

    fit = fitregmodel(V,K,r,'huber','aic')

    assert_confidence_intervals(fit.uncertainty,fit.P)
#============================================================

def test_confinter_global():
#============================================================
    "Check that the confidence intervals are correctly computed even with multiple datasets"

    r,_,V1,V2,K1,K2 = generate_global_dataset()

    fit = fitregmodel([V1,V2],[K1,K2],r,'tikhonov','aic')

    assert_confidence_intervals(fit.uncertainty,fit.P)
#============================================================

def test_renormalize():
#============================================================
    "Check that renormalization of the distribution can be correctly enabled/disabled"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P

    fitraw = fitregmodel(V,K,r,'tikhonov','aic',renormalize = False)
    fitrenorm = fitregmodel(V,K,r,'tikhonov','aic',renormalize = True)

    assert max(abs(fitraw.P - fitrenorm.P*np.trapz(fitraw.P,r))) < 1e-8
#============================================================

def test_scale_agnostic():
#============================================================
    "Check agnosticity of the w.r.t. an arbitrary scale of dipolar signal"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    scale = 1e8
    V = scale*K@P

    fit = fitregmodel(V,K,r,'tikhonov','aic',renormalize = False)

    assert max(abs(P - fit.P/scale)) < 1e-3
#============================================================

def test_scale_fit():
#============================================================
    "Check agnosticity of the w.r.t. an arbitrary scale of dipolar signal"

    t = np.linspace(-2,4,300)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    scale = 1e8
    V = scale*K@P

    fit = fitregmodel(V,K,r,'tikhonov','aic',renormalize = False)

    assert abs(1 - scale/fit.scale) < 1e-4
#============================================================

def test_nonuniform_r():
#============================================================
    "Check that fit works correctly with non-uniform distance vectors"

    t = np.linspace(0,4,300)
    r = np.sqrt(np.linspace(1,7**2,200))
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P

    fit = fitregmodel(V,K,r,'tikhonov','aic')
    Vfit = K@fit.P
    assert abs(Vfit[0] - 1) < 1e-6
#============================================================


def test_obir():
#============================================================
    "Check that the OBIR algorithm runs and returns a correct result"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(2,7,100)
    P = dd_gauss2(r,[3.7,0.2,0.5,4.3,0.1,0.5])
    K = dipolarkernel(t,r)

    V = K@P + whitegaussnoise(t,0.05)

    fit = fitregmodel(V,K,r,'tikhonov','aic',obir=True)
    
    assert ovl(P,fit.P) > 0.75 # more than 80% overlap
#============================================================

def test_obir_global():
#============================================================
    "Check that the OBIR algorithm runs and returns a correct result with multiple datasets"

    r,P,V1,V2,K1,K2 = generate_global_dataset()

    fit = fitregmodel([V1,V2],[K1,K2],r,'tikhonov','aic',obir=True)
    
    assert ovl(P,fit.P) > 0.75 # more than 80% overlap
#============================================================


def test_goodnes_of_fit():
#============================================================
    "Check the goodness-of-fit statistics are correct"
    
    np.random.seed(1)
    t = np.linspace(0,8,500)
    r = np.linspace(2,5,300)
    P = dd_gauss(r,[4.5,0.2])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)

    fit = fitregmodel(V,K,r,'tikhonov','aic')
    stats = fit.stats

    assert abs(stats['chi2red'] - 1) < 5e-2 and abs(stats['R2'] - 1) < 5e-2
#============================================================


def test_goodnes_of_fit_global():
#============================================================
    "Check the goodness-of-fit statistics are correct with multiple signals"

    t1 = np.linspace(0,3,100)
    t2 = np.linspace(0,4,200)
    r = np.linspace(2.5,5,80)
    P = dd_gauss2(r,[3.7,0.2,0.5,4.3,0.1,0.5])
    K1 = dipolarkernel(t1,r)
    K2 = dipolarkernel(t2,r)
    np.random.seed(1)
    V1 = K1@P + whitegaussnoise(t1,0.02)
    np.random.seed(2)
    V2 = K2@P + whitegaussnoise(t2,0.01)

    fit = fitregmodel([V1,V2],[K1,K2],r,'tikhonov','aic')
    stats = fit.stats

    assert (abs(stats[0]['R2'] - 1) < 5e-2) and (abs(stats[1]['R2'] - 1) < 5e-2)
#============================================================

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

    fit = fitregmodel([V1,V2],[K1,K2],r,'tikhonov','aic',weights=[1,1],renormalize=False)

    assert max(abs(np.asarray(scales)/np.asarray(fit.scale) - 1)) < 1e-2 
#============================================================

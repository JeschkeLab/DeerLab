import numpy as np
from deerlab import dipolarkernel,dd_gauss,dd_gauss2,snlls
from deerlab.bg_models import bg_exp
from deerlab.utils import ovl

def assert_multigauss_SNLLS_problem(nonlinearconstr=True, linearconstr=True):
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    K = dipolarkernel(t,r)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    def Kmodel(p,t,r):
        # Unpack parameters
        r1, w1, r2, w2  = p
        # Generate basic kernel
        K0 = dipolarkernel(t,r)
        # Get Gauss basis functions
        P1 = dd_gauss(r,[r1,w1])
        P2 = dd_gauss(r,[r2,w2])
        # Combine all non-linear functions into one
        K = np.zeros((len(t),2))
        K[:,0] = K0@P1
        K[:,1] = K0@P2
        return K

    # Non-linear parameters
    # nlpar = [r1 w1 r2 w2]
    nlpar0 = [3.2, 0.2, 4.2, 0.3]
    if nonlinearconstr:
        lb = [1, 0.1, 1, 0.1]
        ub = [20, 5, 20, 5]
    else:
        lb = []
        ub = []
    # Linear parameters
    if linearconstr:
        lbl = [0, 0]
        ubl = [1, 1]
    else:
        lbl = []
        ubl = []    

    # Separable LSQ fit
    fit = snlls(V,lambda p: Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl, reg=False, uqanalysis=False)
    nonlinfit = fit.nonlin
    linfit = fit.lin
    parout = [nonlinfit[0], nonlinfit[1], linfit[0], nonlinfit[2], nonlinfit[3], linfit[1]] 
    parout = np.asarray(parout)
    parin = np.asarray(parin)
    assert np.all(abs(parout - parin) < 1e-1)


def test_const_nonlin_const_lin():
#=======================================================================
    "Check SNLLS of a nonlinear-constrained + linear-constrained problem"

    assert_multigauss_SNLLS_problem(nonlinearconstr=True, linearconstr=True)
#=======================================================================

def test_const_nonlin_unconst_lin():
#=======================================================================
    "Check SNLLS of a nonlinear-constrained + linear-unconstrained problem"

    assert_multigauss_SNLLS_problem(nonlinearconstr=True, linearconstr=False)
#=======================================================================

def test_unconst_nonlin_const_lin():
#=======================================================================
    "Check SNLLS of a nonlinear-unconstrained + linear-constrained problem"

    assert_multigauss_SNLLS_problem(nonlinearconstr=False, linearconstr=True)
#=======================================================================

def test_unconst_nonlin_unconst_lin():
#=======================================================================
    "Check SNLLS of a nonlinear-unconstrained + linear-unconstrained problem"

    assert_multigauss_SNLLS_problem(nonlinearconstr=False, linearconstr=False)
#=======================================================================

def test_regularized():
#=======================================================================
    "Check SNLLS of a nonlinear-constrained + linear-regularized problem"
    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,lam),nlpar0,lb,ub,lbl,ubl, uqanalysis=False)
    Pfit = fit.lin

    assert  np.max(abs(P - Pfit)) < 1e-2
#=======================================================================


def assert_confidence_intervals(pci50,pci95,pfit,lb,ub):
    
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
    assert not errors, "Errors occured:\n{}".format("\n".join(errors))

def test_confinter_linear():
#=======================================================================
    "Check that the confidence intervals of the linear parameters are correct"
    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = np.full(len(r), np.inf)
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,lam),nlpar0,lb,ub,lbl,ubl)
    Pfit = fit.lin
    uq = fit.uncertainty
    Pci50 = uq.ci(50,'lin')
    Pci95 = uq.ci(95,'lin')

    assert_confidence_intervals(Pci50,Pci95,Pfit,lbl,ubl)
#=======================================================================

def test_confinter_nonlinear():
#=======================================================================
    "Check that the confidence intervals of the non-linear parameters are correct"
    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = np.full(len(r), np.inf)
    # Separable LSQ fit

    fit = snlls(V,lambda lam: dipolarkernel(t,r,lam),nlpar0,lb,ub,lbl,ubl)
    parfit = fit.nonlin
    uq = fit.uncertainty
    parci50 = uq.ci(50,'nonlin')
    parci95 = uq.ci(95,'nonlin')

    assert_confidence_intervals(parci50,parci95,parfit,lb,ub)
#=======================================================================

def test_regularized_global():
#=======================================================================
    "Check global SNLLS of a nonlinear-constrained + linear-regularized problem"
        
    t1 = np.linspace(0,3,150)
    t2 = np.linspace(0,4,200)
    r = np.linspace(2.5,5,80)
    P = dd_gauss2(r,[3.7,0.5,0.5,4.3,0.3,0.5])
    kappa = 0.50
    lam1 = 0.25
    lam2 = 0.35
    K1 = dipolarkernel(t1,r,lam1,bg_exp(t1,kappa))
    K2 = dipolarkernel(t2,r,lam2,bg_exp(t2,kappa))
    V1 = K1@P
    V2 = K2@P

    # Global non-linear model
    def globalKmodel(par):
        # Unpack parameters
        kappa,lam1,lam2 = par
        K1 = dipolarkernel(t1,r,lam1,bg_exp(t1,kappa))
        K2 = dipolarkernel(t2,r,lam2,bg_exp(t2,kappa))
        return K1, K2

    # Non-linear parameters
    # [kappa lambda1 lambda2]
    par0 = [0.5, 0.5, 0.5]
    lb = [0, 0, 0]
    ub = [1, 1, 1]
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []

    # Separable LSQ fit
    fit = snlls([V1,V2],globalKmodel,par0,lb,ub,lbl,ubl, uqanalysis=False)
    Pfit = fit.lin

    assert  ovl(P,Pfit) > 0.9
#=======================================================================


    
def assert_solver(solver):    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,lam),nlpar0,lb,ub,lbl,ubl,nnlsSolver=solver, uqanalysis=False)
    Pfit = fit.lin

    assert  np.max(abs(P - Pfit)) < 1e-2

def test_nnls_cvx():
#=======================================================================
    "Check that the 'cvx' algorithm used for the NNLS subproblem work"

    assert_solver(solver='cvx')
#=======================================================================

def test_nnls_fnnls():
#=======================================================================
    "Check that the 'cvx' algorithm used for the NNLS subproblem work"

    assert_solver(solver='fnnls')
#=======================================================================

def test_nnls_nnlsbpp():
#=======================================================================
    "Check that the 'nnlsbpp' algorithm used for the NNLS subproblem work"

    assert_solver(solver='nnlsbpp')
#=======================================================================

def test_goodnes_of_fit():
#============================================================
    "Check the goodness-of-fit statistics are correct"
        
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,lam)
    parin = [3.5, 0.15, 0.6, 4.5, 0.2, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,lam),nlpar0,lb,ub,lbl,ubl, uqanalysis=False)
    stats = fit.stats

    assert abs(stats['chi2red'] - 1) < 5e-2 and abs(stats['R2'] - 1) < 5e-2
#============================================================

def assert_reg_type(regtype):    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,100)
    lam = 0.25
    K = dipolarkernel(t,r,lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,lam),nlpar0,lb,ub,lbl,ubl,regtype = regtype, uqanalysis=False)
    Pfit = fit.lin
    
    assert  np.max(abs(P - Pfit)) < 4e-2

def test_reg_tikh():
#=======================================================================
    "Check SNNLS using Tikhonov regularization works"

    assert_reg_type(regtype='tikhonov')
#=======================================================================

def test_reg_tv():
#=======================================================================
    "Check SNNLS using Tikhonov regularization works"

    assert_reg_type(regtype='tv')
#=======================================================================

def test_reg_huber():
#=======================================================================
    "Check SNNLS using Tikhonov regularization works"

    assert_reg_type(regtype='huber')
#=======================================================================
import numpy as np
from deerlab import dipolarkernel,dd_gauss,dd_gauss2,snlls,whitegaussnoise
from deerlab.bg_models import bg_exp
from deerlab.utils import ovl, skip_on, assert_docstring
import pytest

def assert_multigauss_SNLLS_problem(nonlinearconstr=True, linearconstr=True):
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    K = dipolarkernel(t,r)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P

    def Kmodel(p,t,r):
        # Unpack parameters
        r1, w1, r2, w2  = p
        # Generate basic kernel
        K0 = dipolarkernel(t,r)
        # Get Gauss basis functions
        P1 = dd_gauss(r,r1,w1)
        P2 = dd_gauss(r,r2,w2)
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
    fit = snlls(V,lambda p: Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl, reg=False, uq=False)
    nonlinfit = fit.nonlin
    linfit = fit.lin
    parout = [nonlinfit[0], nonlinfit[1], linfit[0], nonlinfit[2], nonlinfit[3], linfit[1]] 
    parout = np.asarray(parout)
    parin = np.asarray(parin)
    assert np.max(abs(parout - parin) < 1e-1)

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
    t = np.linspace(0,4,100)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl, uq=False)
    Pfit = fit.lin

    assert  ovl(P,Pfit) > 0.95
#=======================================================================


def assert_confidence_intervals(pci50,pci95,pfit,lb,ub):
    
    p95lb = pci95[:,0]
    p95ub = pci95[:,1]
    p50lb = pci50[:,0]
    p50ub = pci50[:,1]
    errors = []
    if not np.all(p95lb <= pfit) and not np.all(p50lb <= pfit):
        errors.append("Some fitted values are below the lower bound of the confidence intervals.")
    if not np.all(p95ub >= pfit) and not np.all(p50ub >= pfit):
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

def test_confinter_linear():
#=======================================================================
    "Check that the confidence intervals of the linear parameters are correct"
    
    # Prepare test data
    r = np.linspace(1,8,150)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P + whitegaussnoise(t,0.05,seed=1)

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = np.full(len(r), np.inf)
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl)
    Pfit =  np.round(fit.lin,6)
    uq = fit.linUncert
    Pci50 = np.round(uq.ci(50),6)
    Pci95 = np.round(uq.ci(95),6)

    assert_confidence_intervals(Pci50,Pci95,Pfit,lbl,ubl)
#=======================================================================

def test_confinter_nonlinear():
#=======================================================================
    "Check that the confidence intervals of the non-linear parameters are correct"
    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
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
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ubl)
    parfit = fit.nonlin
    uq = fit.nonlinUncert
    parci50 = np.atleast_2d(uq.ci(50))
    parci95 = np.atleast_2d(uq.ci(95))

    assert_confidence_intervals(parci50,parci95,parfit,lb,ub)
#=======================================================================

def test_confinter_model():
#=======================================================================
    "Check that the confidence intervals of the fitted model are correct"
    
    # Prepare test data
    r = np.linspace(1,8,150)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P + whitegaussnoise(t,0.05,seed=1)

    nlpar0 = 0.2
    lb = 0
    ub = 1
    lbl = np.full(len(r), 0)
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl)
    Vfit =  fit.model
    Vuq = fit.modelUncert
    Vci50 = Vuq.ci(50)
    Vci95 = Vuq.ci(95)

    Vlb = np.full(len(t), -np.inf)
    Vub = np.full(len(t), np.inf)

    assert_confidence_intervals(Vci50,Vci95,Vfit,Vlb,Vub)
#=======================================================================

def test_regularized_global():
#=======================================================================
    "Check global SNLLS of a nonlinear-constrained + linear-regularized problem"
        
    t1 = np.linspace(0,3,150)
    t2 = np.linspace(0,4,200)
    r = np.linspace(2.5,5,80)
    P = dd_gauss2(r,3.7,0.5,0.5,4.3,0.3,0.5)
    kappa = 0.50
    lam1 = 0.25
    lam2 = 0.35
    K1 = dipolarkernel(t1,r,mod=lam1,bg=bg_exp(t1,kappa))
    K2 = dipolarkernel(t2,r,mod=lam2,bg=bg_exp(t2,kappa))
    V1 = K1@P
    V2 = K2@P

    # Global non-linear model
    def globalKmodel(par):
        # Unpack parameters
        kappa,lam1,lam2 = par
        K1 = dipolarkernel(t1,r,mod=lam1,bg=bg_exp(t1,kappa))
        K2 = dipolarkernel(t2,r,mod=lam2,bg=bg_exp(t2,kappa))
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
    fit = snlls([V1,V2],globalKmodel,par0,lb,ub,lbl,ubl, uq=False, weights=[1,1])
    Pfit = fit.lin

    assert  ovl(P,Pfit) > 0.9
#=======================================================================

def assert_solver(solver):    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
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
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ubl,nnlsSolver=solver, uq=False)
    Pfit = fit.lin

    assert  ovl(P,Pfit) > 0.95

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

def test_goodness_of_fit():
#============================================================
    "Check the goodness-of-fit statistics are correct"
        
    # Prepare test data
    r = np.linspace(2,5,150)
    t = np.linspace(-0.2,4,100)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.15, 0.6, 4.5, 0.2, 0.4]
    P = dd_gauss2(r,*parin)
    sigma = 0.03
    V = K@P + whitegaussnoise(t,sigma,seed=2,rescale=True)

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ubl, noiselvl=sigma, uq=False)
    stats = fit.stats

    assert abs(stats['chi2red'] - 1) < 0.05
#============================================================


def test_goodness_of_fit_scaled():
#============================================================
    "Check the goodness-of-fit statistics are correct even with arbitrary scaling"
        
    # Prepare test data
    r = np.linspace(2,5,150)
    t = np.linspace(-0.2,4,300)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.15, 0.6, 4.5, 0.2, 0.4]
    P = dd_gauss2(r,*parin)
    V0 = 1e3
    V = V0*K@P 
    sigma = V0*0.03
    V += whitegaussnoise(t,sigma,seed=1,rescale=True)

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = []

    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ubl, noiselvl=sigma, uq=False)
    stats = fit.stats

    assert abs(stats['chi2red'] - 1) < 0.05
#============================================================

def test_size_control():
#============================================================
    "Check that the function checks the size of the input and output"

    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,100)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P
    nlpar0 = 0.2
    lb = 0
    ub = 1
    lbl = np.zeros(len(r))

    twrong = np.linspace(0,4,200)
    with pytest.raises(RuntimeError):
        fit = snlls(V,lambda lam: dipolarkernel(twrong,r,mod=lam),nlpar0,lb,ub,lbl,uq=False)
#============================================================

def test_reg_tikhonov():
#============================================================
    "Check that Tikhonov regularization of linear problem works"

    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,100)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P

    # Non-linear parameters
    # nlpar = [lam]
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,uq=False)
    Pfit = fit.lin
    
    assert  ovl(P,Pfit) > 0.95
#============================================================


@skip_on('_tkinter.TclError', reason="A problem with the Tk backend occured")
def test_plot():
# ======================================================================
    "Check that the plot method works"

    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),par0=0.2,lb=0,ub=1,lbl=lbl, uq=False)
    fig = fit.plot()
    
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ======================================================================

def test_cost_value():
#============================================================
    "Check that the cost value is properly returned"

    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P
    # Non-linear parameters
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    ubl = np.full(len(r), np.inf)
    # Separable LSQ fit
    fit = snlls(V,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ubl)

    assert isinstance(fit.cost,float) and np.round(fit.cost/np.sum(fit.residuals**2),5)==1
#============================================================

def test_confinter_values():
# ======================================================================
    "Check that the values of the confidence intervals are correct"

    np.random.seed(0)
    A0  = np.random.rand(300,2)
    A = lambda p: 1 - p[0] + A0**p[1]
    pnonlin = np.array([0.5,2])
    plin = np.array([0.75,5.5])
    y = A(pnonlin)@plin + 0.3*np.random.randn(300)

    # Reference confidence intervals (calculated using lmfit package, see #116)
    cov = [99.73,95.45,68.27]
    a_ci_ref = [[0.4636991758611843, 0.5411943256617782], 
                [0.47730341575508245, 0.5286896152296477], 
                [0.4905189733270308, 0.5161245529171482]]
    b_ci_ref = [[1.7876921926935445, 2.1089296764536907], 
                [1.8384482433779732, 2.0515346381926847], 
                [1.8899306688797555, 1.9961521871803736]]

    fit = snlls(y,A,[0.5,0.2],reg=False)
    a_ci = [fit.nonlinUncert.ci(cov[i])[0,:] for i in range(3)]
    b_ci = [fit.nonlinUncert.ci(cov[i])[1,:] for i in range(3)]

    ci_match = lambda ci,ci_ref,truth:np.max(abs(np.array(ci) - np.array(ci_ref)))/truth < 0.01
    assert ci_match(a_ci,a_ci_ref,pnonlin[0]) & ci_match(b_ci,b_ci_ref,pnonlin[1])
# ======================================================================

def test_confinter_scaling():
#============================================================
    "Check that the confidence intervals are agnostic w.r.t. scaling"

    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,50)
    lam = 0.25
    K = dipolarkernel(t,r,mod=lam)
    parin = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,*parin)
    V = K@P + whitegaussnoise(t,0.01,seed=1)
    # Non-linear parameters
    nlpar0 = 0.2
    lb = 0
    ub = 1
    # Linear parameters: non-negativity
    lbl = np.zeros(len(r))
    V0_1 = 1
    V0_2 = 1e8

    # Separable LSQ fit
    fit1 = snlls(V*V0_1,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ftol=1e-3)
    fit2 = snlls(V*V0_2,lambda lam: dipolarkernel(t,r,mod=lam),nlpar0,lb,ub,lbl,ftol=1e-3)

    # Assess linear parameter uncertainties
    ci1 = fit1.linUncert.ci(95)
    ci2 = fit2.linUncert.ci(95)
    ci1[ci1==0] = 1e-16
    ci2[ci2==0] = 1e-16

    assert np.max(abs(ci2/V0_2 - ci1)) < 1e-6

    # Assess nonlinear parameter uncertainties
    ci1 = fit1.nonlinUncert.ci(95)
    ci2 = fit2.nonlinUncert.ci(95)
    ci1[ci1==0] = 1e-16
    ci2[ci2==0] = 1e-16

    assert np.max(abs(ci2 - ci1)) < 1e-6    
#============================================================

def test_global_weights():
# ======================================================================
    "Check that the global weights properly work when specified"

    t = np.linspace(-0.3,5,300)
    r = np.linspace(2,6,150)

    P1 = dd_gauss(r,3,0.2)
    P2 = dd_gauss(r,5,0.2)

    K = dipolarkernel(t,r,mod=0.2)

    scales = [1e3, 1e9]
    sigma1 = 0.001
    V1 = K@P1 + whitegaussnoise(t,sigma1,seed=1)
    sigma2 = 0.001
    V2 = K@P2 + whitegaussnoise(t,sigma2,seed=1)

    V1 = scales[0]*V1
    V2 = scales[1]*V2

    Kmodel= lambda lam: [dipolarkernel(t,r,mod=lam)]*2
    fit1 = snlls([V1,V2],Kmodel,par0=[0.2],lb=0,ub=1,lbl=np.zeros_like(r),weights=[1,1e-10])
    fit2 = snlls([V1,V2],Kmodel,par0=[0.2],lb=0,ub=1,lbl=np.zeros_like(r),weights=[1e-10,1])

    assert ovl(P1,fit1.lin) > 0.93 and ovl(P2,fit2.lin) > 0.93
# ======================================================================

def test_global_weights_default():
# ======================================================================
    "Check the correct fit of two signals when one is of very low quality"

    t = np.linspace(0,5,300)
    r = np.linspace(2,6,90)
    param = [4.5, 0.25]
    P = dd_gauss(r,*param)

    K = dipolarkernel(t,r,mod=0.2)
    scales = [1e3,1e9]
    V1 = scales[0]*(K@P + whitegaussnoise(t,0.001,seed=1))
    V2 = scales[1]*(K@P + whitegaussnoise(t,0.1,seed=1))

    Kmodel= lambda lam: [dipolarkernel(t,r,mod=lam)]*2
    fit = snlls([V1,V2],Kmodel,par0=[0.2],lb=0,ub=1,lbl=np.zeros_like(r))

    assert ovl(P,fit.lin) > 0.93
# ======================================================================

def test_extrapenalty():
# ======================================================================
    "Check that an additional penalty can be passed correctly"

    t = np.linspace(0,5,300)
    r = np.linspace(2,6,90)
    P = dd_gauss(r,4.5, 0.25)
    K = dipolarkernel(t,r,mod=0.2)
    V = K@P + whitegaussnoise(t,0.001,seed=1)
    dr = np.mean(np.diff(r))
    beta = 0.05
    compactness_penalty = lambda _,plin: beta*np.sqrt(plin*(r - np.trapz(plin*r,r))**2*dr)
    Kmodel= lambda lam: dipolarkernel(t,r,mod=lam)
    fit = snlls(V,Kmodel,par0=0.2,lb=0,ub=1,lbl=np.zeros_like(r),extrapenalty=compactness_penalty)

    assert ovl(P,fit.lin) > 0.95
# ======================================================================

def test_multiple_penalties():
# ======================================================================
    "Check that multiple additional penaltyies can be passed correctly"

    t = np.linspace(0,5,300)
    r = np.linspace(2,6,90)
    P = dd_gauss(r,4.5, 0.25)
    param = 0.2
    K = dipolarkernel(t,r,mod=param)
    V = K@P + whitegaussnoise(t,0.001,seed=1)
    dr = np.mean(np.diff(r))
    beta = 0.05
    R = 0.5
    compactness_penalty = lambda pnonlin,plin: beta*np.sqrt(plin*(r - np.trapz(plin*r,r))**2*dr)
    radial_penalty = lambda pnonlin,plin: 1/R**2*(np.linalg.norm((pnonlin-param)/param-R))**2

    Kmodel= lambda lam: dipolarkernel(t,r,mod=lam)
    fit0 = snlls(V,Kmodel,par0=0.2,lb=0,ub=1,lbl=np.zeros_like(r),extrapenalty=[compactness_penalty])
    fitmoved = snlls(V,Kmodel,par0=0.2,lb=0,ub=1,lbl=np.zeros_like(r),extrapenalty=[compactness_penalty,radial_penalty])
    
    assert ovl(P,fit0.lin) > ovl(P,fitmoved.lin)
# ======================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(snlls)
# ======================================================================

def test_frozen_param():
# ======================================================================
    "Check that linear and nonlinear parameters can be frozen during the optimization"
    r = np.linspace(0,6,300)
    def Amodel(p):
        mean1,mean2,std1,std2 = p
        return np.atleast_2d([dd_gauss.nonlinmodel(r,mean1,std1), dd_gauss.nonlinmodel(r,mean2,std2)]).T

    x = np.array([0.5,0.6])
    y = Amodel([3,5,0.2,0.3])@x

    nonlin_frozen = [None,5,None,None]
    lin_frozen = [0.5,None]
    fit = snlls(y,Amodel,par0=[3.2,5.2,0.2,0.3],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0])
    fit_frozen = snlls(y,Amodel,par0=[3.2,5.2,0.2,0.3],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0],
            nonlin_frozen=nonlin_frozen, lin_frozen=lin_frozen)
    
    assert np.allclose(fit_frozen.model,y,atol=1e-2) and np.allclose(fit_frozen.model,fit.model,atol=1e-2)
# ======================================================================

def test_frozen_Nparam():
# ======================================================================
    "Check that the correct number of linear and nonlinear parameters are return even when freezing"
    r = np.linspace(0,6,90)
    def Amodel(p):
        mean1,mean2,std1,std2 = p
        return np.atleast_2d([dd_gauss.nonlinmodel(r,mean1,std1), dd_gauss.nonlinmodel(r,mean2,std2)]).T
    x = np.array([0.5,0.6])
    y = Amodel([3,5,0.2,0.3])@x
    nonlin_frozen = [None,5,None,None]
    lin_frozen = [0.5,None]

    fit = snlls(y,Amodel,par0=[2,4,0.2,0.2],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0],
            nonlin_frozen=nonlin_frozen, lin_frozen=lin_frozen)
    
    fit_frozen = snlls(y,Amodel,par0=[2,4,0.2,0.2],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0])
    
    assert len(fit.nonlin)==4 and len(fit.lin)==2 and len(fit_frozen.nonlin)==4 and len(fit_frozen.lin)==2
# ======================================================================

def test_complex_model_complex_data():
# ======================================================================
    "Check the fit of a real-valued model to complex-valued data"

    x = np.linspace(0,7,100)
    def model(p):
        phase, center, std = p
        y = dd_gauss(x,center, std)
        y = y*np.exp(-1j*phase)
        return y

    y = model([np.pi/5, 3, 0.5])

    fitResult = snlls(y,model,par0=[2*np.pi/5,4,0.2],lb=[-np.pi,1,0.05],ub=[np.pi,6,5])

    assert np.allclose(fitResult.model.real,y.real) and np.allclose(fitResult.model.imag, y.imag)
# ======================================================================

def test_complex_model_uncertainty():
# ======================================================================
    "Check the fit of a real-valued model to complex-valued data"

    x = np.linspace(0,7,100)
    def model(p):
        phase, center, std = p
        y = dd_gauss(x,center, std)
        y = y*np.exp(-1j*phase)
        return y

    y = model([np.pi/5, 3, 0.5])     

    y = y + whitegaussnoise(x,0.01)
    y = y + 1j*whitegaussnoise(x,0.05)

    fitResult = snlls(y,model,par0=[2*np.pi/5,4,0.2],lb=[-np.pi,1,0.05],ub=[np.pi,6,5])
    ciwidth = np.sum(fitResult.modelUncert.ci(95)[:,1] - fitResult.modelUncert.ci(95)[:,0])

    assert (ciwidth.real < ciwidth.imag).all()
# ======================================================================
 
def test_masking():
# ======================================================================
    "Check that datapoints can be masked out"

    x = np.linspace(0,7,100)
    def model(p):
        center, std = p
        y = dd_gauss(x,center, std)
        return y

    mask = np.ones_like(x).astype(bool)   
    mask[(x>2.5) & (x<3.5)] = False 

    y = model([3, 0.5])  
    yref = y.copy()
    y[~mask] = 0 

    fitmasked = snlls(y,model,par0=[4,0.2],lb=[1,0.05],ub=[6,5], mask=mask)
    fit = snlls(y,model,par0=[4,0.2],lb=[1,0.05],ub=[6,5])

    assert np.allclose(fitmasked.model,yref) and not np.allclose(fit.model,yref) 
# ======================================================================

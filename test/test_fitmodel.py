
import numpy as np
from numpy import inf
import deerlab as dl
from deerlab import dipolarkernel, whitegaussnoise, fitmodel
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp, bg_hom3d
from deerlab.ex_models import ex_4pdeer, ex_5pdeer, ex_7pdeer, ex_ovl4pdeer
import deerlab as dl
from deerlab.utils import ovl


def assert_experiment_model(model):
# --------------------------------------------------------------------
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,50)
    P = dd_gauss(r,[4.5, 0.25])

    parIn = model.start
    pathways = model(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)
    K = dipolarkernel(t,r,pathways=pathways,bg=Bmodel)
    V = K@P

    fit = fitmodel(V,t,r,'P',bg_exp,model,uq=None)

    assert ovl(P,fit.P) > 0.90
# --------------------------------------------------------------------

def test_4pdeer():
# ======================================================================
    "Check that the fit of a 4-pulse DEER signal is correct"
        
    assert_experiment_model(ex_4pdeer)
# ======================================================================

def test_5pdeer():
# ======================================================================
    "Check that the fit of a 4-pulse DEER signal is correct"
        
    assert_experiment_model(ex_5pdeer)
# ======================================================================

def test_7pdeer():
# ======================================================================
    "Check that the fit of a 4-pulse DEER signal is correct"
        
    assert_experiment_model(ex_7pdeer)
# ======================================================================

def test_ovl4pdeer():
# ======================================================================
    "Check that the fit of a 4-pulse DEER signal is correct"
        
    assert_experiment_model(ex_ovl4pdeer)
# ======================================================================

def test_ridme1():
# ======================================================================
    "Check that the fit of a S=1/2 RIDME signal is correct"
        
    assert_experiment_model(dl.ex_ridme1)
# ======================================================================

def test_ridme3():
# ======================================================================
    "Check that the fit of a S=3/2 RIDME signal is correct"
        
    assert_experiment_model(dl.ex_ridme3)
# ======================================================================

def test_ridme5():
# ======================================================================
    "Check that the fit of a S=5/2 RIDME signal is correct"
        
    assert_experiment_model(dl.ex_ridme5)
# ======================================================================

def test_dipevo_function():
# ======================================================================
    "Check that the fit of a dipolar evolution function is correct"

    t = np.linspace(0,5,100)
    r = np.linspace(2,7,150)
    P = dd_gauss(r,[4.5, 0.25])
    K = dipolarkernel(t,r)
    V = K@P
    fit = fitmodel(V,t,r,'P',None,None)
    assert ovl(P,fit.P) > 0.90
# ======================================================================

def test_form_factor():
# ======================================================================
    "Check that the fit of a simple form factor is correct"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])
    K = dipolarkernel(t,r,mod=0.3)
    V = K@P
    fit = fitmodel(V,t,r,'P',None,ex_4pdeer)
    assert ovl(P,fit.P) > 0.90
# ======================================================================

def test_full_parametric():
# ======================================================================
    "Check that the fit of full parametric model"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.25])
    lam = 0.3
    kappa = 0.2
    B = bg_exp(t,lam*kappa)
    V = dipolarkernel(t,r,mod=lam,bg=B)@P

    fit = fitmodel(V,t,r,dd_gauss,bg_exp,ex_4pdeer,uq=None)

    assert ovl(P,fit.P) > 0.99
# ======================================================================

def test_no_foreground():
# ======================================================================
    "Check that the fit of a pure background works"

    t = np.linspace(0,5,500)
    r = np.linspace(2,6,200)
    k = 0.2
    B = bg_exp(t,k)

    fit = fitmodel(B,t,r,None,bg_exp,None,uq=None)

    assert abs(k-fit.bgparam) < 1
# ======================================================================

def test_start_values():
# ======================================================================
    "Check that starts values can be correctly specified"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])

    lam = 0.4
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=lam,bg=Bmodel)
    V = K@P

    fit = fitmodel(V,t,r,'P',bg_exp,ex_4pdeer,bg_par0=0.5,ex_par0=0.5,uq=None)

    assert ovl(P,fit.P) > 0.95
# ======================================================================

def test_boundaries():
# ======================================================================
    "Check that boundaries can be correctly specified"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])

    lam = 0.4
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=lam,bg=Bmodel)
    V = K@P 

    fit = fitmodel(V,t,r,'P',bg_exp,ex_4pdeer, uq=None,
                    bg_par0=0.4, bg_lb=0.2, bg_ub=0.5,
                    ex_par0=0.4, ex_lb=0.2, ex_ub=0.5)

    assert ovl(P,fit.P) > 0.95
# ======================================================================

def test_boundaries_adjust_bg():
# ======================================================================
    "Check that start values are adjusted when defining bounds leaving default par0 out"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])

    lam,conc = 0.4,80
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    K = dipolarkernel(t,r,mod=lam,bg=Bmodel)
    V = K@P 

    fit = fitmodel(V,t,r,'P',bg_hom3d,ex_4pdeer, uq=None, bg_lb=70, bg_ub=90)

    assert ovl(P,fit.P) > 0.95
# ======================================================================

def test_boundaries_adjust_ex():
# ======================================================================
    "Check that start values are adjusted when defining bounds leaving default par0 out"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])

    lam,conc = 0.6,80
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    K = dipolarkernel(t,r,mod=lam,bg=Bmodel)
    V = K@P 

    fit = fitmodel(V,t,r,'P',bg_hom3d,ex_4pdeer, uq=None, ex_lb=0.55, ex_ub=0.65)

    assert ovl(P,fit.P) > 0.95
# ======================================================================

def test_boundaries_adjust_dd():
# ======================================================================
    "Check that start values are adjusted when defining bounds leaving default par0 out"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.35])

    lam,conc = 0.4,80
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    K = dipolarkernel(t,r,mod=lam,bg=Bmodel)
    V = K@P 

    fit = fitmodel(V,t,r,dd_gauss,bg_hom3d,ex_4pdeer, uq=None, dd_lb=[4,0.3],dd_ub=[6,0.5])

    assert ovl(P,fit.P) > 0.95
# ======================================================================


def test_global_4pdeer():
# ======================================================================
    "Check the correct fit of two 4-DEER signals"

    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.3])

    parIn = ex_4pdeer.start
    pathways = ex_4pdeer(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)

    t1 = np.linspace(0,5,100)
    V1 = dipolarkernel(t1,r,pathways=pathways,bg=Bmodel)@P
    t2 = np.linspace(0,5,250)
    V2 = dipolarkernel(t2,r,pathways=pathways,bg=Bmodel)@P
    
    fit = fitmodel([V1,V2],[t1,t2],r,'P',bg_exp,ex_4pdeer,uq=None)

    assert ovl(P,fit.P) > 0.90
# ======================================================================

def test_global_full_parametric():
# ======================================================================
    "Check global fitting with fully parametric models"

    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.25])
    lam = 0.3
    kappa = 0.4
    B = lambda t: bg_exp(t,kappa)
    t1 = np.linspace(0,5,100)
    t2 = np.linspace(0,3,200)
    V1 = dipolarkernel(t1,r,mod=lam,bg=B)@P
    V2 = dipolarkernel(t2,r,mod=lam,bg=B)@P

    fit = fitmodel([V1,V2],[t1,t2],r,dd_gauss,bg_exp,ex_4pdeer,uq=None)

    assert ovl(P,fit.P) > 0.99
# ======================================================================

def test_global_mixed_backgrounds():
# ======================================================================
    "Check global fitting with different background models"

    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.25])
    lam = 0.3
    kappa = 0.4
    B = lambda t: bg_exp(t,kappa)
    t1 = np.linspace(0,5,100)
    t2 = np.linspace(0,3,200)
    V1 = dipolarkernel(t1,r,mod=lam,bg=B)@P
    V2 = dipolarkernel(t2,r,mod=lam)@P

    fit = fitmodel([V1,V2],[t1,t2],r,dd_gauss,[bg_exp,None],ex_4pdeer,uq=None)

    assert ovl(P,fit.P) > 0.90
# ======================================================================

def test_global_mixed_experiments():
# ======================================================================
    "Check global fitting with different experiment models"

    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.25])
    lam = 0.3
    t1 = np.linspace(0,5,100)
    t2 = np.linspace(0,3,200)
    V1 = dipolarkernel(t1,r,mod=lam)@P
    V2 = dipolarkernel(t2,r)@P

    fit = fitmodel([V1,V2],[t1,t2],r,dd_gauss,None,[ex_4pdeer,None],uq=None)

    assert ovl(P,fit.P) > 0.9
# ======================================================================


def assert_confidence_intervals(pci50,pci95,pfit,lb,ub):
#----------------------------------------------------------------------\
    p95lb = pci95[:,0]
    p95ub = pci95[:,1]
    p50lb = pci50[:,0]
    p50ub = pci50[:,1]
    p95lb,p95ub,p50lb,p50ub,pfit = [np.round(x,5) for x in (p95lb,p95ub,p50lb,p50ub,pfit)]
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

#----------------------------------------------------------------------
exmodel = ex_4pdeer
ddmodel = dd_gauss
bgmodel = bg_exp

r = np.linspace(2,6,40)
P = ddmodel(r,[4.5, 0.25])

parIn = exmodel.start
pathways = exmodel(parIn)

kappa = 0.4
Bmodel = lambda t,lam: bgmodel(t,kappa)

t = np.linspace(0,5,100)
np.random.seed(0)
V = dipolarkernel(t,r,pathways=pathways,bg=Bmodel)@P + whitegaussnoise(t,0.01)

fit = fitmodel(V,t,r,ddmodel,bgmodel,exmodel,uq='covariance')
#----------------------------------------------------------------------

def assert_confinter_param(subset):
#----------------------------------------------------------------------
    if subset == 'ex':
        model = exmodel
        pfit = fit.exparam
        ci50 = fit.exparamUncert.ci(50)
        ci95 = fit.exparamUncert.ci(95)
    if subset == 'dd':
        model = ddmodel
        pfit = fit.ddparam
        ci50 = fit.ddparamUncert.ci(50)
        ci95 = fit.ddparamUncert.ci(95)
    if subset == 'bg':
        model = bgmodel
        pfit = fit.bgparam
        ci50 = fit.bgparamUncert.ci(50)
        ci95 = fit.bgparamUncert.ci(95)
    lb = model.lower
    ub = model.upper

    assert_confidence_intervals(ci50,ci95,pfit,lb,ub)
#----------------------------------------------------------------------


def test_confinter_exparam():
# ======================================================================
    "Check that the confidence inervals for the experiment parameter are correct"
    assert_confinter_param('ex')
# ======================================================================

def test_confinter_bgparam():
# ======================================================================
    "Check that the confidence inervals for the experiment parameter are correct"
    assert_confinter_param('bg')
# ======================================================================

def test_confinter_ddparam():
# ======================================================================
    "Check that the confidence inervals for the experiment parameter are correct"
    assert_confinter_param('dd')
# ======================================================================

#----------------------------------------------------------------------
exmodel = ex_4pdeer
bgmodel = bg_exp

r = np.linspace(2,6,40)
P = dd_gauss(r,[4.5, 0.25])

parIn = exmodel.start
pathways = exmodel(parIn)

kappa = 0.4
Bmodel = lambda t: bgmodel(t,kappa)

t = np.linspace(0,5,100)
np.random.seed(0)
V = dipolarkernel(t,r,pathways=pathways,bg=Bmodel)@P + whitegaussnoise(t,0.03)

fit_Pparam = fitmodel(V,t,r,ddmodel,bgmodel,exmodel,uq='covariance')
fit_Pfree = fitmodel(V,t,r,'P',bgmodel,exmodel,uq='covariance')
#----------------------------------------------------------------------

def assert_confinter_models(subset):
#----------------------------------------------------------------------

    if subset=='Pfitfree':
        fit = fit_Pfree
    else:
        fit = fit_Pparam

    if subset == 'Pfit' or subset == 'Pfitfree':
        modelfit = fit.P
        lb = np.zeros_like(r)
        ub = np.full_like(r,inf)
        ci50 = fit.Puncert.ci(50)
        ci95 = fit.Puncert.ci(95)
    if subset == 'Bfit':
        modelfit = fit.B
        lb = np.full_like(t,-inf)
        ub = np.full_like(t,inf)
        ci50 = fit.Buncert.ci(50)
        ci95 = fit.Buncert.ci(95)
    if subset == 'Vfit':
        modelfit = fit.V
        lb = np.full_like(t,-inf)
        ub = np.full_like(t,inf)
        ci50 = fit.Vuncert.ci(50)
        ci95 = fit.Vuncert.ci(95)

    assert_confidence_intervals(ci50,ci95,modelfit,lb,ub)
#----------------------------------------------------------------------

def test_confinter_Pfit():
# ======================================================================
    "Check that the confidence inervals for fitted parametric distribution are correct"
    assert_confinter_models('Pfit')
# ======================================================================

def test_confinter_Pfitfree():
# ======================================================================
    "Check that the confidence inervals for fitted distribution are correct"
    assert_confinter_models('Pfitfree')
# ======================================================================

def test_confinter_Vfit():
# ======================================================================
    "Check that the confidence inervals for fitted distribution are correct"
    assert_confinter_models('Vfit')
# ======================================================================

def test_confinter_Bfit():
# ======================================================================
    "Check that the confidence inervals for fitted distribution are correct"
    assert_confinter_models('Bfit')
# ======================================================================


def assert_confinter_noforeground():
# ======================================================================
    "Check that the confidence inervals for a pure background fit are correct"

    bgmodel = bg_exp
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.25])

    kappa = 0.4
    lam = 0.3
    Bmodel = bgmodel(t,kappa)

    np.random.seed(0)
    V = dipolarkernel(t,r,mod=lam,bg=Bmodel)@P + whitegaussnoise(t,0.01)
    
    fit = fitmodel(V,t,r,None,bgmodel,ex_4pdeer,uq='covariance')

    Bfit = fit.B
    Buq = fit.Buncert
    Bci50 = Buq.ci(50)
    Bci95 = Buq.ci(95)
    lb = np.full_like(t,-inf)
    ub = np.full_like(t,inf)

    assert_confidence_intervals(Bci50,Bci95,Bfit,lb,ub)
# ======================================================================

def assert_confinter_dipevofun():
# ======================================================================
    "Check that the confidence inervals for a dipolar evolution function fit are correct"

    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.25])

    t = np.linspace(0,5,100)
    np.random.seed(0)
    V = dipolarkernel(t,r)@P + whitegaussnoise(t,0.01)
    
    fit = fitmodel(V,t,r,'P',None,None,uq='covariance')

    Pfit = fit.P
    Puq = fit.Puncert
    Pci50 = Puq.ci(50)
    Pci95 = Puq.ci(95)
    lb = np.zeros_like(r,0)
    ub = np.full_like(r,inf)

    assert_confidence_intervals(Pci50,Pci95,Pfit,lb,ub)
# ======================================================================

def test_global_scale_4pdeer():
# ======================================================================
    "Check the correct fit of two 4-DEER signals"
    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.25])

    parIn = ex_4pdeer.start
    pathways = ex_4pdeer(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)
    scales = [1e3,1e9]
    t1 = np.linspace(0,5,100)
    V1 = scales[0]*dipolarkernel(t1,r,pathways=pathways,bg=Bmodel)@P
    t2 = np.linspace(0,5,250)
    V2 = scales[1]*dipolarkernel(t2,r,pathways=pathways,bg=Bmodel)@P
    
    fit = fitmodel([V1,V2],[t1,t2],r,'P',bg_exp,ex_4pdeer,uq=None)

    assert max(abs(np.asarray(scales)/np.asarray(fit.scale) - 1)) < 5e-2 
# ======================================================================

def test_global_weights_default():
# ======================================================================
    "Check the correct fit of two 4-DEER signals when one is of very low quality"
    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.25])

    parIn = ex_4pdeer.start
    pathways = ex_4pdeer(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)
    scales = [1e3,1e9]
    t1 = np.linspace(0,5,100)
    V1 = scales[0]*dipolarkernel(t1,r,pathways=pathways,bg=Bmodel)@P + whitegaussnoise(t1,0.001,seed=1)
    t2 = np.linspace(0,5,250)
    V2 = scales[1]*dipolarkernel(t2,r,pathways=pathways,bg=Bmodel)@P + whitegaussnoise(t2,0.1,seed=1)
    
    fit = fitmodel([V1,V2],[t1,t2],r,'P',bg_exp,ex_4pdeer,uq=None)

    assert ovl(P,fit.P)>0.95
# ======================================================================

def test_V_scale_parametric():
# ======================================================================
    "Check that the signal is properly scaled in fully parametric mode"
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.2])
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=0.3,bg=Bmodel)
    scale = 1.54e6
    V = scale*(K@P)

    fit = fitmodel(V,t,r,dd_gauss,bg_exp,ex_4pdeer,uq=None)

    assert isinstance(fit.scale,float) and  abs(1 - scale/fit.scale) < 1e-2
# ======================================================================

def test_V_scale():
# ======================================================================
    "Check that the signal is properly scaled in SNLLS mode"
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.2])
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=0.3,bg=Bmodel)
    scale =1.54e6
    V = scale*(K@P)

    fit = fitmodel(V,t,r,'P',bg_exp,ex_4pdeer,uq=None)

    assert isinstance(fit.scale,float) and abs(1 - scale/fit.scale) < 1e-2
# ======================================================================

def test_V_scale_regularized():
# ======================================================================
    "Check that the signal is properly scaled in regularization mode"
    t = np.linspace(0,5,200)
    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.2])
    K = dipolarkernel(t,r)
    scale = 1.54e6
    V = scale*(K@P)

    fit = fitmodel(V,t,r,'P',None,None,uq=None)

    assert isinstance(fit.scale,float) and  abs(1 - scale/fit.scale) < 1e-2
# ======================================================================

def test_plot():
# ======================================================================
    "Check that the plot method works"
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])

    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=0.4,bg=Bmodel)
    V = K@P

    fit = fitmodel(V,t,r,'P',bg_exp,ex_4pdeer,uq=None)
    
    fig = fit.plot(show=False)
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ======================================================================


def test_physical_bg_model():
# ======================================================================
    "Check that the background parameters of a physical model are fitted correctly"

    t = np.linspace(-0.1,7,200)
    r = np.linspace(3,5,50)
    P = dl.dd_gauss(r,[4,0.2])
    V0 = 3000
    K = dl.dipolarkernel(t,r,mod=0.4,bg=lambda t,lam:dl.bg_hom3d(t,50,lam))
    V = K@P
    V = V0*V

    fit = dl.fitmodel(V,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,uq=None)

    assert abs(fit.bgparam - 50)<1e-1 and abs(fit.exparam - 0.4)<1e-1
# ======================================================================

def test_phenomenological_bg_model():
# ======================================================================
    "Check that the background parameters of a phenomenological model are fitted correctly"

    t = np.linspace(-0.1,7,200)
    r = np.linspace(3,5,50)
    P = dl.dd_gauss(r,[4,0.2])
    V0 = 3000
    K = dl.dipolarkernel(t,r,pathways=dl.ex_4pdeer(0.4),bg=lambda t: dl.bg_exp(t,0.3))
    V = K@P
    V = V0*V

    fit = dl.fitmodel(V,t,r,'P',dl.bg_exp,dl.ex_4pdeer,uq=None)

    assert abs(fit.bgparam - 0.3)<1e-1 and abs(fit.exparam - 0.4)<1e-1
# ======================================================================


def test_Vunmod():
# ======================================================================
    "Check that the background scaling is correct if requested"

    t = np.linspace(-0.1,7,200)
    r = np.linspace(3,5,50)
    P = dl.dd_gauss(r,[4,0.2])
    lam = 0.4
    K = dl.dipolarkernel(t,r,pathways=dl.ex_4pdeer(lam),bg=lambda t,lam: dl.bg_hom3d(t,50,lam))
    Bscaled = (1-lam)*dl.bg_hom3d(t,50,lam)
    V = K@P

    fit = dl.fitmodel(V,t,r,'P',dl.bg_exp,dl.ex_4pdeer,uq=None)

    assert max(abs(Bscaled - fit.Vunmod))<1e-4
# ======================================================================

def test_cost_value():
# ======================================================================
    "Check that starts values can be correctly specified"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,50)
    P = dd_gauss(r,[4.5, 0.25])

    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=0.4,bg=Bmodel)
    V = K@P

    fit = fitmodel(V,t,r,'P',bg_exp,ex_4pdeer,uq=None)

    assert isinstance(fit.cost,float) and np.round(fit.cost/np.sum(fit.residuals**2),5)==1
# ======================================================================

# ----------------------------------------------------------------------
    exmodel = ex_4pdeer
    bgmodel = bg_exp

    r = np.linspace(2,6,40)
    P = dd_gauss(r,[4.5, 0.25])

    parIn = exmodel.start
    pathways = exmodel(parIn)

    kappa = 0.4
    Bmodel = lambda t: bgmodel(t,kappa)

    t = np.linspace(0,5,100)
    np.random.seed(0)
    V = dipolarkernel(t,r,pathways=pathways,bg=Bmodel)@P + whitegaussnoise(t,0.03)

    fit = fitmodel(V,t,r,ddmodel,bgmodel,exmodel,uq=['bootstrap',2])
# ----------------------------------------------------------------------

def assert_boot_ci(quantity):
# ----------------------------------------------------------------------

    if quantity=='P':
        ci = fit.Puncert.ci(95)
    elif quantity=='V':
        ci = fit.Vuncert.ci(95)
    elif quantity=='Vmod':
        ci = fit.VmodUncert.ci(95)
    elif quantity=='Vunmod':
        ci = fit.VunmodUncert.ci(95)
    elif quantity=='B':
        ci = fit.Buncert.ci(95)
    elif quantity=='ddparam':
        ci = fit.ddparamUncert.ci(95)
    elif quantity=='bgparam':
        ci = fit.bgparamUncert.ci(95)
    elif quantity=='exparam':
        ci = fit.exparamUncert.ci(95)

    assert np.all(ci[:,0]<=ci[:,1])
# ----------------------------------------------------------------------


def test_bootci_P():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('P')
# ======================================================================

def test_bootci_V():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('V')
# ======================================================================

def test_bootci_B():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('B')
# ======================================================================

def test_bootci_Vmod():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('Vmod')
# ======================================================================

def test_bootci_Vunmod():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('Vunmod')
# ======================================================================

def test_bootci_ddparam():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('ddparam')
# ======================================================================

def test_bootci_exparam():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('exparam')
# ======================================================================

def test_bootci_bparam():
# ======================================================================
    "Check that the bootstrapped confidence intervals work"
    assert_boot_ci('bgparam')
# ======================================================================

def test_convergence_criteria():
# ======================================================================
    "Check that convergence criteria can be specified without crashing"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,80)
    P = dd_gauss(r,[4.5, 0.25])

    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,mod=0.4,bg=Bmodel)
    V = K@P

    fit = fitmodel(V,t,r,'P',bg_exp,ex_4pdeer,uq=None,tol=1e-3,maxiter=1e2)

    assert ovl(P,fit.P) > 0.90
# ======================================================================

def test_global_weights():
# ======================================================================
    "Check that the global weights properly work when specified"

    t = np.linspace(-0.3,5,300)
    r = np.linspace(2,6,80)

    P1 = dl.dd_gauss(r,[3,0.2])
    P2 = dl.dd_gauss(r,[5,0.2])

    K = dl.dipolarkernel(t,r,mod=0.2)

    scales = [1e3, 1e9]
    sigma1 = 0.001
    V1 = K@P1 + dl.whitegaussnoise(t,sigma1,seed=1)
    sigma2 = 0.001
    V2 = K@P2 + dl.whitegaussnoise(t,sigma2,seed=1)

    V1 = scales[0]*V1
    V2 = scales[1]*V2

    Kmodel= lambda lam: [dl.dipolarkernel(t,r,mod=lam)]*2
    fit1 = dl.fitmodel([V1,V2],[t,t],r,'P',None,dl.ex_4pdeer,weights=[1,0])
    fit2 = dl.fitmodel([V1,V2],[t,t],r,'P',None,dl.ex_4pdeer,weights=[0,1])

    assert ovl(P1,fit1.P) > 0.95 and ovl(P2,fit2.P) > 0.95
# ======================================================================


def test_global_weights_param():
# ======================================================================
    "Check that the global weights properly work when specified"

    t = np.linspace(-0.3,5,300)
    r = np.linspace(2,6,80)

    P1 = dl.dd_gauss(r,[3,0.2])
    P2 = dl.dd_gauss(r,[5,0.2])

    K = dl.dipolarkernel(t,r,mod=0.2)

    scales = [1e3, 1e9]
    sigma1 = 0.001
    V1 = K@P1 + dl.whitegaussnoise(t,sigma1,seed=1)
    sigma2 = 0.001
    V2 = K@P2 + dl.whitegaussnoise(t,sigma2,seed=1)

    V1 = scales[0]*V1
    V2 = scales[1]*V2

    Kmodel= lambda lam: [dl.dipolarkernel(t,r,mod=lam)]*2
    fit1 = dl.fitmodel([V1,V2],[t,t],r,dd_gauss,None,dl.ex_4pdeer,weights=[1,0])
    fit2 = dl.fitmodel([V1,V2],[t,t],r,dd_gauss,None,dl.ex_4pdeer,weights=[0,1])

    assert ovl(P1,fit1.P) > 0.95 and ovl(P2,fit2.P) > 0.95
# ======================================================================

import numpy as np
from numpy import inf
from deerlab import dipolarkernel, whitegaussnoise, fitsignal
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp, bg_hom3d
from deerlab.ex_models import ex_4pdeer, ex_5pdeer, ex_7pdeer, ex_ovl4pdeer
from deerlab.utils import ovl


def assert_experiment_model(model):
# --------------------------------------------------------------------
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.255])

    info = model()
    parIn = info['Start']
    pathways = model(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)
    K = dipolarkernel(t,r,pathways,Bmodel)
    V = K@P

    fit = fitsignal(V,t,r,'P',bg_exp,model,uqanalysis=False)

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

def test_dipevo_function():
# ======================================================================
    "Check that the fit of a dipolar evolution function is correct"

    t = np.linspace(0,5,100)
    r = np.linspace(2,7,150)
    P = dd_gauss(r,[4.5, 0.25])
    K = dipolarkernel(t,r)
    V = K@P
    fit = fitsignal(V,t,r,'P',None,None,uqanalysis=False)
    assert ovl(P,fit.P) > 0.90
# ======================================================================

def test_form_factor():
# ======================================================================
    "Check that the fit of a simple form factor is correct"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])
    K = dipolarkernel(t,r,0.3)
    V = K@P
    fit = fitsignal(V,t,r,'P',None,ex_4pdeer,uqanalysis=False)
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
    V = dipolarkernel(t,r,lam,B)@P

    fit = fitsignal(V,t,r,dd_gauss,bg_exp,ex_4pdeer,uqanalysis=False)

    assert ovl(P,fit.P) > 0.99
# ======================================================================

def test_no_foreground():
# ======================================================================
    "Check that the fit of a pure background works"

    t = np.linspace(0,5,500)
    r = np.linspace(2,6,200)
    conc = 50
    B = bg_hom3d(t,conc)

    fit = fitsignal(B,t,r,None,bg_hom3d,ex_4pdeer,uqanalysis=False)

    assert max(abs(B - fit.B)) < 1e-2 and abs(conc-fit.bgparam*fit.exparam) < 1
# ======================================================================

def test_start_values():
# ======================================================================
    "Check that starts values can be correctly specified"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.25])

    lam = 0.4
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,lam,Bmodel)
    V = K@P

    par0 = [[],0.5,0.5]
    fit = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,par0,uqanalysis=False)

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
    K = dipolarkernel(t,r,lam,Bmodel)
    V = K@P 

    par0 = [[],0.4,0.4]
    lb = [[],0.2,0.2]
    ub = [[],0.5,0.5]
    fit = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,par0,lb,ub,uqanalysis=False)

    assert ovl(P,fit.P) > 0.95
# ======================================================================


def test_global_4pdeer():
# ======================================================================
    "Check the correct fit of two 4-DEER signals"

    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.3])

    info = ex_4pdeer()
    parIn = info['Start']
    pathways = ex_4pdeer(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)

    t1 = np.linspace(0,5,100)
    V1 = dipolarkernel(t1,r,pathways,Bmodel)@P
    t2 = np.linspace(0,5,250)
    V2 = dipolarkernel(t2,r,pathways,Bmodel)@P
    
    fit = fitsignal([V1,V2],[t1,t2],r,'P',bg_exp,ex_4pdeer,uqanalysis=False)

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
    V1 = dipolarkernel(t1,r,lam,B)@P
    V2 = dipolarkernel(t2,r,lam,B)@P

    fit = fitsignal([V1,V2],[t1,t2],r,dd_gauss,bg_exp,ex_4pdeer,uqanalysis=False)

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
    V1 = dipolarkernel(t1,r,lam,B)@P
    V2 = dipolarkernel(t2,r,lam)@P

    fit = fitsignal([V1,V2],[t1,t2],r,dd_gauss,[bg_exp,None],ex_4pdeer,uqanalysis=False)

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
    V1 = dipolarkernel(t1,r,lam)@P
    V2 = dipolarkernel(t2,r)@P

    fit = fitsignal([V1,V2],[t1,t2],r,dd_gauss,None,[ex_4pdeer,None],uqanalysis=False)

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
    assert not errors, "Errors occured:\n{}".format("\n".join(errors))
#----------------------------------------------------------------------

def assert_confinter_param(subset):
#----------------------------------------------------------------------
    exmodel = ex_4pdeer
    ddmodel = dd_gauss
    bgmodel = bg_exp

    r = np.linspace(2,6,40)
    P = ddmodel(r,[4.5, 0.25])

    info = exmodel()
    parIn = info['Start']
    pathways = exmodel(parIn)

    kappa = 0.4
    Bmodel = lambda t,lam: bgmodel(t,kappa)

    t = np.linspace(0,5,100)
    np.random.seed(0)
    V = dipolarkernel(t,r,pathways,Bmodel)@P + whitegaussnoise(t,0.01)
    
    fit = fitsignal(V,t,r,ddmodel,bgmodel,exmodel,uqanalysis=True)

    if subset == 'ex':
        info = exmodel()
        pfit = fit.exparam
        ci50 = fit.exparamUncert.ci(50)
        ci95 = fit.exparamUncert.ci(95)
    if subset == 'dd':
        info = ddmodel()
        pfit = fit.ddparam
        ci50 = fit.ddparamUncert.ci(50)
        ci95 = fit.ddparamUncert.ci(95)
    if subset == 'bg':
        info = bgmodel()
        pfit = fit.bgparam
        ci50 = fit.bgparamUncert.ci(50)
        ci95 = fit.bgparamUncert.ci(95)
    lb = info['Lower']
    ub = info['Upper']

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

def assert_confinter_models(subset):
#----------------------------------------------------------------------
    exmodel = ex_4pdeer
    if subset == 'Pfitfree':
        ddmodel = 'P'
        subset = 'Pfit'
    else:
        ddmodel= dd_gauss
    bgmodel = bg_exp

    r = np.linspace(2,6,40)
    P = dd_gauss(r,[4.5, 0.25])

    info = exmodel()
    parIn = info['Start']
    pathways = exmodel(parIn)

    kappa = 0.4
    Bmodel = lambda t: bgmodel(t,kappa)

    t = np.linspace(0,5,100)
    np.random.seed(0)
    V = dipolarkernel(t,r,pathways,Bmodel)@P + whitegaussnoise(t,0.03)
    
    fit = fitsignal(V,t,r,ddmodel,bgmodel,exmodel,uqanalysis=True)

    if subset == 'Pfit':
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
    V = dipolarkernel(t,r,lam,Bmodel)@P + whitegaussnoise(t,0.01)
    
    fit = fitsignal(V,t,r,None,bgmodel,ex_4pdeer,uqanalysis=True)

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
    
    fit = fitsignal(V,t,r,'P',None,None,uqanalysis=True)

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

    info = ex_4pdeer()
    parIn = info['Start']
    pathways = ex_4pdeer(parIn)

    kappa = 0.4
    Bmodel = lambda t: bg_exp(t,kappa)
    scales = [1e3,1e9]
    t1 = np.linspace(0,5,100)
    V1 = scales[0]*dipolarkernel(t1,r,pathways,Bmodel)@P
    t2 = np.linspace(0,5,250)
    V2 = scales[1]*dipolarkernel(t2,r,pathways,Bmodel)@P
    
    fit = fitsignal([V1,V2],[t1,t2],r,'P',bg_exp,ex_4pdeer,uqanalysis=False)

    assert max(abs(np.asarray(scales)/np.asarray(fit.scale) - 1)) < 1e-2 
# ======================================================================


def test_V_scale_parametric():
# ======================================================================
    "Check that the signal is properly scaled in fully parametric mode"
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.2])
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,0.3,Bmodel)
    scale = 1.54e6
    V = scale*(K@P)

    fit = fitsignal(V,t,r,dd_gauss,bg_exp,ex_4pdeer,uqanalysis=False)

    assert max(abs(1 - V/fit.V)) < 1e-4
# ======================================================================

def test_V_scale():
# ======================================================================
    "Check that the signal is properly scaled in SNLLS mode"
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.2])
    Bmodel = lambda t: bg_exp(t,0.4)
    K = dipolarkernel(t,r,0.3,Bmodel)
    scale =1.54e6
    V = scale*(K@P)

    fit = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,uqanalysis=False)

    assert max(abs(1 - V/fit.V)) < 1e-4
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

    fit = fitsignal(V,t,r,'P',None,None,uqanalysis=False)

    assert max(abs(1 - V/fit.V)) < 1e-4
# ======================================================================

import numpy as np
from numpy import pi, inf, NaN
from deerlab import dipolarkernel, whitegaussnoise, fitsignal
from deerlab.dd_models import dd_gauss, dd_rice
from deerlab.bg_models import bg_exp
from deerlab.ex_models import ex_4pdeer, ex_5pdeer, ex_7pdeer, ex_ovl4pdeer
from deerlab.utils import ovl


def assert_experiment_model(model):
# --------------------------------------------------------------------
    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.6])

    info = model()
    parIn = info['Start']
    pathinfo = model(parIn)

    kappa = 0.4
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)
    K = dipolarkernel(t,r,pathinfo,Bmodel)
    V = K@P

    _,Pfit,_,_,_,_,_ = fitsignal(V,t,r,'P',bg_exp,model)

    assert ovl(P,Pfit) > 0.90
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
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.6])
    K = dipolarkernel(t,r)
    V = K@P
    _,Pfit,_,_,_,_,_ = fitsignal(V,t,r,'P','none','none')
    assert ovl(P,Pfit) > 0.90
# ======================================================================

def test_form_factor():
# ======================================================================
    "Check that the fit of a simple form factor is correct"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.6])
    K = dipolarkernel(t,r,0.3)
    V = K@P
    _,Pfit,_,_,_,_,_ = fitsignal(V,t,r,'P','none',ex_4pdeer)
    assert ovl(P,Pfit) > 0.90
# ======================================================================

def test_full_parametric():
# ======================================================================
    "Check that the fit of full parametric model"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.6])
    lam = 0.3
    kappa = 0.2
    B = bg_exp(t,lam*kappa)
    V = dipolarkernel(t,r,lam,B)@P

    _,Pfit,_,_,_,_,_ = fitsignal(V,t,r,dd_gauss,bg_exp,ex_4pdeer)

    assert ovl(P,Pfit) > 0.99
# ======================================================================

def test_no_foreground():
# ======================================================================
    "Check that the fit of a pure background works"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,200)
    B = bg_exp(t,0.2)

    _,_,Bfit,_,_,_,_ = fitsignal(B,t,r,'none',bg_exp,ex_4pdeer)

    assert max(abs(B - Bfit)) < 1e-2
# ======================================================================

def test_start_values():
# ======================================================================
    "Check that starts values can be correctly specified"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.6])

    lam = 0.4
    Bmodel = lambda t,lam: bg_exp(t,0.4,lam)
    K = dipolarkernel(t,r,lam,Bmodel)
    V = K@P

    par0 = [[],0.5,0.5]
    _,Pfit,_,_,_,_,_ = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,par0)

    assert ovl(P,Pfit) > 0.95
# ======================================================================

def test_boundaries():
# ======================================================================
    "Check that boundaries can be correctly specified"

    t = np.linspace(0,5,100)
    r = np.linspace(2,6,150)
    P = dd_gauss(r,[4.5, 0.6])

    lam = 0.4
    Bmodel = lambda t,lam: bg_exp(t,0.4,lam)
    K = dipolarkernel(t,r,lam,Bmodel)
    V = K@P 

    par0 = [[],0.4,0.4]
    lb = [[],0.2,0.2]
    ub = [[],0.5,0.5]
    _,Pfit,_,_,_,_,_ = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,par0,lb,ub)

    assert ovl(P,Pfit) > 0.95
# ======================================================================


def test_global_4pdeer():
# ======================================================================
    "Check the correct fit of two 4-DEER signals"
    r = np.linspace(2,6,90)
    P = dd_gauss(r,[4.5, 0.6])

    info = ex_4pdeer()
    parIn = info['Start']
    pathinfo = ex_4pdeer(parIn)

    kappa = 0.4
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)

    t1 = np.linspace(0,5,100)
    V1 = dipolarkernel(t1,r,pathinfo,Bmodel)@P
    t2 = np.linspace(0,5,250)
    V2 = dipolarkernel(t2,r,pathinfo,Bmodel)@P
    
    _,Pfit,_,_,_,_,_ = fitsignal([V1,V2],[t1,t2],r,'P',bg_exp,ex_4pdeer)

    assert ovl(P,Pfit) > 0.95
# ======================================================================

def test_global_full_parametric():
# ======================================================================
    "Check global fitting with fully parametric models"

    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.6])
    lam = 0.3
    kappa = 0.4
    B = lambda t,lam: bg_exp(t,lam*kappa,lam)
    t1 = np.linspace(0,5,100)
    t2 = np.linspace(0,3,200)
    V1 = dipolarkernel(t1,r,lam,B)@P
    V2 = dipolarkernel(t2,r,lam,B)@P

    _,Pfit,_,_,_,_,_ = fitsignal([V1,V2],[t1,t2],r,dd_gauss,bg_exp,ex_4pdeer)

    assert ovl(P,Pfit) > 0.99
# ======================================================================

def test_global_mixed_backgrounds():
# ======================================================================
    "Check global fitting with different background models"

    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.6])
    lam = 0.3
    kappa = 0.4
    B = lambda t,lam: bg_exp(t,lam*kappa,lam)
    t1 = np.linspace(0,5,100)
    t2 = np.linspace(0,3,200)
    V1 = dipolarkernel(t1,r,lam,B)@P
    V2 = dipolarkernel(t2,r,lam)@P

    _,Pfit,_,_,_,_,_ = fitsignal([V1,V2],[t1,t2],r,dd_gauss,[bg_exp,'none'],ex_4pdeer)

    assert ovl(P,Pfit) > 0.99
# ======================================================================

def test_global_mixed_experiments():
# ======================================================================
    "Check global fitting with different experiment models"

    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4.5, 0.6])
    lam = 0.3
    t1 = np.linspace(0,5,100)
    t2 = np.linspace(0,3,200)
    V1 = dipolarkernel(t1,r,lam)@P
    V2 = dipolarkernel(t2,r)@P

    _,Pfit,_,_,_,_,_ = fitsignal([V1,V2],[t1,t2],r,dd_gauss,'none',[ex_4pdeer,'none'])

    assert ovl(P,Pfit) > 0.99
# ======================================================================

import numpy as np
from deerlab.bg_models import bg_exp
from deerlab.dipolarbackground import dipolarbackground


def test_basic():
#==================================================================================
    "Basic functionality test on dipolarbackground"

    t = np.linspace(0,5,150)
    kappa = 0.3
    lam = 0.5

    #Reference
    Bref = bg_exp(t,lam*kappa)

    #Output
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)
    path = np.zeros((2,2))
    path[0,:] = [1-lam, np.NaN]
    path[1,:] = [lam, 0]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================


def test_singletime():
#==================================================================================
    "Check that one can generate a background with one single time-domain point"

    t = 2.5
    kappa = 0.3
    lam = 0.5

    #Reference
    Bref = bg_exp(t,lam*kappa)

    #Output
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)
    path = np.zeros((2,2))
    path[0,:] = [1-lam, np.NaN]
    path[1,:] = [lam, 0]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================


def test_harmonics():
#==================================================================================
    "Check that one can generate a background with multiple higher pathway harmonics"

    t = np.linspace(0,5,150)
    kappa = 0.3
    lam = 0.5
    n = 2
    #Reference
    Bref = bg_exp(n*t,lam*kappa)

    #Output
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)
    path = np.zeros((2,3))
    path[0,:] = [1-lam, np.NaN, 1]
    path[1,:] = [lam, 0, n]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================


def test_lambda():
#==================================================================================
    "Check that dipolar background with modulation depth works"

    t = np.linspace(0,5,150)
    kappa = 0.3
    lam = 0.5

    #Reference
    Bref = bg_exp(t,lam*kappa)

    #Output
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)
    B = dipolarbackground(t,lam,Bmodel)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================

def test_multipath_renorm():
#==================================================================================
    "Check that multi-pathway kernels are properly generated without renormalization"

    t1 = np.linspace(0,10,300)
    kappa = 0.3
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = np.array([1-prob, prob**2, prob*(1-prob)])
    T0 = [np.NaN, 0, tau2-t2]
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)

    #Reference
    unmodulated = np.isnan(T0)
    Bref = 1
    Bnorm = 1
    for p in range(len(lam)):
        if not unmodulated[p]:
            Bref = Bref*Bmodel((t-T0[p]),lam[p])
            Bnorm = Bnorm*Bmodel(-T0[p],lam[p])
    
    #Output
    B = dipolarbackground(t,np.array([lam, T0]).T,Bmodel,renormalize=False)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================

def test_multipath_raw():
#==================================================================================
    "Check that multi-pathway kernels are properly generated with normalization"

    t1 = np.linspace(0,10,300)
    kappa = 0.3
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = np.array([1-prob, prob**2, prob*(1-prob)])
    T0 = [np.NaN, 0, tau2-t2]
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)

    #Reference
    unmodulated = np.isnan(T0)
    Bref = 1
    Bnorm = 1
    for p in range(len(lam)):
        if not unmodulated[p]:
            Bref = Bref*Bmodel((t-T0[p]),lam[p])
            Bnorm = Bnorm*Bmodel(-T0[p],lam[p])
    Bref = Bref/Bnorm
    
    #Output
    B = dipolarbackground(t,np.array([lam, T0]).T,Bmodel)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================

def test_overtones():

    t = np.linspace(0,5,150)
    kappa = 0.3
    lam = 0.5
    overtones = [0.6, 0.3, 0.1]
    # Reference
    Bref = 1
    for i in range(len(overtones)):
        n = i+1
        Bref = Bref*bg_exp(t,n*overtones[i]*lam*kappa)
    
    # Output
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)
    path = np.zeros((2,2))
    path[0,:] = [1-lam, np.NaN]
    path[1,:] = [lam, 0]

    B = dipolarbackground(t,path,Bmodel,overtonecoeff=overtones)

    assert max(abs(B-Bref)) < 1e-8


import numpy as np
from deerlab.bg_models import bg_hom3d, bg_exp
from deerlab.dipolarbackground import dipolarbackground


def test_basic():
#==================================================================================
    "Basic functionality test on dipolarbackground"

    t = np.linspace(0,5,150)
    conc = 50
    lam = 0.5

    #Reference
    Bref = bg_hom3d(t,lam*conc)

    #Output
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    path = [[],[]]
    path[0] = [1-lam]
    path[1] = [lam, 0]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================


def test_singletime():
#==================================================================================
    "Check that one can generate a background with one single time-domain point"

    t = 2.5
    conc = 50
    lam = 0.5

    #Reference
    Bref = bg_hom3d(t,lam*conc)

    #Output
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    path = [[],[]]
    path[0] = [1-lam]
    path[1] = [lam, 0]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================


def test_harmonics():
#==================================================================================
    "Check that one can generate a background with multiple higher pathway harmonics"

    t = np.linspace(0,5,150)
    conc = 50
    lam = 0.5
    n = 2
    #Reference
    Bref = bg_hom3d(n*t,lam*conc)

    #Output
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    path = [[],[]]
    path[0] = [1-lam]
    path[1] = [lam, 0,2]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================

def test_multipath_renorm():
#==================================================================================
    "Check that multi-pathway kernels are properly generated without renormalization"

    t1 = np.linspace(0,10,300)
    conc = 50
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = [prob**2, prob*(1-prob)]
    T0 = [0, tau2-t2]
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)

    # Reference
    Bref = 1
    for p in range(len(lam)):
            Bref = Bref*Bmodel((t-T0[p]),lam[p])
    
    paths = []
    paths.append(1-prob)
    paths.append([prob**2,0])
    paths.append([prob*(1-prob),tau2-t2])

    # Output
    B = dipolarbackground(t,paths,Bmodel)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================

def test_multipath_raw():
#==================================================================================
    "Check that multi-pathway kernels are properly generated with normalization"

    t1 = np.linspace(0,10,300)
    conc = 50
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = [prob**2, prob*(1-prob)]
    T0 = [0, tau2-t2]
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)

    #Reference
    Bref = 1
    for p in range(len(lam)):
            Bref = Bref*Bmodel((t-T0[p]),lam[p])
    
    paths = []
    paths.append(1-prob)
    paths.append([prob**2,0])
    paths.append([prob*(1-prob),tau2-t2])

    #Output
    B = dipolarbackground(t,paths,Bmodel)

    assert max(abs(B-Bref)) < 1e-8
#==================================================================================

def test_physical():
#==================================================================================
    "Check that physical background models are propely used"

    t = np.linspace(0,5,150)
    conc = 50
    lam = 0.5

    #Reference
    Bref = bg_hom3d(t,conc,lam)

    #Output
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    path = [[],[]]
    path[0] = [1-lam]
    path[1] = [lam, 0]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================

def test_phenomenological():
#==================================================================================
    "Check that phenomenological background models are propely used"

    t = np.linspace(0,5,150)
    kappa = 0.3
    lam = 0.5

    #Reference
    Bref = bg_exp(t,0.3)

    #Output
    Bmodel = lambda t: bg_hom3d(t,kappa)
    path = [[],[]]
    path[0] = [1-lam]
    path[1] = [lam, 0]
    B = dipolarbackground(t,path,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================

import numpy as np
from deerlab.bg_models import bg_hom3d, bg_exp, bg_homfractal
from deerlab.dipolarbackground import dipolarbackground
from deerlab.utils import assert_docstring

def test_basic():
#==================================================================================
    "Basic functionality test on dipolarbackground"

    t = np.linspace(0,5,150)
    conc = 50
    lam = 0.5

    #Reference
    Bref = bg_hom3d(t,conc,lam)

    #Output
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================


def test_singletime():
#==================================================================================
    "Check that one can generate a background with one single time-domain point"

    t = 2.5
    conc = 50
    lam = 0.5

    #Reference
    Bref = bg_hom3d(t,conc,lam)

    #Output
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,Bmodel)


    assert max(abs(B-Bref) < 1e-8)
#==================================================================================


def test_fractalharmonics():
#==================================================================================
    "Check that one can generate a background with higher pathway harmonics"

    t = np.linspace(0, 5, 150)

    conc = 500
    lam = 0.423
    delta = 2
    tref = 0.14
    dim = 2.7

    # Reference
    Bref = bg_homfractal(delta*(t-tref), conc, dim, lam)

    # Output
    Bmodel = lambda t, lam: bg_homfractal(t, conc, dim, lam)
    pathways = [{'amp': 1-lam},
                {'reftime': tref, 'amp': lam, 'harmonic': delta}]
    B = dipolarbackground(t,pathways,Bmodel)

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

    Bref = 1
    for p in range(len(lam)):
        Bref = Bref*Bmodel((t-T0[p]),lam[p])
    
    pathways = [{'amp': (1-prob)**2},
                {'reftime': 0, 'amp': prob**2},
                {'reftime': tau2-t2, 'amp': prob*(1-prob)}]
    B = dipolarbackground(t,pathways,Bmodel)

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

    Bref = 1
    for p in range(len(lam)):
            Bref = Bref*Bmodel((t-T0[p]),lam[p])
    
    pathways = [{'amp': (1-prob)**2},
                {'reftime': 0, 'amp': prob**2},
                {'reftime': tau2-t2, 'amp': prob*(1-prob)}]
    B = dipolarbackground(t,pathways,Bmodel)

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
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,Bmodel)

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
    Bmodel = lambda t: bg_hom3d(t,kappa,1)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,Bmodel)

    assert max(abs(B-Bref) < 1e-8)
#==================================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(dipolarbackground)
# ======================================================================

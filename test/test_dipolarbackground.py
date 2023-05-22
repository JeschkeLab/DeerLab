
import numpy as np
from deerlab.bg_models import _hom3d, bg_exp, _homfractal
from deerlab.dipolarbackground import dipolarbackground
from deerlab.utils import assert_docstring
import pytest 

# Fixtures 
# -----------------------------------------------------------------------------------
t = np.linspace(0,5,150)
t1 = t
t2 = np.linspace(0,3,150)
conc = 50
kappa = 0.3
lam = 0.5
tref1 = 0.14
delta = 1
dim = 2.7
lam1,lam2,lam3 = lam,0.2,0.3

@pytest.fixture(scope='module')
def BmodelA():
    return lambda t,lam: _hom3d(t,conc,lam)

@pytest.fixture(scope='module')
def BmodelB():
    return lambda t,lam: _homfractal(t, conc, dim, lam)
# -----------------------------------------------------------------------------------


# ==================================================================================
def test_basic(BmodelA):
    "Basic functionality test on dipolarbackground"
    Bref = BmodelA(t,lam)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,BmodelA)
    assert np.allclose(B,Bref)
# ==================================================================================

# ==================================================================================
def test_singletime(BmodelA):
    "Check that one can generate a background with one single time-domain point"
    t = 2
    Bref = BmodelA(t,lam)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,BmodelA)
    assert np.allclose(B,Bref)
# ==================================================================================

# ==================================================================================
def test_fractalharmonics(BmodelB):
    "Check that one can generate a background with higher pathway harmonics"
    Bref = BmodelB(delta*(t-tref1), lam)
    pathways = [{'amp': 1-lam},
                {'reftime': tref1, 'amp': lam, 'harmonic': delta}]
    B = dipolarbackground(t,pathways,BmodelB)
    assert np.allclose(B,Bref)
# ==================================================================================

# ==================================================================================
def test_multipath_renorm(BmodelA):
    "Check that multi-pathway kernels are properly generated without renormalization"
    Bref = 1
    for tref,lam in zip([0,tref1], [lam2,lam3]):
        Bref = Bref*BmodelA((t-tref),lam)
    pathways = [{'amp': lam1},
                {'reftime': 0, 'amp': lam2},
                {'reftime': tref, 'amp': lam3}]
    B = dipolarbackground(t,pathways,BmodelA)
    assert np.allclose(B,Bref)
# ==================================================================================

# ==================================================================================
def test_physical(BmodelA):
    "Check that physical background models are propely used"
    Bref = BmodelA(t,lam)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,BmodelA)
    assert np.allclose(B,Bref)
# ==================================================================================

# ==================================================================================
def test_phenomenological():
    "Check that phenomenological background models are propely used"
    Bref = bg_exp(t,kappa)
    pathways = [{'amp': 1-lam},
                {'reftime': 0, 'amp': lam}]
    B = dipolarbackground(t,pathways,lambda t: bg_exp(t,kappa))
    assert np.allclose(B,Bref)
# ==================================================================================

def test_twodimensional(BmodelA):
#==================================================================================
    "Check the construction of a two-pathway bidimensional background function"
    Bref = _hom3d(t1[:,None],conc,lam)*_hom3d(t2[None,:],conc,lam)
    pathways = [{'amp': 1-lam},
                {'reftime': [0,0], 'amp': lam, 'harmonic': [1,0]},
                {'reftime': [0,0], 'amp': lam, 'harmonic': [0,1]},]
    B = dipolarbackground([t1,t2],pathways,BmodelA)
    assert np.allclose(B,Bref)
#==================================================================================

# ==================================================================================
def test_threespin(BmodelA):
    "Check the construction of a three-spin background contribution"
    Bref = _hom3d(t,conc,lam1)*_hom3d(t,conc,lam2)*_hom3d(t,conc,lam3)
    pathways = [{'amp': 1-lam1-lam2-lam3},
                {'reftime': (0,None,None), 'amp': lam1, 'harmonic': (1,0,0)},
                {'reftime': (None,0,None), 'amp': lam2, 'harmonic': (0,1,0)},
                {'reftime': (None,None,0), 'amp': lam3, 'harmonic': (0,0,1)}]
    B = dipolarbackground(t,pathways,BmodelA)
    assert np.allclose(B,Bref)
# ==================================================================================

# ==================================================================================
def test_threespin_twodimensional(BmodelA):
    "Check the construction of a three-spin background contribution"
    Bref = _hom3d(t1[:,None],conc,lam1)*_hom3d(t2[None,:],conc,lam2)*_hom3d(t1[:,None],conc,lam3)
    pathways = [{'amp': 1-lam1-lam2-lam3},
                {'reftime': ([0,0],[None,None],[None,None]), 'amp': lam1, 'harmonic': ([1,0],[0,0],[0,0])},
                {'reftime': ([None,None],[0,0],[None,None]), 'amp': lam2, 'harmonic': ([0,0],[0,1],[0,0])},
                {'reftime': ([None,None],[None,None],[0,0]), 'amp': lam3, 'harmonic': ([0,0],[0,0],[1,0])}]
    B = dipolarbackground([t1,t2],pathways,BmodelA)
    assert np.allclose(B,Bref)
# ==================================================================================

# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(dipolarbackground)
# ======================================================================

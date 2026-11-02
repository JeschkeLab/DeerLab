
import numpy as np
from math import factorial
from deerlab.bg_models import _hom3d, bg_exp, _homfractal, bg_hom3d
from deerlab.dipolarbackground import dipolarbackground, dipolarbackgroundmodel
from deerlab.dipolarmodel import ex_3pdeer, ex_4pdeer
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


# ==================================================================================
# Tests for dipolarbackgroundmodel
# ==================================================================================

tau = 0.5
tau1, tau2 = 0.5, 3.0

# ==================================================================================
def test_dipolarbackgroundmodel_twospin_single_pathway():
    "Test dipolarbackgroundmodel for a two-spin system with a single pathway"
    mod = 0.4
    tref = tau
    Bref = (1-mod) * _hom3d(t-tref, conc, mod)
    experiment = ex_3pdeer(tau, pathways=[1])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d)
    Bval = Bmodel(t, conc, mod, tref)
    assert np.allclose(Bval, Bref)
# ==================================================================================

# ==================================================================================
def test_dipolarbackgroundmodel_twospin_single_pathway_no_experiment():
    "Test dipolarbackgroundmodel for a two-spin system with a single pathway with no experiment"
    mod = 0.4
    tref = tau
    Bref = (1-mod) * _hom3d(t-tref, conc, mod)
    Bmodel = dipolarbackgroundmodel(None, bg_hom3d)
    Bval = Bmodel(t, conc, mod, tref)
    assert np.allclose(Bval, Bref)
# ==================================================================================


# ==================================================================================
def test_dipolarbackgroundmodel_twospin_multi_pathway():
    "Test dipolarbackgroundmodel for a two-spin system with multiple pathways"
    lam1, lam2 = 0.3, 0.1
    tref1, tref2 = tau1, tau1+tau2
    Bref = (1-lam1-lam2) * _hom3d(t-tref1, conc, lam1) * _hom3d(t-tref2, conc, lam2)
    experiment = ex_4pdeer(tau1, tau2, pathways=[1,2])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d)
    Bval = Bmodel(t, conc, lam1, lam2, tref1, tref2)
    assert np.allclose(Bval, Bref)
# ==================================================================================

# ==================================================================================
def test_dipolarbackgroundmodel_twospin_delays():
    "Test dipolarbackgroundmodel for a two-spin system with delays parametrization"
    lam1, lam2 = 0.3, 0.1
    tref1, tref2 = tau1, tau1+tau2
    Bref = (1-lam1-lam2) * _hom3d(t-tref1, conc, lam1) * _hom3d(t-tref2, conc, lam2)
    experiment = ex_4pdeer(tau1, tau2, pathways=[1,2])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d, parametrization='delays')
    Bval = Bmodel(t, conc, tau1, tau2, lam1, lam2)
    assert np.allclose(Bval, Bref)
# ==================================================================================

# ==================================================================================
def _multispin_reference_pathways(lams, trefs, harmonics, lamu, Q=3):
    """
    Reference multi-spin pathways following eqn (15) of Fábregas-Ibáñez et al., PCCP 24 (2022) 22645,
    with lams[q][p] the amplitude of pair pathway p on interaction q (three-spin system, all interactions joined)
    """
    Λ0 = 1.0
    paths, threespin = [], []
    P = len(trefs)
    # Two-spin pathways: (Q-1)! λu^(Q-1) λp, pair pathway p on interaction q
    for p in range(P):
        for q in range(Q):
            λ2k = factorial(Q-1)*lams[q][p]*lamu**(Q-1)
            paths.append({'amp': λ2k,
                          'reftime': tuple(trefs[p] if n==q else None for n in range(Q)),
                          'harmonic': tuple(harmonics[p] if n==q else 0 for n in range(Q))})
            Λ0 -= λ2k
    # Three-spin pathways: 2(Q-2)! λu^(Q-2) λp λp', pair pathways p<=p' on interactions q!=q'
    for p in range(P):
        for p2 in range(p, P):
            for q1 in range(Q):
                for q2 in range(Q):
                    if q1==q2:
                        continue
                    λ3k = 2*factorial(Q-2)*lams[q1][p]*lams[q2][p2]*lamu**(Q-2)
                    threespin.append({'amp': λ3k,
                                      'reftime': tuple(trefs[p] if n==q1 else trefs[p2] if n==q2 else None for n in range(Q)),
                                      'harmonic': tuple(harmonics[p] if n==q1 else harmonics[p2] if n==q2 else 0 for n in range(Q))})
                    Λ0 -= λ3k
    paths.append({'amp': Λ0})
    return paths + threespin
# ==================================================================================

# ==================================================================================
def test_dipolarbackgroundmodel_threespin_single_pathway():
    "Test dipolarbackgroundmodel for a three-spin system with a single pathway"
    lam, lamu = 0.4, 0.8
    paths_ref = _multispin_reference_pathways([[lam]]*3, [tau], [1], lamu)
    # A single pathway also gives three-spin pathways (same pair pathway on two interactions)
    assert len(paths_ref) == 3 + 1 + 6
    Bref = dipolarbackground(t, paths_ref, lambda td, l: _hom3d(td, conc, l))
    experiment = ex_3pdeer(tau, pathways=[1])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d, spins=3)
    Bval = Bmodel(t, conc, lam, tau, lamu)
    assert np.allclose(Bval, Bref)
# ==================================================================================

# ==================================================================================
def test_dipolarbackgroundmodel_threespin_multi_pathway():
    "Test dipolarbackgroundmodel for a three-spin system with multiple pathways"
    lam1, lam2, lamu = 0.3, 0.1, 0.7
    paths_ref = _multispin_reference_pathways([[lam1,lam2]]*3, [tau1, tau1+tau2], [1, 1], lamu)
    Bref = dipolarbackground(t, paths_ref, lambda td, l: _hom3d(td, conc, l))
    experiment = ex_4pdeer(tau1, tau2, pathways=[1,2])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d, spins=3)
    Bval = Bmodel(t, conc, lam1, lam2, tau1, tau1+tau2, lamu)
    assert np.allclose(Bval, Bref)
# ==================================================================================

# ==================================================================================
def test_dipolarbackgroundmodel_threespin_nosamespins():
    "Test dipolarbackgroundmodel for a three-spin system with samespins=False"
    lam_q = [0.3, 0.4, 0.2]   # different amplitudes per interaction
    lamu = 0.8
    paths_ref = _multispin_reference_pathways([[lam] for lam in lam_q], [tau], [1], lamu)
    Bref = dipolarbackground(t, paths_ref, lambda td, l: _hom3d(td, conc, l))

    experiment = ex_3pdeer(tau, pathways=[1])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d, spins=3, samespins=False)
    Bval = Bmodel(t, conc, *lam_q, tau, lamu)
    assert np.allclose(Bval, Bref)
# ==================================================================================

# ==================================================================================
def test_dipolarbackgroundmodel_threespin_delays():
    "Test dipolarbackgroundmodel for a three-spin system with delays parametrization"
    lam, lamu = 0.4, 0.8
    paths_ref = _multispin_reference_pathways([[lam]]*3, [tau], [1], lamu)
    Bref = dipolarbackground(t, paths_ref, lambda td, l: _hom3d(td, conc, l))

    experiment = ex_3pdeer(tau, pathways=[1])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d, spins=3, parametrization='delays')
    # signature: t(const), conc, delay1, lam1, lamu
    Bval = Bmodel(t, conc, tau, lam, lamu)
    assert np.allclose(Bval, Bref)
# ==================================================================================

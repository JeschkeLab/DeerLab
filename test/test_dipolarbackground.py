
import numpy as np
from itertools import permutations
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
def test_dipolarbackgroundmodel_threespin_single_pathway():
    "Test dipolarbackgroundmodel for a three-spin system with a single pathway"
    lam, lamu = 0.4, 0.8
    Q = 3
    # Construct reference pathways manually
    Λ0 = 1.0
    paths_ref = []
    for perm in set(permutations([0]+[None]*(Q-1))):
        λ2k = factorial(Q-1)*lam*lamu**(Q-1)
        tref_path = tuple([tau if n==0 else None for n in perm])
        δs_path = tuple([1 if n==0 else 0 for n in perm])
        paths_ref.append({'amp': λ2k, 'reftime': tref_path, 'harmonic': δs_path})
        Λ0 -= λ2k
    paths_ref.append({'amp': Λ0})
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
    Q = 3
    harmonics_exp = [1, 1]
    trefs_exp = [tau1, tau1+tau2]

    # Construct reference pathways manually
    Λ0 = 1.0
    paths_ref, threespin_ref = [], []
    for idx, (tref_idx, δ_idx) in enumerate(zip(trefs_exp, harmonics_exp)):
        lam_idx = [lam1, lam2][idx]
        for perm in set(permutations([0]+[None]*(Q-1))):
            q = int(np.where(np.array(perm)==0)[0][0])
            λp = lam_idx
            λ2k = factorial(Q-1)*λp*lamu**(Q-1)
            paths_ref.append({'amp': λ2k,
                               'reftime': tuple([tref_idx if n==0 else None for n in perm]),
                               'harmonic': tuple([δ_idx if n==0 else 0 for n in perm])})
            Λ0 -= λ2k
            for idx2 in range(idx+1, 2):
                tref_idx2 = trefs_exp[idx2]
                δ_idx2 = harmonics_exp[idx2]
                for perm2 in set(permutations([0,1]+[None]*(Q-2))):
                    perm2_arr = np.array(perm2)
                    q2 = int(np.where(perm2_arr==1)[0][0])
                    λm = lam_idx
                    λ3k = 2*factorial(Q-2)*λp*λm*lamu**(Q-2)
                    threespin_ref.append({'amp': λ3k,
                                          'reftime': tuple([tref_idx if n==0 else tref_idx2 if n==1 else None for n in perm2]),
                                          'harmonic': tuple([δ_idx if n==0 else δ_idx2 if n==1 else 0 for n in perm2])})
                    Λ0 -= λ3k
    paths_ref.append({'amp': Λ0})
    Bref = dipolarbackground(t, paths_ref+threespin_ref, lambda td, l: _hom3d(td, conc, l))

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
    Q = 3
    # Construct reference pathways with per-interaction amplitudes
    Λ0 = 1.0
    paths_ref = []
    for perm in set(permutations([0]+[None]*(Q-1))):
        q = int(np.where(np.array(perm)==0)[0][0])
        λ2k = factorial(Q-1)*lam_q[q]*lamu**(Q-1)
        tref_path = tuple([tau if n==0 else None for n in perm])
        δs_path = tuple([1 if n==0 else 0 for n in perm])
        paths_ref.append({'amp': λ2k, 'reftime': tref_path, 'harmonic': δs_path})
        Λ0 -= λ2k
    paths_ref.append({'amp': Λ0})
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
    Q = 3
    Λ0 = 1.0
    paths_ref = []
    for perm in set(permutations([0]+[None]*(Q-1))):
        λ2k = factorial(Q-1)*lam*lamu**(Q-1)
        tref_path = tuple([tau if n==0 else None for n in perm])
        δs_path = tuple([1 if n==0 else 0 for n in perm])
        paths_ref.append({'amp': λ2k, 'reftime': tref_path, 'harmonic': δs_path})
        Λ0 -= λ2k
    paths_ref.append({'amp': Λ0})
    Bref = dipolarbackground(t, paths_ref, lambda td, l: _hom3d(td, conc, l))

    experiment = ex_3pdeer(tau, pathways=[1])
    Bmodel = dipolarbackgroundmodel(experiment, bg_hom3d, spins=3, parametrization='delays')
    # signature: t(const), conc, delay1, lam1, lamu
    Bval = Bmodel(t, conc, tau, lam, lamu)
    assert np.allclose(Bval, Bref)
# ==================================================================================

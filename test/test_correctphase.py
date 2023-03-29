import numpy as np
from numpy import pi
from deerlab import correctphase
from deerlab.utils import assert_docstring

# Fixtures
dummy_data = np.arange(100)
dummy_phase = 4/7*pi
dummy_phased_data = dummy_data*np.exp(-1j*dummy_phase)
    
def test_basics():
#============================================================
    "Check the basic behavior correcting a phased signal"
    corrected_data = correctphase(dummy_phased_data)
    assert max(np.real(dummy_data) - np.real(corrected_data)) < 1e-4
#============================================================

def test_phase_fit():
#============================================================
    "Check that the phase shift is correctly fitted"
    _,_,corrected_phase = correctphase(dummy_phased_data,full_output=True)
    assert abs(dummy_phase - corrected_phase) < 1e-4
#============================================================

def test_multiple_datasets():
#============================================================
    "Check that the phase correction works when passing multiple datasets"
    dummy_datasets = np.tile(np.arange(100),(20,1)).T
    dummy_phased_datasets = dummy_datasets*np.exp(-1j*np.mod(np.linspace(-3*pi/4,pi/2,20),pi))
    corrected_datasets = correctphase(dummy_phased_datasets)
    assert np.max(np.real(corrected_datasets) - np.real(dummy_datasets)) < 1e-4
#============================================================

def test_posrealsum():
#============================================================
    "Check that real part of the phase signal has a positive average"
    corrected_data = correctphase(dummy_phased_data)
    assert max(np.real(dummy_data) - np.real(corrected_data)) < 1e-4
#============================================================

def test_offset():
#============================================================
    "Check the basic behavior correcting a phased signal"
    corrected_data = correctphase(dummy_phased_data + 1j*500,offset=True)
    assert max(np.real(dummy_data) - np.real(corrected_data)) < 1e-4
#============================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(correctphase)
# ======================================================================


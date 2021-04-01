import numpy as np
from numpy import pi
from deerlab import correctphase

def test_basics():
#============================================================
    "Check the basic behaviour correcting a phased signal"

    V = np.arange(100)
    Vphased = V*np.exp(-1j*4/7*pi)
    Vcorr = correctphase(Vphased)

    assert max(np.real(V) - np.real(Vcorr)) < 1e-4
#============================================================

def test_phase_fit():
#============================================================
    "Check that the phase shift is correctly fitted"

    V = np.arange(100)
    phase = 4/7*pi
    Vphased = V*np.exp(-1j*phase)
    _,_,phase_fit = correctphase(Vphased,full_output=True)

    assert abs(phase - phase_fit) < 1e-4
#============================================================

def test_multiple_datasets():
#============================================================
    "Check that the phase correction works when passing multiple datasets"

    V = np.tile(np.arange(100),(20,1)).T
    phases = np.mod(np.linspace(-3*pi/4,pi/2,20),pi)
    Vphased = V*np.exp(-1j*phases)

    Vcorr = correctphase(Vphased)
    
    assert np.max(np.real(V) - np.real(Vcorr)) < 1e-4
#============================================================

def test_mode_posrealint():
#============================================================
    "Check the basic behaviour correcting a phased signal with maximal real integral"

    V = np.arange(100)
    Vphased = V*np.exp(-1j*4/7*pi)
    Vcorr = correctphase(Vphased,phase='posrealint')

    assert max(np.real(V) - np.real(Vcorr)) < 1e-4
#============================================================

def test_mode_negrealint():
#============================================================
    "Check the basic behaviour correcting a phased signal with maximal negative real integral"

    V = np.arange(100)
    Vphased = V*np.exp(-1j*8/7*pi)
    Vcorr = correctphase(Vphased,phase='negrealint')

    assert max(np.real(V) - np.real(-Vcorr)) < 1e-4
#============================================================

def test_mode_close():
#============================================================
    "Check the basic behaviour correcting a phased signal with phase closest to average"

    V = np.arange(100)
    Vphased = V*np.exp(-1j*4/7*pi)
    Vcorr = correctphase(Vphased,phase='close')

    assert max(np.real(V) - np.real(Vcorr)) < 1e-4
#============================================================
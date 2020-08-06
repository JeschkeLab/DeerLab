import numpy as np
from numpy import pi
from deerlab import correctphase


def test_basics():
#============================================================
    "Check the basic behaviour correcting a phased signal"

    V = np.arange(100)
    Vphased = V*np.exp(1j*4/7*pi)
    Vcorr = correctphase(Vphased)

    assert max(abs(V - Vcorr)) < 1e-4
#============================================================


def test_phase_fit():
#============================================================
    "Check that the phase shift is correctly fitted"

    V = np.arange(100)
    phase = 4/7*pi
    Vphased = V*np.exp(1j*phase)
    _,_,phase_fit,_ = correctphase(Vphased,full_output=True)

    assert abs(phase - phase_fit) < 1e-4
#============================================================


def test_phase_manual():
#============================================================
    "Check that the phase shift can be correctly specified manually" 

    V = np.arange(100)
    phase = 4/7*pi
    Vphased = V*np.exp(1j*phase)
    Vcorr = correctphase(Vphased,phase)

    assert max(abs(V - Vcorr)) < 1e-4
#============================================================


def test_fit_imag_offset():
#============================================================
    "Check the fitting of imaginary offsets"

    V = np.arange(100)
    phiIn = pi/4
    ImOffsetIn = 40j
    Vphased = V*np.exp(1j*phiIn) + ImOffsetIn

    _,_,_,ImOffsetOut = correctphase(Vphased,[],imagoffset=True,full_output=True)

    assert abs(ImOffsetIn - ImOffsetOut) < 1e-4
#============================================================


def test_multiple_datasets():
#============================================================
    "Check that the phase correcion works when passing multiple datasets"

    V = np.tile(np.arange(100),(20,1)).T
    phases = np.mod(np.linspace(-3*pi/4,pi/2,20),pi)
    Vphased = V*np.exp(1j*phases)

    Vcorr = correctphase(Vphased)
    
    assert np.max(np.abs(Vcorr - V)) < 1e-4
#============================================================


def test_multiple_datasets_manual():
#============================================================
    "Check that manual phase correcion works when passing multiple datasets"

    V = np.tile(np.arange(100),(20,1)).T
    phases = np.mod(np.linspace(-3*pi/4,pi/2,20),pi)
    Vphased = V*np.exp(1j*phases)

    Vcorr = correctphase(Vphased,phases)
    
    assert np.max(np.abs(Vcorr - V)) < 1e-4
#============================================================
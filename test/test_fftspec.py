import numpy as np
from deerlab import fftspec


def test_basic():
# ======================================================================
    "Check that the magnitude spectrum can be computed from the real/imag spectra"

    t = np.linspace(0,1/np.mean(2*10)*100,100)
    S = np.exp(-t)*np.cos(2*np.pi*5*t)

    specRef = abs(np.fft.fftshift(np.fft.fft(S,2*len(S))))
    _,spec = fftspec(S,t,apodization=False)


    assert max(abs(specRef - spec)) < 1e-10
# ======================================================================

def test_modes():
# ======================================================================
    "Check that the magnitude spectrum can be computed from the real/imag spectra"

    t = np.linspace(0,1/np.mean(2*10)*100,100)
    S = np.exp(-t)*np.cos(2*np.pi*5*t)

    _,specAbs = fftspec(S,t,mode='abs')
    _,specRe = fftspec(S,t,mode='real')
    _,specIm = fftspec(S,t,mode='imag')

    assert max(np.abs(np.sqrt(specRe**2 + specIm**2) - specAbs)) < 1e-10
# ======================================================================

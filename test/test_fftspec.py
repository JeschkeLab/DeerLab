import numpy as np
from deerlab import fftspec
from deerlab.utils import assert_docstring
from pytest import fixture

# Fixtures 
# -----------------------------------------------------------------------
@fixture(scope='module')
def mock_axis():
    return np.linspace(0,1/np.mean(2*10)*100,100)

@fixture(scope='module')
def mock_data(mock_axis):
    return np.exp(-mock_axis)*np.cos(2*np.pi*5*mock_axis)
# -----------------------------------------------------------------------
    
def test_basic(mock_axis,mock_data):
# ======================================================================
    "Check that the magnitude spectrum can be computed from the real/imag spectra"

    specref = abs(np.fft.fftshift(np.fft.fft(mock_data,2*len(mock_data))))
    _,spec = fftspec(mock_data,mock_axis,apodization=False)

    assert max(abs(specref - spec)) < 1e-10
# ======================================================================

def test_modes(mock_axis,mock_data):
# ======================================================================
    "Check that the magnitude spectrum can be computed from the real/imag spectra"

    _,specAbs = fftspec(mock_data,mock_axis,mode='abs')
    _,specRe = fftspec(mock_data,mock_axis,mode='real')
    _,specIm = fftspec(mock_data,mock_axis,mode='imag')

    assert max(np.abs(np.sqrt(specRe**2 + specIm**2) - specAbs)) < 1e-10
# ======================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(fftspec)
# ======================================================================

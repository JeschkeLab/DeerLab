from deerlab.utils.utils import assert_docstring
import numpy as np
import pytest
from deerlab import whitegaussnoise
from deerlab.utils import assert_docstring

def test_length():
# ======================================================================
    "Check that the array length is correct"

    N = 151
    t = np.linspace(0,3,N)
    noise = whitegaussnoise(t,0.2)

    assert len(noise) == N
# ======================================================================


def test_rescale():
# ======================================================================
    "Check that the rescale keyword is working"

    N = 10
    t = np.linspace(0,3,N)
    sig = 1.1234
    noise = whitegaussnoise(t,sig,rescale=True)
    
    assert np.isclose(np.std(noise),sig)
# ======================================================================


def test_seed():
# ======================================================================
    "Check that the seed keyword is working"

    N = 151
    t = np.linspace(0,3,N)
    sig = 0.2
    s = 432
    noise1 = whitegaussnoise(t,sig,seed=s)
    noise2 = whitegaussnoise(t,sig,seed=s)
    
    assert np.array_equal(noise1,noise2)
# ======================================================================


def test_noseed():
# ======================================================================
    "Check that the seed keyword is working if set to None"

    N = 151
    t = np.linspace(0,3,N)
    sig = 0.2
    noise1 = whitegaussnoise(t,sig,seed=None)
    noise2 = whitegaussnoise(t,sig,seed=None)
    
    assert not np.array_equal(noise1,noise2)
# ======================================================================


def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(whitegaussnoise)
# ======================================================================


def test_requiredstd():
# ======================================================================
    "Check that the standard deviation is a required argument"

    N = 10
    t = np.linspace(0,3,N)
    
    with pytest.raises(TypeError):
        noise = whitegaussnoise(t)
    
# ======================================================================



from deerlab.utils import assert_docstring
import numpy as np
import pytest
from deerlab import whitegaussnoise
from deerlab.utils import assert_docstring

# Fixtures 
# ----------------------------------------------------------------------
N = 151
t = np.linspace(0,3,N)
noiselvl = 0.234
# ----------------------------------------------------------------------

def test_length():
# ======================================================================
    "Check that the array length is correct"
    noise = whitegaussnoise(t,0.2)
    assert len(noise) == N
# ======================================================================

def test_rescale():
# ======================================================================
    "Check that the rescale keyword is working"
    noise = whitegaussnoise(t,noiselvl,rescale=True)
    assert np.isclose(np.std(noise),noiselvl)
# ======================================================================

def test_seed():
# ======================================================================
    "Check that the seed keyword is working"
    noise1 = whitegaussnoise(t,noiselvl,seed=432)
    noise2 = whitegaussnoise(t,noiselvl,seed=432)
    assert np.array_equal(noise1,noise2)
# ======================================================================

def test_noseed():
# ======================================================================
    "Check that the seed keyword is working if set to None"
    noise1 = whitegaussnoise(t,noiselvl,seed=None)
    noise2 = whitegaussnoise(t,noiselvl,seed=None)
    assert not np.array_equal(noise1,noise2)
# ======================================================================

def test_requiredstd():
# ======================================================================
    "Check that the standard deviation is a required argument"
    with pytest.raises(TypeError):
        noise = whitegaussnoise(t)
# ======================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(whitegaussnoise)
# ======================================================================


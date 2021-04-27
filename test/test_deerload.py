import os
import pytest
import numpy as np 
from deerlab import deerload

#==================================================================================
@pytest.fixture
def rootdir():
    """ Pytest fixture to be able to define data files relative to this test's location """
    return os.path.dirname(os.path.abspath(__file__))
#==================================================================================


def test_1D_load(rootdir):
#==================================================================================
    """ Check that a 1-dimensional dataset if properly loded """

    t,V = deerload(os.path.join(rootdir, 'data/data_1d.DTA'))
    
    assert isinstance(t,np.ndarray) and isinstance(V,np.ndarray) and len(t)==len(V)
#==================================================================================


def test_2D_load(rootdir):
#==================================================================================
    """ Check that a 2-dimensional dataset if properly loded """

    t,V = deerload(os.path.join(rootdir, 'data/data_2d.DTA'))
    
    assert len(t[0]) == np.shape(V)[0] and len(t[1]) == np.shape(V)[1]
#==================================================================================


def test_sparse(rootdir):
#==================================================================================
    """ Check that a sparsely non-uniformly sampled dataset can be properly loded """

    t,V = deerload(os.path.join(rootdir, 'data/data_sparse.DTA'))
    
    assert isinstance(t,np.ndarray) and isinstance(V,np.ndarray) and len(t)==len(V)
#==================================================================================

def test_multiple_labs(rootdir):
#==================================================================================
    """ Check that data from different labs/spectrometers can be loaded succesfully """

    successfull = True
    for lab in np.arange(1,6,1):
        print(lab)
        t,V = deerload(os.path.join(rootdir, f'data/data_multi_lab{lab}.DTA'))
    
        successfull = successfull and (isinstance(t,np.ndarray) and isinstance(V,np.ndarray))
        successfull = successfull and len(t)==len(V)

    assert np.all(successfull)
#==================================================================================





import numpy as np
from dipolarkernel import dipolarkernel

def test_matrixsize_fresnel():
    "Check non-square kernel matrix construction"
    Ntime = 50
    Ndist = 70
    t = np.linspace(-0.5,5,50)
    r = np.linspace(2,6,70)
    K = dipolarkernel(t,r)
    assert K.shape == (Ntime, Ndist)

def test_singletime():
    "Check that one can generate a kernel with one single time-domain point"

    t = 0.5
    r = np.linspace(1,5,100)
    K = dipolarkernel(t,r)

    assert np.shape(K)[0]==1
    
def test_singledist():
    "Check that one can generate a kernel with one single distance-domain point"

    t = np.linspace(0,5,50)
    r = 3.5
    K = dipolarkernel(t,r)

    assert np.shape(K)[1]==1

def test_negative_time_fresnel():
    "Check that kernel is constructed properly for negative times using fresnel method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = np.linspace(2,6,70)
    Kneg = dipolarkernel(tneg,r)
    Kpos = dipolarkernel(tpos,r)

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)


def test_value_fresnel():
    "Test whether kernel matrix element (calculated using Fresnel integrals) is correct."

    # Generate kernel numerically
    t = 1 # us
    r = 1 # nm
    K = dipolarkernel(t,r)

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-14
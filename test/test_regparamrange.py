
import numpy as np
from deerlab import dipolarkernel, regoperator, regparamrange

def test_resolution():
#=======================================================================
    "Check resolution can be set up manually"
    
    t = np.linspace(0,3,100)
    r = np.linspace(2,6,100)
    K = dipolarkernel(t,r)
    L = regoperator(r,2)

    alpha = regparamrange(K,L,logres=0.1)

    assert len(alpha)==84
#=======================================================================


def test_values():
#=======================================================================
    "Check that returned values are correct"
    
    t = np.linspace(0,3,100)
    r = np.linspace(2,6,5)
    K = dipolarkernel(t,r,integralop=False)
    L = regoperator(r,2)

    alpharef = [0.39810, 0.50119,  0.63096,  0.79432,  1.00000, 1.25893,
                1.58489, 1.99526,  2.51189,  3.16228,  3.98107, 5.01187,
                6.30957, 7.94328, 10.00000, 12.58925, 15.84893]
    alpha = regparamrange(K,L)
    print(alpha, alpharef)
    assert np.all(abs(alpha - alpharef) < 1e-4)
#=======================================================================

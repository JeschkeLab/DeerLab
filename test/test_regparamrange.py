
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

    # Reference: using MATLAB gsvd()
    alpharef = [0.3981,0.5012,0.6310,0.7943,1.0000,1.2589,1.5849,1.9953]
    alpha = regparamrange(K,L)
    print(alpha, alpharef)
    assert np.all(abs(alpha - alpharef) < 1e-4)
#=======================================================================

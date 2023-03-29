import numpy as np
from numpy import pi
from deerlab import regoperator
from deerlab.utils import assert_docstring
import pytest 

# Fixtures 
# ----------------------------------------------------------------------
n = 100
r = np.linspace(2, 6, n)
# ----------------------------------------------------------------------

@pytest.mark.parametrize('includeedges',[True,False])
@pytest.mark.parametrize('order',[0,1,2,3])
def test_regoperator_shape(order,includeedges):
#=======================================================================
    "Check that regularization operator is returned with correct size"
    L = regoperator(r,order,includeedges=includeedges)
    if includeedges:
        expected_shape = (n,n)
    else:   
        expected_shape = (n-order, n)
    assert L.shape == expected_shape
#=======================================================================

# ======================================================================
@pytest.mark.parametrize('type',['uniform','nonuniform'])
@pytest.mark.parametrize('order',[0,1,2,3])
def test_numerical_vs_analytical_derivative(order,type):
    "Check that the numerical derivatives obtained with the operators match the analytical ones accurately"
    if type == 'nonuniform':
        r = np.linspace(0,4.5**9,2500)**(1/9)
    else:
        r = np.linspace(2,5,2500)
    dr = np.mean(np.diff(r))

    # Distance axis for the derivative
    if order==0:
        rn = r
    elif order==1:
        rn = r - dr/2
    elif order==2:
        rn = r 
    elif order==3:
        rn = r - dr/2

    # Test function
    sig = 0.2
    rmean = 3.5
    P = 1/sig/np.sqrt(2*pi)*np.exp(-(r - rmean)**2/2/sig**2)

    # Analytical derivatives of the Gaussian function
    if order==0:
        dPnref = P
    elif order==1:
        dPnref = -(rn - rmean)/sig**3/np.sqrt(2*pi)*np.exp(-(rn - rmean)**2/2/sig**2)
    elif order==2:
        dPnref = -(sig**2-(rn - rmean)**2)/sig**5/np.sqrt(2*pi)*np.exp(-(rn - rmean)**2/2/sig**2)
    elif order==3:
        dPnref = -((rn - rmean)**3 - 3*((rn - rmean))*sig**2)/sig**7/np.sqrt(2*pi)*np.exp(-(rn - rmean)**2/2/sig**2)

    # Numerical finite-difference operators
    Ln = regoperator(r, order, includeedges=True)

    # Estimated numerical derivative of the Gaussian function
    dPn = Ln@P

    if type == 'nonuniform':
        dPn = dPn[:-1]
        dPnref = dPnref[:-1]

        if order==0:
            accuracy = 1e-5
        elif order==1:
            accuracy = 5e-2
        elif order==2:
            accuracy = 5e-1
        elif order==3:
            accuracy = 5e-0
    else:
        accuracy = 10**(-5+order)
        
    assert max(abs(dPn - dPnref)) < accuracy
# ======================================================================

def test_edge_smoothing():
#=======================================================================
    "Check that L applies to edges and is NxN"
    r = np.arange(5)
    n = len(r)
    L = regoperator(r, 2, includeedges=True)
    Lref = (np.eye(n, n)*(-2) + np.eye(n, n, k=-1) + np.eye(n, n, k=1))
    assert np.max(np.abs(L - Lref)) < 1e-10
#=======================================================================

def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(regoperator)
# ======================================================================

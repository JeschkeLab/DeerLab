
import numpy as np
from numpy import pi
from deerlab import regoperator

def test_L0shape():
#=======================================================================
    "Check that L0 is returned with correct size"
    
    n = 100
    r = np.linspace(2,6,n)
    L = regoperator(r,0)

    assert L.shape==(n, n)
#=======================================================================

def test_L1shape():
#=======================================================================
    "Check that L1 is returned with correct size"
    
    n = 100
    r = np.linspace(2,6,n)
    L = regoperator(r,1)

    assert L.shape==(n, n)
#=======================================================================

def test_L2shape():
#=======================================================================
    "Check that L2 is returned with correct size"
    
    n = 100
    r = np.linspace(2,6,n)
    L = regoperator(r,2)

    assert L.shape==(n, n)
#=======================================================================

def test_edge_smoothing():
#=======================================================================
    "Check that L applies to edges and is NxN"

    r = np.arange(5)
    n = len(r)
    L = regoperator(r,2)
    Lref = (np.eye(n, n)*(-2) + np.eye(n, n, k=-1) + np.eye(n, n, k=1))

    assert np.max(np.abs(L - Lref)) < 1e-10
#=======================================================================

def compare_analytical_derivative(n,nonuniform=False):
#=======================================================================

    if nonuniform:
        r = np.linspace(0,4.5**9,2500)**(1/9)
    else:
        r = np.linspace(2,5,2500)
    dr = np.mean(np.diff(r))

    # Distance axis for the derivative
    if n==0:
        rn = r
    elif n==1:
        rn = r - dr/2
    elif n==2:
        rn = r 
    elif n==3:
        rn = r - dr/2

    # Test function
    sig = 0.2
    rmean = 3.5
    P = 1/sig/np.sqrt(2*pi)*np.exp(-(r - rmean)**2/2/sig**2)

    # Analytical derivatives of the Gaussian function
    if n==0:
        dPnref = P
    elif n==1:
        dPnref = -(rn - rmean)/sig**3/np.sqrt(2*pi)*np.exp(-(rn - rmean)**2/2/sig**2)
    elif n==2:
        dPnref = -(sig**2-(rn - rmean)**2)/sig**5/np.sqrt(2*pi)*np.exp(-(rn - rmean)**2/2/sig**2)
    elif n==3:
        dPnref = -((rn - rmean)**3 - 3*((rn - rmean))*sig**2)/sig**7/np.sqrt(2*pi)*np.exp(-(rn - rmean)**2/2/sig**2)

    # Numerical finitie difference operators
    Ln = regoperator(r,n)

    # Estimated derivatives of the Gaussian function
    dPn = Ln@P

    if nonuniform:
        dPn = dPn[:-1]
        dPnref = dPnref[:-1]

    if nonuniform and n==0:
        accuracy = 1e-5
    elif nonuniform and n==1:
        accuracy = 5e-2
    elif nonuniform and n==2:
        accuracy = 5e-1
    elif nonuniform and n==3:
        accuracy = 5e-0
    else:
        accuracy = 10**(-5+n)
        
    assert max(abs(dPn - dPnref)) < accuracy

def test_0th_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 0th derivative of a Gaussian on a uniform grid
    """
    
    compare_analytical_derivative(0)
#=======================================================================

def test_1st_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 1st derivative of a Gaussian on a uniform grid
    """
    
    compare_analytical_derivative(1)
#=======================================================================

def test_2nd_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 2nd derivative of a Gaussian on a uniform grid
    """
    
    compare_analytical_derivative(2)
#=======================================================================

def test_3rd_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 3rd derivative of a Gaussian on a uniform grid
    """
    
    compare_analytical_derivative(3)
#=======================================================================

def test_nonuniform_0th_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 0th derivative of a Gaussian on a non-uniform grid
    """
    
    compare_analytical_derivative(0,nonuniform=True)
#=======================================================================

def test_nonuniform_1st_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 1st derivative of a Gaussian on a non-uniform grid
    """
    
    compare_analytical_derivative(1,nonuniform=True)
#=======================================================================

def test_nonuniform_2nd_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 2nd derivative of a Gaussian on a non-uniform grid
    """
    
    compare_analytical_derivative(2,nonuniform=True)
#=======================================================================

def test_nonuniform_3rd_derivative():
#=======================================================================
    """
    Check that returned finite difference matrix correctly estimates 
    the 3rd derivative of a Gaussian on a non-uniform grid
    """
    
    compare_analytical_derivative(3,nonuniform=True)
#=======================================================================
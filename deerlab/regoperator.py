
# regoperator.py - Regularization Operators
# --------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np

def regoperator(r,d=2):
    r""" Computes the discrete approximation to the derivative operators used as regularization operators.

    Parameters
    ----------
    r : array_like with shape(n,)
        Distance axis, in nanometers.
    d : int scalar
        Derivative order, the default is 2.

    Returns
    -------
    L : ndarray with shape(n-2,n)
        Regularization operator (or finite-difference derivative operator).
        
    Notes
    -----
    All finite-difference matrix elements are computed via Fornberg's method [1]_. Therefore, the 
    resulting operator is agnostic with respect to the increments in the distance vector, allowing 
    uniformly as well as non-uniformly spaced distance vectors. 

    References
    ----------
    .. [1] B. Fornberg, "Calculation of weights in finite difference formulas", SIAM Review 40 (1998), pp. 685-691.
    """
    r = np.atleast_1d(r)

    # Length of axis
    n = len(r)

    # Construct finite difference matrix
    L = np.zeros((n-d,n))

    # Compute non-zero finite forward difference coefficients via Fornberg's method
    for i in range(n-d):
        cols = np.arange(i,i+d+1,1)
        L[i,cols] = _fdcoeffF(d,r[i],r[cols])

    # Introduce missing rows to account for edges of axis
    for __ in range(int(np.ceil(d/2))):
        L = np.concatenate([np.atleast_2d(np.append(L[0,1:],0)), L])
    for __ in range(int(np.floor(d/2))):
        L = np.concatenate([L, np.atleast_2d(np.insert(L[-1,:-1],0,0))])

    return L


def _fdcoeffF(k,xbar,x):
    """
    Compute coefficients for finite difference approximation for the
    derivative of order k at xbar based on grid values at points in x.

    This function returns a row vector c oAf dimension 1 by n, where n=length(x),
    containing coefficients to approximate u^{(k)}(xbar),
    the k'th derivative of u evaluated at xbar,  based on n values
    of u at x(1), x(2), ... x(n).

    If U is a column vector containing u(x) at these n points, then
    c*U will give the approximation to u^{(k)}(xbar).

    Note for k=0 this can be used to evaluate the interpolating polynomial
    itself.

    Requires length(x) > k.
    Usually the elements x(i) are monotonically increasing
    and x(1) <= xbar <= x(n), but neither condition is required.
    The x values need not be equally spaced but must be distinct.

    This program should give the same results as fdcoeffV.m, but for large
    values of n is much more stable numerically.

    Based on the program "weights" in
    B. Fornberg, "Calculation of weights in finite difference formulas",
    SIAM Review 40 (1998), pp. 685-691.

    Note: Forberg's algorithm can be used to simultaneously compute the
    coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
    This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
    the coefficients for the k'th derivative.

    In this version we set m=k and only compute the coefficients for
    derivatives of order up to order k, and then return only the k'th column
    of the resulting C matrix (converted to a row vector).
    This routine is then compatible with fdcoeffV.
    It can be easily modified to return the whole array if desired.

    From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)
    """

    n = len(x)
    if k >= n:
        raise TypeError('Numer of elements in x must be larger than k')


    m = k
    c1 = 1
    c4 = x[0] - xbar
    C = np.zeros((n,m+1))
    C[0,0] = 1
    for i in range(n-1):
        i1 = i+1
        mn = min(i,m)
        c2 = 1
        c5 = c4
        c4 = x[i1] - xbar
        for j in range(i+1):
            j1 = j
            c3 = x[i1] - x[j1]
            c2 = c2*c3
            if j==i:
                for s in range(mn+1,0,-1):
                    s1 = s
                    C[i1,s1] = c1*(s*C[i1-1,s1-1] - c5*C[i1-1,s1])/c2
                C[i1,0] = -c1*c5*C[i1-1,0]/c2
            for s in range(mn+1,0,-1):
                s1 = s
                C[j1,s1] = (c4*C[j1,s1] - s*C[j1,s1-1])/c3
            C[j1,0] = c4*C[j1,0]/c3
        c1 = c2
    # Last column of c gives desired row vector
    c = C[:,-1]           
    return c

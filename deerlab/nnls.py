# nnls.py - Non-linear least squares solvers
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import cvxopt as cvx
import math as m

def fnnls(AtA,Atb,tol=[],maxiter=[],verbose=False):
#=====================================================================================
    r"""
    FNNLS   Fast non-negative least-squares algorithm.
    x = fnnls(AtA,Atb) solves the problem min ||b - Ax|| if
        AtA = A'*A and Atb = A'*b.
    A default tolerance of TOL = MAX(SIZE(AtA)) * NORM(AtA,1) * EPS
    is used for deciding when elements of x are less than zero.
    This can be overridden with x = fnnls(AtA,Atb,TOL).

    [x,w] = fnnls(AtA,Atb) also returns dual vector w where
        w(i) < 0 where x(i) = 0 and w(i) = 0 where x(i) > 0.
    
    For the FNNLS algorithm, see
        R. Bro, S. De Jong
        A Fast Non-Negativity-Constrained Least Squares Algorithm
        Journal of Chemometrics 11 (1997) 393-401
    The algorithm FNNLS is based on is from
        Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.
    """

    unsolvable = False
    count = 0

    # Use all-zero starting vector
    N = np.shape(AtA)[1]
    x = np.zeros(N)

    # Calculate tolerance and maxiter if not given.
    if np.size(np.atleast_1d(tol))==0:
        eps = np.finfo(float).eps
        tol = 10*eps*np.linalg.norm(AtA,1)*max(np.shape(AtA))
    if np.size(np.atleast_1d(maxiter))==0:
        maxiter = 5*N


    passive = x>0       # initial positive/passive set (points where constraint is not active)
    x[~passive] = 0
    w = Atb - AtA @ x     # negative gradient of error functional 0.5*||A*x-y||^2
    
    # Outer loop: Add variables to positive set if w indicates that fit can be improved.
    outIteration = 0
    maxIterations = 5*N    
    while np.any(w>tol) and np.any(~passive):
        outIteration += 1
        
        # Add the most promising variable (with largest w) to positive set.
        t = np.argmax(w)
        passive[t] = True
        
        # Solve unconstrained problem for new augmented positive set.
        # This gives a candidate solution with potentially new negative variables.
        x_ =np.zeros(N)
        if np.sum(passive)==1:
            x_[passive] = Atb[passive]/AtA[passive,passive]
        else:
            x_[passive] = np.linalg.solve(AtA[np.ix_(passive,passive)],Atb[passive])
        
        # Inner loop: Iteratively eliminate negative variables from candidate solution.
        iIteration = 0
        while any((x_<=tol) & passive) and iIteration<maxIterations:
            iIteration += 1
            
            # Calculate maximum feasible step size and do step.
            negative = (x_<=tol) & passive
            alpha = min(x[negative]/(x[negative]-x_[negative]))
            x += alpha*(x_-x)
            
            # Remove all negative variables from positive set.
            passive[x<tol] = False
            
            # Solve unconstrained problem for reduced positive set.
            x_ = np.zeros(N)
            if np.sum(passive)==1:
                x_[passive] = Atb[passive]/AtA[passive,passive]
            else:
                x_[passive] = np.linalg.solve(AtA[np.ix_(passive,passive)],Atb[passive])
            
        # Accept non-negative candidate solution and calculate w.
        if all(x == x_):
            count += 1
        else:
            count = 0
        if count > 5:
            unsolvable = True
            break
        x = x_
        
        w = Atb - AtA@x
        w[passive] = -m.inf
        if verbose:
            print('{:10.0f}{:15.0f}{:20.4e}\n'.format(outIteration,iIteration,max(w)) )

    if verbose:
        if unsolvable:
            print('Optimization stopped because the solution cannot be further changed. \n')
        elif any(~passive):
            print('Optimization stopped because the active set has been completely emptied. \n')
        elif w>tol:
            print('Optimization stopped because the gradient (w) is inferior than the tolerance value TolFun = #.6e. \n' %tol)
        else:
            print('Solution found. \n')
    
    return x

#=====================================================================================




def cvxnnls(AtA, Atb, tol=[], maxiter=[]):
#=====================================================================================
    """
    NNLS problem solved via CVXOPT
    
        
    References:
    -----------
    [1] 
    Rein, Lewe, Andrade, Kacprzak and Weber. J. Magn. Reson., 295 (2018) 17–26.
    Global analysis of complex PELDOR time traces
    https://doi.org/10.1016/j.jmr.2018.07.015
    """

    N = np.shape(AtA)[1]
    if np.size(np.atleast_1d(tol))==0:
        eps = np.finfo(float).eps
        tol = 10*eps*np.linalg.norm(AtA,1)*max(np.shape(AtA))
    if np.size(np.atleast_1d(maxiter)):
        maxiter = 5*N

    x0 = np.zeros(N)

    cAtA = cvx.matrix(AtA)
    cAtb = -cvx.matrix(Atb)
       
    lb = cvx.matrix(np.zeros(N))
    I = -cvx.matrix(np.eye(N, N))
    
    # Set optimization stop criteria
    cvx.solvers.options['show_progress'] = False
    cvx.solvers.options['max_iters'] = maxiter
    cvx.solvers.options['abstol'] = tol
    cvx.solvers.options['reltol'] = tol

    P = cvx.solvers.qp(cAtA, cAtb, I, lb, initvals=cvx.matrix(x0))['x']
    P = np.squeeze(np.asarray(P))
    return P
#=====================================================================================


def nnlsbpp(AtA,AtB,x0=None):
#=====================================================================================
    """
    Non-Negative Least Squares using Block Principal Pivoting
    ==========================================================

    Usage:
    ------
        x = nnlsbpp(AtA,AtB,x0)

    Arguments: 
    ----------
    AtA (NxN-element matrix)

    AtB (Nx1-element array or Nxk matrix)
    
    x0  (Nx1 vector or Nxk matrix) 
        Initial guess(es) vector
    
    Returns:
    --------
    x (Nx1-element array or Nxk matrix)
        Solution vector(s)

    References:
    -----------
    [1] 
    Portugal, Judice, Vicente, Mathematics of Computation, 1994, 63, 625-643
    A comparison of block pivoting and interior-point algorithms for linear
    least squares problems with nonnegative variables
    https://doi.org/10.1090/S0025-5718-1994-1250776-4
    
    [2]
    Kim, Park,SIAM J. Sci. Comput. 2011, 33(6), 3261-3281
    Fast Nonnegative Matrix Factorization: An Active-Set-Like Method and Comparisons
    https://doi.org/10.1137/110821172
    """

    # Size checks
    #-------------------------------------------------------------------------------
    n1,n2 = np.shape(AtA)

    if n1 is not n2:
        raise TypeError('AtA must be a square matrix. You gave a ',n1,'x',n2,' matrix.')
    n = np.size(AtB)
    k = 1
    if n is not n1:
        raise TypeError('AtB must have the same number of rows as AtA. You gave ',n,' instead of ',n1)

    if x0 is None:
        x0 = np.linalg.solve(AtA,AtB)

    # Loop over multiple right-hand sides
    #-------------------------------------------------------------------------------
    if k > 1:
        x = np.zeros((n1,k))
        for k_ in reversed(range(k)):
            x[:,k_] = nnlsbpp(AtA,AtB[:,k_],x0[:,k_])
        return  x

    # Calculate initial solution
    #-------------------------------------------------------------------------------
    x = np.zeros(n)
    if not np.all(x0):
        Fset = np.full((n,1),False)
        y = -AtB
    else:
        Fset = x0 > 0
        x = np.zeros(n)
        x[Fset] = np.linalg.solve(AtA[np.ix_(Fset,Fset)],AtB[Fset])
        y = AtA@x - AtB
 
    # Determine infeasible variables ( = variables with negative values)
    xFnegative = (x < 0) & Fset
    yGnegative = (y < 0) & ~Fset
    nInfeasible = sum(xFnegative | yGnegative)

    # Iterative algorithm
    #-------------------------------------------------------------------------------
    p = 3
    t = np.inf
    while nInfeasible > 0: # iterate until no infeasible variables left
    
        # Swap full blocks in/out of F set or swap a single element as backup plan.
        if nInfeasible < t:
            t = nInfeasible
            p = 3
            Fset[xFnegative] = False
            Fset[yGnegative] = True
        else:
            if p >= 1:
                p = p - 1
                Fset[xFnegative] = False
                Fset[yGnegative] = True
            else:
                idx_ = np.where(xFnegative | yGnegative)[-1]
                Fset[idx_] = ~Fset[idx_]
        
        # Solve linear system over F set
        x = np.zeros(n)
        x[Fset] = np.linalg.solve(AtA[np.ix_(Fset,Fset)],AtB[Fset])
        y = AtA@x - AtB
        
        # Determine infeasible variables
        xFnegative = (x < 0) & Fset
        yGnegative = (y < 0) & ~Fset
        nInfeasible = sum(xFnegative | yGnegative)

    return x
#=====================================================================================

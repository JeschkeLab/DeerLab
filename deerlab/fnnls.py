import numpy as np
import cvxopt as cvo
import math as m

def fnnls(AtA,Atb,tol=-1,verbose=False):
    """
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
    x = np.zeros(np.shape(AtA)[1])

    # Calculate tolerance if not given.
    if tol==-1:
        eps = 2**-52
        tol = 10*eps*np.linalg.norm(AtA,1)*max(np.shape(AtA))
    
    N = np.size(x)
    
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
            print("%10i%15i%20.4e\n" % (outIteration,iIteration,max(w)))

    if verbose:
        if unsolvable:
            print('Optimization stopped because the solution cannot be further changed. \n')
        elif any(~passive):
            print('Optimization stopped because the active set has been completely emptied. \n')
        elif w>tol:
            print('Optimization stopped because the gradient (w) is inferior than the tolerance value TolFun = %.6e. \n' % tol)
        else:
            print('Solution found. \n')
    
    return x

def cvxnnls(AtA, Atb, K, V):
    """
    CVXOPT Based NNLS solver modified from Rein et al. 2018 (GloPel) 
    """
    m, n = K.shape
    
    P = np.linalg.inv(AtA) @ K.T @ V
    P = P.clip(min=0)

    cAtA = cvo.matrix(AtA)
    cAtb = -cvo.matrix(Atb)

    lb = cvo.matrix(np.zeros(n))
    I = -cvo.matrix(np.eye(n, n))
    
    cvo.solvers.options['show_progress'] = False
    cvo.solvers.options['abstol'] = 1e-16
    cvo.solvers.options['reltol'] = 1e-15 

    P = cvo.solvers.qp(cAtA, cAtb, I, lb, initvals=cvo.matrix(P))['x']
    P = np.squeeze(np.asarray(P))
    return P

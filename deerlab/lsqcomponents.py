import numpy as np

def lsqcomponents(V, K, L=None, alpha=0, weights=None, regtype='tikhonov', huberparam=1.35):
# ==============================================================================================
    """
    Calculate the components needed for the least-squares (LSQ) solvers
    """
    
    if weights is None:
        weights = np.ones_like(V)
    else:
        weights = np.atleast_1d(weights)
        
    # Compute components of the LSQ normal equations
    KtK = K.T@(weights[:,np.newaxis]*K)
    KtV = K.T@(weights*V)
    
    # No regularization term -> done
    if L is None:
        return KtK, KtV
    
    def optimizeregterm(regfun, threshold, maxIter):
    # ============================================================================
        nr = K.shape[1]
        P = np.zeros(nr)
        for _ in range(maxIter):
            Pprev = P
            # Compute punconstrained distribution iteratively
            KtKreg_ = KtK + alpha**2*regfun(P)
            P = P + np.linalg.solve(KtKreg_,KtV)
            
            # Normalize distribution by its integral to stabilize convergence
            P = P/np.sum(abs(P))
            # Monitor largest changes in the distribution
            change = np.max(abs(P - Pprev))
            # Stop if result is stable
            if change < threshold:
                break
        # Compute the regularization term from the optimized result
        regterm = regfun(P)
        return regterm
    # ============================================================================
    
    # Compute the regularization term
    if regtype.lower() == 'tikhonov':
        regterm = L.T@L

    elif regtype.lower() == 'tv':
        maxIter = 500
        changeThreshold = 1e-1
        TVFcn = lambda _P: L.T@((L/np.sqrt((L@_P)**2 + np.finfo(float).eps)[:,np.newaxis]))
        regterm = optimizeregterm(TVFcn,changeThreshold,maxIter)
        
    elif regtype.lower() == 'huber':
        maxIter = 500
        changeThreshold = 1e-2
        HuberFcn = lambda _P: 1/(huberparam**2)*(L.T@(L/np.sqrt((L@_P/huberparam)**2 + 1)[:,np.newaxis]))
        regterm = optimizeregterm(HuberFcn,changeThreshold,maxIter)
        
    else:
        raise ValueError(f"Regularization type {regtype} not recognized. Must be 'tikhonov', 'tv' or 'huber'")
    
    KtKreg = KtK + alpha**2*regterm
    
    return KtKreg, KtV
# ==============================================================================================

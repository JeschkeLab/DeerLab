import numpy as np

def lsqcomponents(V,K,L,alpha,weights, regtype = 'tikh',huberparam = 1.35):
# ==============================================================================================

    # Prepare
    nr = np.shape(K)[1]

    # Compute residual terms
    KtV = weights*K.T@V
    KtK = weights*K.T@K

    def optimizeregterm(fun,threshold):
    # ============================================================================
        P = np.zeros(nr)
        for _ in range(maxIter):
            Pprev = P
            # Compute pseudoinverse and unconst. distribution recursively
            KtKreg_ = KtK + alpha**2*fun(P)
            P = P + np.linalg.solve(KtKreg_,weights*K.T@V)
            
            # Normalize distribution by its integral to stabilize convergence
            P = P/np.sum(abs(P))
            # Monitor largest changes in the distribution
            change = np.max(abs(P - Pprev))
            # Stop if result is stable
            if change < threshold:
                break
        # Compute the regularization term from the optimized result
        regterm = fun(P)
        return regterm
    # ============================================================================

    # Compute then the LSQ components needed by NNLS optimizers    
    if regtype.lower() == 'tikh' or regtype.lower() == 'tikhonov':
        regterm = L.T@L
        
    elif regtype.lower() == 'tv':
        maxIter = 500
        changeThreshold = 1e-1
        TVFcn = lambda p: L.T@((L/np.sqrt((L@p)**2 + np.finfo(float).eps)[:,np.newaxis]))
        regterm = optimizeregterm(TVFcn,changeThreshold)
        
    elif regtype.lower() == 'huber':
        maxIter = 500
        changeThreshold = 1e-2
        HuberFcn = lambda p: 1/(huberparam**2)*(L.T@(L/np.sqrt((L@p/huberparam)**2 + 1)[:,np.newaxis]))
        regterm = optimizeregterm(HuberFcn,changeThreshold)
    else:
        raise ValueError("Regularization type not found. Must be 'tikh', 'tv' or 'huber'")

    KtKreg = KtK + alpha**2*regterm

      
    return KtKreg, KtV
# ==============================================================================================



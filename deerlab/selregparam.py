# selregparam.py - Regularization parameter selection
# -----------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np 
import scipy.optimize as opt
import math as m
import deerlab as dl

def selregparam(V, K, r, regtype='tikhonov', method='aic', algorithm='brent',
                nonnegativity=True, noiselvl=-1, regorder=2, weights=1, full_output=False,
                huberparam=1.35, candidates=None):
    r"""
    Selection of optimal regularization parameter based on a selection criterion.

    Parameters 
    ----------
    V : array_like or list of array_like
        Dipolar signal, multiple datasets can be globally evaluated by passing a list of signals.
    K : 2D-array_like or list of 2D-array_like
        Dipolar kernel, if a list of signals is specified, a corresponding list of kernels must be passed as well.
    r : array-like
        Distance axis, in nanometers.
    regtype : string
        Regularization functional type: 
    
        * ``'tikhonov'`` - Tikhonov regularizaton
        * ``'tv'`` - Total variation regularization
        * ``'huber'`` - Huber regularization
        The default is ``'tikhonov'``.

    method : string
        Method for the selection of the optimal regularization parameter.

        * ``'lr'`` - L-curve minimum-radius method (LR)
        * ``'lc'`` - L-curve maximum-curvature method (LC)
        * ``'cv'`` - Cross validation (CV)
        * ``'gcv'`` - Generalized Cross Validation (GCV)
        * ``'rgcv'`` - Robust Generalized Cross Validation (rGCV)
        * ``'srgcv'`` - Strong Robust Generalized Cross Validation (srGCV)
        * ``'aic'`` - Akaike information criterion (AIC)
        * ``'bic'`` - Bayesian information criterion (BIC)
        * ``'aicc'`` - Corrected Akaike information criterion (AICC)
        * ``'rm'`` - Residual method (RM)
        * ``'ee'`` - Extrapolated Error (EE)          
        * ``'ncp'`` - Normalized Cumulative Periodogram (NCP)
        * ``'gml'`` - Generalized Maximum Likelihood (GML)
        * ``'mcl'`` - Mallows' C_L (MCL)
    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting, 
        the default is all weighted equally.
    candidates : list, optional
        List or array of candidate regularization parameter values to be evaluated. 
        If not specified, these are automatically computed from the GSVD of the 
        dipolar kernel and regularization operator.
    regorder : int scalar, optional
        Order of the regularization operator, the default is 2.
    algorithm : string, optional
        Search algorithm: 
        
        * ``'grid'`` - Grid-search, slow.
        * ``'brent'`` - Brent-algorithm, fast.
        The default is ``'brent'``.
        
    full_output : boolean, optional
        If enabled the function will return additional output arguments in a tuple, the default is False.
    nonnegativity : boolean, optional
        Enforces the non-negativity constraint on computed distance distributions, by default enabled.
    noiselvl : float scalar, optional
        Estimate of the noise standard deviation, if not specified it is estimated form the fit residuals.
        Used for the MCL selection method.  
    huberparam : float scalar, optional
        Value of the Huber parameter used in Huber regularization, the default is 1.35.

    Returns
    -------
    alphaopt : scalar
        Optimal regularization parameter.
    alphas : ndarray
        Regularization parameter values candidates evaluated during the search. 
        Returned if full_output is True.
    functional : ndarray
        Values of the selection functional specified by (method) evaluated during the search. 
        Returned if full_output is True.
    residuals : ndarray
        Values of the residual norms evaluated during the search. 
        Returned if full_output is True.
    residuals : ndarray
        Values of the penalty norms evaluated during the search. 
        Returned if full_output is True.
    """
#=========================================================           
    # Ensure the use of numpy-arrays
    r = np.atleast_1d(r)

    # If multiple datasets are passed, concatenate the signals and kernels
    V, K, weights,_ = dl.utils.parse_multidatasets(V, K, weights)

    # The L-curve criteria require a grid-evaluation
    if method == 'lr' or method == 'lc':
        algorithm = 'grid'

    # Get regularization operator
    L = dl.regoperator(r,regorder)
    # Get range of potential alpha values candidates
    if candidates is None:
        alphaCandidates = dl.regparamrange(K,L)
    else: 
        alphaCandidates = np.atleast_1d(candidates)

    # Create function handle
    evalalpha = lambda alpha: _evalalpha(alpha, V, K, L, method, nonnegativity, noiselvl, regtype, weights, huberparam)

    # Evaluate functional over search range, using specified search method
    if algorithm == 'brent':
                    
        # Search boundaries
        lga_min = m.log10(min(alphaCandidates))
        lga_max = m.log10(max(alphaCandidates))

        # Create containers for non-local variables
        functional,residuals,penalties,alphas_evaled = (np.array(0) for _ in range(4))
        def register_ouputs(optout):
        #========================
            nonlocal functional,residuals,penalties,alphas_evaled
            # Append the evaluated outpus at a iteration
            functional = np.append(functional,optout[0])
            residuals = np.append(residuals,optout[1])
            penalties = np.append(penalties,optout[2])
            alphas_evaled = np.append(alphas_evaled,optout[3])
            # Return the last element, current evaluation
            return functional[-1]
        #========================

        # Optimize alpha via Brent's method (implemented in scipy.optimize.fminbound)
        lga_opt = opt.fminbound(lambda lga: register_ouputs(evalalpha(10**lga)), lga_min, lga_max, xtol=0.01)
        alphaOpt = 10**lga_opt

    elif algorithm=='grid':
        
        # Evaluate the full grid of alpha-candidates 
        functional,residuals,penalties,alphas_evaled = tuple(zip(*[evalalpha(alpha) for alpha in alphaCandidates]))

        # If an L-curve method is requested evaluate it now with the full grid:
        
        # L-curve minimum-radius method (LR)
        if method == 'lr':
            Eta = np.log(penalties)
            Rho = np.log(residuals)
            dd = lambda x: (x-np.min(x))/(np.max(x)-np.min(x))
            functional = dd(Rho)**2 + dd(Eta)**2         

        # L-curve maximum-curvature method (LC)
        elif method == 'lc': 
            d1Residual = np.gradient(np.log(residuals))
            d2Residual = np.gradient(d1Residual)
            d1Penalty = np.gradient(np.log(penalties))
            d2Penalty = np.gradient(d1Penalty)
            functional = (d1Residual*d2Penalty - d2Residual*d1Penalty)/(d1Residual**2 + d1Penalty**2)**(3/2)

        # Find minimum of the selection functional              
        alphaOpt = alphaCandidates[np.argmin(functional)]

    else: 
        raise KeyError("Search method not found. Must be either 'brent' or 'grid'.")

    if full_output:
        return alphaOpt,alphas_evaled,functional,residuals,penalties
    else:
        return alphaOpt
#=========================================================


#=========================================================
def _evalalpha(alpha,V,K,L,selmethod,nonneg,noiselvl,regtype,weights,HuberParameter):
    "Evaluation of the selection functional at a given regularization parameter value"

    # Prepare LSQ components
    KtKreg, KtV = dl.lsqcomponents(V,K,L,alpha,weights, regtype=regtype, huberparam=HuberParameter)
    # Solve linear LSQ problem
    if nonneg:
        # Non-negative solution
        P = dl.cvxnnls(KtKreg,KtV)
    else:
        # Unconstrained solution
        P = np.linalg.solve(KtKreg,KtV)

    # Moore-PeNose pseudoinverse
    pK = np.linalg.inv(KtKreg)@K.T
    # Influence matrix
    H = K@pK

    # Residual term
    Residual = np.linalg.norm(K@P - V)
    # Regularization penalty term
    if regtype.lower() == 'tikhonov':
        Penalty = np.linalg.norm(L@P)
    elif regtype.lower() == 'tv':
        Penalty = np.sum(np.sqrt((L@P)**2 + np.finfo(float).eps))
    elif regtype.lower() == 'huber':
        Penalty = np.sum(np.sqrt((L@P/HuberParameter)**2 + 1 ) - 1)
    else:
        raise KeyError("Regularization type not valid. Must be 'tikhonov','tv' or 'huber'")
    #-----------------------------------------------------------------------
    #  Selection methods for optimal regularization parameter
    #-----------------------------------------------------------------------
    
    functional = 0
    N = len(V)

    # Cross validation (CV)
    if  selmethod =='cv': 
        f_ = sum(abs((V - K@P)/(np.ones(N) - np.diag(H)))**2)
            
    # Generalized Cross Validation (GCV)
    elif  selmethod =='gcv': 
        f_ = Residual**2/((1 - np.trace(H)/N)**2)
            
    # Robust Generalized Cross Validation (rGCV)
    elif  selmethod =='rgcv': 
        tuning = 0.9
        f_ = Residual**2/((1 - np.trace(H)/N)**2)*(tuning + (1 - tuning)*np.trace(H**2)/N)
            
    # Strong Robust Generalized Cross Validation (srGCV)
    elif  selmethod =='srgcv':
        tuning = 0.8
        f_ = Residual**2/((1 - np.trace(H)/N)**2)*(tuning + (1 - tuning)*np.trace(pK.T@pK)/N)
            
    # Akaike information criterion (AIC)
    elif  selmethod =='aic': 
        crit = 2
        f_ = N*np.log(Residual**2/N) + crit*np.trace(H)
            
    # Bayesian information criterion (BIC)
    elif  selmethod =='bic':  
        crit = np.log(N)
        f_ = N*np.log(Residual**2/N) + crit*np.trace(H)
            
    # Corrected Akaike information criterion (AICC)
    elif  selmethod =='aicc': 
        crit = 2*N/(N-np.trace(H)-1)
        f_ = N*np.log(Residual**2/N) + crit*np.trace(H)
    
    # Residual method (RM)
    elif  selmethod =='rm':
        scale = K.T@(np.eye(np.shape(H)[0],np.shape(H)[1]) - H)
        f_ = Residual**2/np.sqrt(np.trace(scale.T@scale))

    # Extrapolated Error (EE)          
    elif  selmethod =='ee': 
        f_ = Residual**2/np.linalg.norm(K.T@(K@P - V))

    # Normalized Cumulative Periodogram (NCP)
    elif selmethod == 'ncp': 
        resPeriodogram = abs(np.fft.fft(K@P - V))**2
        wnoisePeriodogram = np.zeros(len(resPeriodogram))
        respowSpectrum = np.zeros(len(resPeriodogram))
        for j in range(len(resPeriodogram)-1):
            respowSpectrum[j]  = np.linalg.norm(resPeriodogram[1:j+1],1)/np.linalg.norm(resPeriodogram[1:-1],1)
            wnoisePeriodogram[j] = j/(len(resPeriodogram) - 1)
        f_ = np.linalg.norm(respowSpectrum - wnoisePeriodogram)

    # Generalized Maximum Likelihood (GML)
    elif  selmethod == 'gml': 
        Treshold = 1e-9
        eigs,_ = np.linalg.eig(np.eye(np.shape(H)[0],np.shape(H)[1]) - H)
        eigs[eigs < Treshold] = 0
        nzeigs = np.real(eigs[eigs!=0])
        f_ = V.T@(V - K@P)/np.prod(nzeigs)**(1/len(nzeigs))

    # Mallows' C_L (MCL)
    elif  selmethod == 'mcl':  
        if noiselvl==-1:
            noiselvl = np.std(V - K@P)
        f_ = Residual**2 + 2*noiselvl**2*np.trace(H) - 2*N*noiselvl**2
        
    elif selmethod == 'lr' or selmethod == 'lc':
        f_ = 0
    else:
        raise ValueError('Selection method \'{}\' is not known.'.format(selmethod))

    functional = functional + f_

    return functional, Residual, Penalty, alpha  
#=========================================================           

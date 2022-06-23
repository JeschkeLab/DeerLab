# selregparam.py - Regularization parameter selection
# -----------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np 
import scipy.optimize as opt
import math as m
import deerlab as dl

def selregparam(y, A, solver, method='aic', algorithm='brent', noiselvl=None,
                searchrange=[1e-8,1e2],regop=None, weights=None, full_output=False, candidates=None):
    r"""
    Selection of optimal regularization parameter based on a selection criterion.

    Parameters 
    ----------
    y : array_like or list of array_like
        Dipolar signal, multiple datasets can be globally evaluated by passing a list of signals.

    A : 2D-array_like or list of 2D-array_like
        Dipolar kernel, if a list of signals is specified, a corresponding list of kernels must be passed as well.

    solver : callable
        Linear least-squares solver. Must be a callable function with signature ``solver(AtA,Aty)``.

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
        Array of weighting coefficients for the individual signals in global fitting.
        If not specified all datasets are weighted inversely proportional to their noise levels.

    algorithm : string, optional
        Search algorithm: 
        
        * ``'grid'`` - Grid-search, slow.
        * ``'brent'`` - Brent-algorithm, fast.
        
        The default is ``'brent'``.

    searchrange : two-element list, optional 
        Search range for the optimization of the regularization parameter with the ``'brent'`` algorithm.
        If not specified the default search range defaults to ``[1e-8,1e2]``.

    candidates : list, optional
        List or array of candidate regularization parameter values to be evaluated with the ``'grid'`` algorithm. 
        If not specified, these are automatically computed from the GSVD of the 
        dipolar kernel and regularization operator. 

    regop : 2D array_like, optional
        Regularization operator matrix, the default is the second-order differential operator.

        
    full_output : boolean, optional
        If enabled the function will return additional output arguments in a tuple, the default is False.

    nonnegativity : boolean, optional
        Enforces the non-negativity constraint on computed distance distributions, by default enabled.

    noiselvl : float scalar, optional
        Estimate of the noise standard deviation, if not specified it is estimated automatically.
        Used for the MCL selection method.  

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

    # If multiple datasets are passed, concatenate the signals and kernels
    y, A, weights,_,__, noiselvl = dl.utils.parse_multidatasets(y, A, weights, noiselvl)

    # The L-curve criteria require a grid-evaluation
    if method == 'lr' or method == 'lc':
        algorithm = 'grid'

    if regop is None:
        L = dl.regoperator(np.arange(np.shape(A)[1]),2)
    else: 
        L = regop

    # Create function handle
    evalalpha = lambda alpha: _evalalpha(alpha, y, A, L, solver, method, noiselvl, weights)

    # Evaluate functional over search range, using specified search method
    if algorithm == 'brent':
                    
        # Search boundaries
        lga_min = m.log10(searchrange[0])
        lga_max = m.log10(searchrange[1])

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
        
        # Get range of potential alpha values candidates
        if candidates is None:
            alphaCandidates = 10**np.linspace(np.log10(searchrange[0]),np.log10(searchrange[1]),60)
        else: 
            alphaCandidates = np.atleast_1d(candidates)

        # Evaluate the full grid of alpha-candidates 
        functional,residuals,penalties,alphas_evaled = tuple(zip(*[evalalpha(alpha) for alpha in alphaCandidates]))

        # If an L-curve method is requested evaluate it now with the full grid:
        
        # L-curve minimum-radius method (LR)
        if method == 'lr':
            Eta = np.log(np.asarray(penalties)+1e-20)
            Rho = np.log(np.asarray(residuals)+1e-20)
            dd = lambda x: (x-np.min(x))/(np.max(x)-np.min(x))
            functional = dd(Rho)**2 + dd(Eta)**2         
            
        # L-curve maximum-curvature method (LC)
        elif method == 'lc': 
            d1Residual = np.gradient(np.log(np.asarray(residuals)+1e-20))
            d2Residual = np.gradient(d1Residual)
            d1Penalty = np.gradient(np.log(np.asarray(penalties)+1e-20))
            d2Penalty = np.gradient(d1Penalty)
            functional = (d1Residual*d2Penalty - d2Residual*d1Penalty)/(d1Residual**2 + d1Penalty**2)**(3/2)
            functional = -functional # Maximize instead of minimize 

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
def _evalalpha(alpha,y,A,L,solver,selmethod,noiselvl,weights):
    "Evaluation of the selection functional at a given regularization parameter value"

    # Prepare LSQ components
    AtAreg, Aty = dl.solvers._lsqcomponents(y,A,L,alpha,weights)
    wA = weights[:,np.newaxis]*A
    # Solve linear LSQ problem
    P = solver(AtAreg,Aty)

    # Moore-PeNose pseudoinverse
    pA = np.linalg.inv(AtAreg)@wA.T
    # Influence matrix
    H = wA@pA

    # Residual term
    residuals = weights*(A@P - y)
    Residual = np.linalg.norm(residuals)
    # Regularization penalty term
    Penalty = np.linalg.norm(L@P)
    #-----------------------------------------------------------------------
    #  Selection methods for optimal regularization parameter
    #-----------------------------------------------------------------------
    
    functional = 0
    N = len(y)

    # Cross validation (CV)
    if  selmethod =='cv': 
        f_ = sum(abs(residuals/(np.ones(N) - np.diag(H)))**2)
            
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
        f_ = Residual**2/((1 - np.trace(H)/N)**2)*(tuning + (1 - tuning)*np.trace(pA.T@pA)/N)
            
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
        scale = A.T@(np.eye(np.shape(H)[0],np.shape(H)[1]) - H)
        f_ = Residual**2/np.sqrt(np.trace(scale.T@scale))

    # Extrapolated Error (EE)          
    elif  selmethod =='ee': 
        f_ = Residual**2/np.linalg.norm(A.T@(residuals))

    # Normalized Cumulative Periodogram (NCP)
    elif selmethod == 'ncp': 
        resPeriodogram = abs(np.fft.fft(residuals))**2
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
        f_ = y.T@(-residuals)/np.prod(nzeigs)**(1/len(nzeigs))

    # Mallows' C_L (MCL)
    elif  selmethod == 'mcl':  
        f_ = Residual**2 + 2*noiselvl**2*np.trace(H) - 2*N*noiselvl**2
        
    elif selmethod == 'lr' or selmethod == 'lc':
        f_ = 0
    else:
        raise ValueError(f'Selection method \'{selmethod}\' is not known.')

    functional = functional + f_

    return functional, Residual, Penalty, alpha  
#=========================================================           

# fitregmodel.py - Fiting of regularization problems
# --------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import matplotlib.pyplot as plt
from deerlab.nnls import fnnls, cvxnnls, nnlsbpp
import deerlab as dl
import copy
from deerlab.utils import hccm, goodness_of_fit
from deerlab.classes import UncertQuant, FitResult

def fitregmodel(V,K,r, regtype='tikhonov', regparam='aic', regorder=2, solver='cvx', 
                weights=1, huberparam=1.35, nonnegativity=True, obir = False, 
                uqanalysis=True, renormalize=True, noiselevelaim = -1):
    r"""
    Fits a non-parametric distance distribution to one (or several) signals using regularization aproaches.

    Parameters 
    ----------
    V : array_like or list of array_like  
        Dipolar signal, multiple datasets can be globally evaluated by passing a list of signals.
    
    K : 2D-array_like array or list of 2D-array_like
        Dipolar kernel, if a list of signals is specified, a corresponding list of kernels must be passed as well.
    
    r : array_like
        Distance axis, in nanometers.
    
    regtype : string, optional
        Regularization functional type: 
    
        * ``'tikhonov'`` - Tikhonov regularizaton
        * ``'tv'``  - Total variation regularization
        * ``'huber'`` - Huber regularization
        The default is ``'tikhonov'``.   
    
    regparam : string or scalar, optional
        Method for the automatic selection of the optimal regularization parameter:
    
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
        The regularization parameter can be manually specified by passing a scalar value instead of a string.
        The default ``'aic'``.
    
    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting, the default is all weighted equally.
    
    regorder : int scalar, optional
        Order of the regularization operator, the default is 2.
    
    solver : string, optional
        Optimizer used to solve the non-negative least-squares problem: 

        * ``'cvx'`` - Optimization of the NNLS problem using cvxopt
        * ``'fnnls'`` - Optimization using the fast NNLS algorithm.
        * ``'nnlsbpp'`` - Optimization using the block principal pivoting NNLS algorithm.
        The default is ``'cvx'``.

    uqanalysis : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled. 
    
    nonnegativity : boolean, optional
        Enforces the non-negativity constraint on computed distance distributions, by default it is enabled.
    
    huberparam : scalar, optional
        Value of the Huber parameter used in Huber regularization, the default is 1.35.
    
    renormalize : boolean, optional
        Enable/disable renormalization of the fitted distribution, by default it is enabled.
    
    obir : boolean, optional
        Enable/disable the use of the Osher-Bregman iterated regularization algorithm, by default it is disabled.
    
    noiselevelaim : scalar, optional
        Noise level at which to stop the OBIR algorithm. If not specified it is automatically estimated from the fit residuals.


    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    P : ndarray
        Fitted distance distribution
    
    uncertainty : :ref:`UncertQuant`
        Covariance-based uncertainty quantification of the fitted distance distribution
    
    alpha : float int
        Regularization parameter used in the optimization
    
    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    
    plot : callable
        Function to display the results. It will 
        display the fitted signals and distance distributions with
        confidence intervals. If requested, the function returns 
        the `matplotlib.axes` object as output. 
    
    stats : dict
        Goodness of fit statistical estimators (if ``full_output=True``):

        * ``stats['chi2red']`` - Reduced \chi^2 test
        * ``stats['r2']`` - R^2 test
        * ``stats['rmsd']`` - Root-mean squared deviation (RMSD)
        * ``stats['aic']`` - Akaike information criterion
        * ``stats['aicc']`` - Corrected Akaike information criterion
        * ``stats['bic']`` - Bayesian information criterion
    
    success : bool
        Whether or not the optimizer exited successfully.
    
    cost : float
        Value of the cost function at the solution.
    
    residuals : ndarray
        Vector of residuals at the solution.

    """
    # Prepare signals, kernels and weights if multiple are passed
    V, K, weights, subsets, prescales = dl.utils.parse_multidatasets(V, K, weights,precondition=True)

    # Compute regularization matrix
    L = dl.regoperator(r,regorder)

    # Determine the type of problem to solve
    if obir:
        problem = 'obir'
    elif nonnegativity:
        problem = 'nnls'
    else:
        problem = 'unconstrained'

    # If the regularization parameter is not specified, get the optimal choice
    if type(regparam) is str:
        alpha = dl.selregparam(V,K,r,regtype,regparam,regorder=regorder,weights=weights, nonnegativity=nonnegativity,huberparam=huberparam)
    else:
        alpha = regparam

    # Prepare components of the LSQ-problem
    [KtKreg, KtV] = dl.lsqcomponents(V,K,L,alpha,weights, regtype=regtype, huberparam=huberparam)

    # Unconstrained LSQ problem
    if problem == 'unconstrained':
        Pfit = np.linalg.solve(KtKreg,KtV)

    # Osher-Bregman iterated regularization
    elif problem == 'obir':
        Pfit = _obir(V,K,L,regtype,alpha,weights,noiselevelaim=noiselevelaim,huberparam = huberparam, solver = solver)
    # Non-negative LSQ problem
    elif problem == 'nnls':

        if solver == 'fnnls':
            Pfit = fnnls(KtKreg,KtV)
        elif solver == 'nnlsbpp':
            Pfit = nnlsbpp(KtKreg,KtV,np.linalg.solve(KtKreg,KtV))
        elif solver == 'cvx':
            Pfit = cvxnnls(KtKreg, KtV)
        else:
            raise KeyError(f'{solver} is not a known non-negative least squares solver')

    # Get fit final status
    Vfit = K@Pfit
    success = ~np.all(Pfit==0)
    res = V - Vfit
    fval = np.linalg.norm(V - Vfit)**2 + alpha**2*np.linalg.norm(L@Pfit)**2
        
    # Uncertainty quantification
    # ----------------------------------------------------------------
    if uqanalysis:
        # Construct residual parts for for the residual and regularization terms
        res = weights*(V - K@Pfit)

        # Construct Jacobians for the residual and penalty terms
        Jres = K*weights[:,np.newaxis]
        res,J = _augment(res,Jres,regtype,alpha,L,Pfit,huberparam)

        # Calculate the heteroscedasticity consistent covariance matrix 
        covmat = hccm(J,res,'HC1')
        
        # Construct confidence interval structure for P
        NonNegConst = np.zeros(len(r))
        Puq = UncertQuant('covariance',Pfit,covmat,NonNegConst,[])
    else:
        Puq = UncertQuant('void')

    
    # Re-normalization of the distributions
    # --------------------------------------
    postscale = np.trapz(Pfit,r)
    if renormalize:
        Pfit = Pfit/postscale
        if uqanalysis:
            Puq_ = copy.deepcopy(Puq) # need a copy to avoid infite recursion on next step
            Puq.ci = lambda p: Puq_.ci(p)/postscale


    # Goodness-of-fit
    # --------------------------------------
    H = K@(np.linalg.pinv(KtKreg)@K.T)
    stats = []
    for subset in subsets: 
        Ndof = len(V[subset]) - np.trace(H)
        stats.append(goodness_of_fit(V[subset],Vfit[subset],Ndof))
    if len(stats)==1: 
        stats = stats[0]


    # Signal amplitude scales
    # ---------------------------------------
    scales = []
    for i in range(len(subsets)): 
        scales.append(prescales[i]*postscale)
    if len(scales)==1:
        scales = scales[0]

    # Results display
    # ---------------------------------------
    plotfcn = lambda: _plot(subsets,V,Vfit,r,Pfit,Puq)

    return FitResult(P=Pfit, V=Vfit, uncertainty=Puq, regparam=alpha, scale=scales, stats=stats, 
                     plot=plotfcn, cost=fval, residuals=res, success=success)
# ===========================================================================================


def _augment(res,J,regtype,alpha,L,P,eta):
# ===========================================================================================
    """ 
    LSQ residual and Jacobian augmentation
    =======================================

    Augments the residual and the Jacobian of a LSQ problem to include the
    regularization penalty. The residual and Jacobian contributions of the 
    specific regularization methods are analytically introduced. 
    """
    eps = np.finfo(float).eps
    # Compute the regularization penalty augmentation for the residual and the Jacobian
    if regtype == 'tikhonov':
        resreg = L@P
        Jreg = L
    elif regtype == 'tv':
        resreg =((L@P)**2 + eps)**(1/4)
        Jreg = 2/4*((( ( (L@P)**2 + eps)**(-3/4) )*(L@P))[:, np.newaxis]*L)
    elif regtype == 'huber':
        resreg = np.sqrt(np.sqrt((L@P/eta)**2 + 1) - 1)
        Jreg = 0.5/(eta**2)*( (((np.sqrt((L@P/eta)**2 + 1) - 1 + eps)**(-1/2)*(((L@P/eta)**2 + 1+ eps)**(-1/2)))*(L@P))[:, np.newaxis]*L )

    # Include regularization parameter
    resreg = alpha*resreg
    Jreg = alpha*Jreg

    # Augment jacobian and residual
    res = np.concatenate((res,resreg))
    J = np.concatenate((J,Jreg))

    return res,J
# ===========================================================================================


def _obir(V,K,L, regtype, alpha, weights, noiselevelaim=-1, huberparam=1.35 , solver = 'cvx'):
# ===========================================================================================
    """
    Osher's Bregman-iterated regularization method
    ==============================================

    P = OBIR(V,K,r,'type',alpha)

    OBIR of the N-point signal (V) to a M-point distance
    distribution (P) given a M-point distance axis (r) and NxM point kernel
    (K). The regularization parameter (alpha) controls the regularization
    properties.

    The type of regularization employed in OBIR is set by the 'type'
    input argument. The regularization models implemented in OBIR are:
        'tikhonov' -   Tikhonov regularization
        'tv'       -   Total variation regularization
        'huber'    -   pseudo-Huber regularization

    P = OBIR(...,'Property',Value)
    Additional (optional) arguments can be passed as name-value pairs.
    """

    if noiselevelaim == -1:
        noiselevelaim = dl.noiselevel(V)

    MaxOuterIter = 5000
    stopDivergent = False

    # Preparation
    #-------------------------------------------------------------------------------

    # Initialize
    nr = np.shape(L)[1]
    subGrad = np.zeros(nr)
    Counter = 1
    Iteration = 1
    Pfit = np.zeros(nr)

    # Osher's Bregman Iterated Algorithm
    #-------------------------------------------------------------------------------
    diverged = np.zeros(len(V))
    semiconverged = np.zeros(len(V))

    # Precompute the KtK and KtS input arguments
    [KtKreg,KtV] = dl.lsqcomponents(V,K,L,alpha, weights, regtype=regtype, huberparam=huberparam)

    while Iteration <= MaxOuterIter:

        # Store previous iteration distribution
        Pprev = Pfit
        
        #Update
        KtVsg = KtV - subGrad
        
        #Get solution of current Bregman iteration
        if solver == 'fnnls':
            Pfit = fnnls(KtKreg,KtVsg)
        elif solver == 'nnlsbpp':
            Pfit = nnlsbpp(KtKreg,KtVsg,np.linalg.solve(KtKreg,KtVsg))
        elif solver == 'cvx':
            Pfit = cvxnnls(KtKreg, KtVsg)
        else:
            raise KeyError(f'{solver} is not a known non-negative least squares solver')
        
        # Update subgradient at current solution
        subGrad = subGrad + weights*K.T@(K@Pfit - V)
        
        # Iteration control
        #--------------------------------------------------------------------------
        if Iteration == 1:
            # If at first iteration, the residual deviation is already below the
            # noise deviation then impose oversmoothing and remain at first iteration
            if noiselevelaim  > np.std(K@Pfit - V):
                alpha = alpha*2**Counter
                Counter = Counter + 1
                
                # Recompute the KtK and KtS input arguments with new alpha
                KtKreg,KtV = dl.lsqcomponents(V,K,L,alpha, weights, regtype=regtype, huberparam=huberparam)
            else:
                # Once the residual deviation is above the threshold, then proceed
                # further with the Bregman iterations
                Iteration  = Iteration + 1
        else:
            # For the rest of the Bregman iterations control the condition and stop
            # when fulfilled
            diverged = np.std(K@Pprev - V) < np.std(K@Pfit - V)
            semiconverged = noiselevelaim > np.std(K@Pfit - V)

            if semiconverged:
                break
            else:
                Iteration  = Iteration + 1

            # If residual deviation starts to diverge, stop
            if stopDivergent and diverged:
                Pfit = Pprev
                break

    return Pfit
# ===========================================================================================

def _plot(subsets,Vexp,Vfit,r,Pfit,Puq):
# ===========================================================================================
    nSignals = len(subsets)
    _,axs = plt.subplots(nSignals+1,figsize=[7,3+3*nSignals])
    for i in range(nSignals): 
        subset = subsets[i]
        # Plot the experimental signal and fit
        axs[i].plot(Vexp[subset],'.',color='grey',alpha=0.5)
        axs[i].plot(Vfit[subset],'tab:blue')
        axs[i].grid(alpha=0.3)
        axs[i].set_xlabel('Array Elements')
        axs[i].set_ylabel('V[{}]'.format(i))
        axs[i].legend(('Data','Fit'))

    # Get confidence intervals for the distance distribution
    Pci95 = Puq.ci(95)
    Pci50 = Puq.ci(50)
    # Plot the distribution
    axs[nSignals].plot(r,Pfit,'tab:blue')
    axs[nSignals].fill_between(r,Pci95[:,0], Pci95[:,1],facecolor='tab:blue',linestyle='None',alpha=0.2)
    axs[nSignals].fill_between(r,Pci50[:,0], Pci50[:,1],facecolor='tab:blue',linestyle='None',alpha=0.4)
    axs[nSignals].set_xlabel('Distance [nm]')
    axs[nSignals].set_ylabel('P [nm⁻¹]')
    axs[nSignals].legend(('Fit','95%-CI','50%-CI'))
    axs[nSignals].grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
    return axs
# ===========================================================================================




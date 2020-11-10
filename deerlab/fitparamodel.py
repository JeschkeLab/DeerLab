# fitparamodel.py - Parametric dipolar signal model fit function
# ---------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.utils import isempty, multistarts, hccm, parse_multidatasets, goodness_of_fit, Jacobian
from deerlab.classes import UncertQuant, FitResult
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import warnings

def fitparamodel(V, model, par0=[],lb=[],ub=[], weights = 1,
                 multistart=1, tol=1e-10, maxeval=5000, maxiter = 3000,
                 rescale=True, uqanalysis=True, covmatrix=[]):
    r""" Fits the dipolar signal(s) to a parametric model using non-linear least-squares.

    Parameters
    ----------
    V : array_like or list of array_like
        Dipolar signal(s) to be fitted.
    
    Vmodel : callable
        Function taking an array of parameters and returning the dipolar signal(s) as an array or a list of arrays.
    
    par0  : array_like
        Start values of the model parameters.
    
    lb : array_like      
        Lower bounds for the model parameters. If not specified, it is left unbounded.
    
    ub : array_like     
        Upper bounds for the model parameters. If not specified, it is left unbounded.
    
    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting, the default is all weighted equally.
    
    rescale : boolean, optional
        Enable/disable optimization of the signal scale, by default it is enabled.
    
    multiStarts : scalar, optional
        Number of starting points for global optimization, the default is 1.
    
    uqanalysis : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.    
    
    tol : scalar, optional
        Optimizer function tolerance, the default is 1e-10.
    
    maxeval : scalar, optional
        Maximum number of optimizer iterations, the default is 5000.
    
    maxiter : scalar, optional
        Maximum number of optimizer iterations, the default is 3000.
    
    covmatrix : array_like with shape(n,n), optional
        Covariance matrix of the noise in the dataset(s). If not specified it is automatically computed.


    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    param : ndarray
        Fitted model parameters
    
    uncertainty : :ref:`UncertQuant`
        Covariance-based uncertainty quantification of the fitted parameters.
    
    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    
    plot : callable
        Function to display the results. It will 
        display the fitted signals. If requested, the function returns 
        the ``matplotlib.axes`` object as output. 
    
    stats :  dict
        Goodness of fit statistical estimators:

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


    Examples
    --------
    A simple example for fitting a fully parametric 4-pulse DEER signal::

        K0 = dl.dipolarkernel(t,r)
        def Vmodel(par):
            # Extract parameters
            rmean,fwhm,lam,conc = par
            B = dl.bg_hom3d(t,conc)
            P = dl.dd_gauss(r,[rmean,fwhm])
            V = (1-lam + lam*K0@P)*B
            return V
        
        # par: rmean fwhm lam conc
        par0 = [3.5, 0.5, 0.2, 250] # start values
        lb   = [ 1,  0.1,  0,  50 ] # lower bounds
        ub   = [20,   5,   1,  800] # upper bounds
        fit = dl.fitparamodel(V,Vmodel,par0,lb,ub)


    If multiple datasets are to be fitted globally, this can be easily achieved by adapting the example above. Differentiating between local and global parameters is also simple. This example shows the definition of a full parametric model of two signals with same underlying distance distribution and background but different modulation depths::
    
        K1 = dl.dipolarkernel(t1,r)
        K1 = dl.dipolarkernel(t2,r)
        def Vmodel(par):
            # Extract parameters
            rmean,fwhm,lam1,lam,2conc = par
            B = dl.bg_hom3d(t,conc)
            P = dl.dd_gauss(r,[rmean,fwhm])
            V1 = (1-lam1 + lam1*K1@P)*B
            V2 = (1-lam2 + lam2*K2@P)*B
            return V
        
        # par: rmean fwhm lam1 lam2 conc
        par0 = [3.5, 0.5, 0.2, 0.3, 250] # start values
        lb   = [ 1,  0.1,  0,   0,  50 ] # lower bounds
        ub   = [20,   5,   1,   1,  800] # upper bounds
        fit = dl.fitparamodel([V1,V2],Vmodel,par0,lb,ub)

    """

    if isempty(par0):
        raise KeyError('The start values par0 of the parameters must be specified')

    lb,ub = np.atleast_1d(lb,ub)
    V, model, weights, Vsubsets = parse_multidatasets(V, model, weights)
    Nsignals = len(Vsubsets)
    scales = [1]*Nsignals

    if not np.all(np.isreal(V)):
        raise ValueError('The input signal(s) cannot be complex-valued.')

    if isempty(lb):
        lb = np.full_like(par0, -np.inf)
    if isempty(ub):
        ub = np.full_like(par0, +np.inf)

    # Prepare upper/lower bounds on parameter search range
    unboundedparams = np.any(np.isinf(lb)) or np.any(np.isinf(ub))

    if  len(par0)!=len(ub) or len(par0)!=len(lb):
        raise ValueError('The inital guess and upper/lower boundaries must have equal length.')
    if np.any(ub<lb):
        raise ValueError('Lower bound values cannot be larger than upper bound values.')

    # Preprare multiple start global optimization if requested
    if multistart>1 and unboundedparams:
        raise ValueError('multistart optimization cannot be used with unconstrained parameters.')
    multistarts_par0 = multistarts(multistart,par0,lb,ub)

    def lsqresiduals(p):
    # -------------------------------------------------------------------------------    
        """
        Residual vector function
        ------------------------------------------------------------------
        Function that provides vector of residuals, which is the objective
        function for the least-squares solvers
        """
        nonlocal scales

        Vsim = model(p)

        # Check if there are invalid values...
        if any(np.isnan(Vsim)) or any(np.isinf(Vsim)):
            res = np.zeros_like(Vsim) # ...can happen when Jacobian is evaluated outside of bounds
            return res

        # Otherwise if requested, compute the scale of the signal via linear LSQ   
        if rescale:
            scales = []
            for subset in Vsubsets:
                Vsim_,V_ = (V[subset] for V in [Vsim, V]) # Rescale the subsets corresponding to each signal
                scale = np.squeeze(np.linalg.lstsq(np.atleast_2d(Vsim_).T,np.atleast_2d(V_).T,rcond=None)[0])
                Vsim[subset] = scale*Vsim_  
                scales.append(scale) # Store the optimized scales of each signal

        # Compute the residual
        res = weights*(V-Vsim)
        
        return res
    # -------------------------------------------------------------------------------    
    fvals = []
    parfits = []
    sols = []
    for par0 in multistarts_par0:
        # Solve the non-linear least squares (NLLS) problem 
        sol = least_squares(lsqresiduals ,par0, bounds=(lb,ub), max_nfev=int(maxiter), ftol=tol, method='dogbox')
        sols.append(sol)
        parfits.append(sol.x)
        fvals.append(sol.cost)        

    # Find global minimum from multiple runs
    globmin = np.argmin(fvals)
    parfit = parfits[globmin]
    sol = sols[globmin]
    # Issue warnings if fitted parameter values are at search range boundaries
    tol = 1e-5
    atLower = abs(parfit-lb)<tol
    atUpper = abs(parfit-ub)<tol
    fixedPar = lb==ub
    for p in range(len(parfit)):
        if fixedPar[p]:
            # Skip if parameter is fixed
            continue
        if atLower[p]:
            warnings.warn('The fitted value of parameter #{}, is at the lower bound ({}).'.format(p,lb[p]))
        if atUpper[p]:
            warnings.warn('The fitted value of parameter #{},  is at the upper bound ({}).'.format(p,ub[p]))

    # Calculate parameter confidence intervals
    if uqanalysis:
        
        # Compute residual vector and estimate variance from that
        residuals = lsqresiduals(parfit)
        
        # Calculate numerical estimates of the Jacobian and Hessian
        # of the negative log-likelihood
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            J = Jacobian(lsqresiduals,parfit,lb,ub)

        # Estimate the heteroscedasticity-consistent covariance matrix
        if isempty(covmatrix):
            # Use estimated data covariance matrix
            covmatrix = hccm(J,residuals,'HC1')
        else:
            # Use user-given data covariance matrix
            covmatrix = hccm(J,covmatrix)
        # Construct confidence interval structure
        paruq = UncertQuant('covariance',parfit,covmatrix,lb,ub)
    else:
        paruq = UncertQuant('void')

    # Calculate goodness of fit
    stats = []
    Vfit = model(parfit)
    for subset in Vsubsets: 
        Ndof = len(V[subset]) - len(par0)
        stats.append(goodness_of_fit(V[subset],Vfit[subset],Ndof))
    if Nsignals==1: 
        stats = stats[0]
        scales = scales[0]

    # Get plot function
    plotfcn = lambda: _plot(Vsubsets,V,Vfit)

    return FitResult(
            param=parfit, uncertainty=paruq, scale=scales, stats=stats, cost=fvals,
            plot=plotfcn, residuals=sol.fun, success=sol.success)

def _plot(Vsubsets,V,Vfit):
    nSignals = len(Vsubsets)
    _,axs = plt.subplots(nSignals,figsize=[7,3*nSignals])
    axs = np.atleast_1d(axs)
    for i in range(nSignals): 
        subset = Vsubsets[i]
        # Plot the experimental signal and fit
        axs[i].plot(V[subset],'.',color='grey',alpha=0.5)
        axs[i].plot(Vfit[subset],'tab:blue')
        axs[i].grid(alpha=0.3)
        axs[i].set_xlabel('Array Elements')
        axs[i].set_ylabel('V[{}]'.format(i))
        axs[i].legend(('Data','Fit'))

    plt.tight_layout()
    plt.show()
    return axs
# fitparamodel.py - Parametric dipolar signal model fit function
# ---------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.utils import isempty, multistarts, jacobianest, hccm, parse_multidatasets, goodness_of_fit
from deerlab import uqst
from scipy.optimize import least_squares
import warnings

def fitparamodel(V, model, par0=[],lb=[],ub=[], weights = 1, MultiStart=1, tolFun=1e-10, maxFunEvals=5000, maxIter = 3000, rescale=True, uqanalysis=True, covmatrix=[]):
    """ Fits the dipolar signal(s) to a parametric model using non-linear least-squares.

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

    Returns
    -------
    param : ndarray
        Fitted model parameters
    paramuq : obj
        Covariance-based uncertainty quantification of the fitted parameters.
    stats :  dict
        Goodness of fit statistical estimators:

        * stats['chi2red'] - Reduced \chi^2 test
        * stats['r2'] - R^2 test
        * stats['rmsd'] - Root-mean squared deviation (RMSD)
        * stats['aic'] - Akaike information criterion
        * stats['aicc'] - Corrected Akaike information criterion
        * stats['bic'] - Bayesian information criterion

    Other parameters
    ----------------
    weights : array_like 
        Array of weighting coefficients for the individual signals in global fitting, the default is all weighted equally.
    rescale : boolean
        Enable/disable optimization of the signal scale, by default it is enabled.
    multiStarts : scalar
        Number of starting points for global optimization, the default is 1.
    uqanalysis : boolean
        Enable/disable the uncertainty quantification analysis, by default it is enabled.    
    tolFun : scalar
        Optimizer function tolerance, the default is 1e-10.
    maxFunEvals : scalar
        Maximum number of optimizer iterations, the default is 5000.
    maxIter : scalar
        Maximum number of optimizer iterations, the default is 3000.
    covmatrix : array_like with shape(n,n)
        Covariance matrix of the noise in the dataset(s). If not specified it is automatically computed.

    Examples
    --------
    A simple example for fitting a fully parametric 4-pulse DEER signal::

        K0 = dipolarkernel(t,r)
        def Vmodel(par):
            # Extract parameters
            rmean,fwhm,lam,conc = par
            B = bg_hom3d(t,conc)
            P = dd_gauss(r,[rmean,fwhm])
            V = (1-lam + lam*K0@P)*B[:,np,newaxis]
            return V
        
        # par: rmean fwhm lam conc
        par0 = [3.5, 0.5, 0.2, 250] # start values
        lb   = [ 1,  0.1,  0,  50 ] # lower bounds
        ub   = [20,   5,   1,  800] # upper bounds
        parfit,paruq,stats = fitparamodel(V,Vmodel,par0,lb,ub)


    If multiple datasets are to be fitted globally, this can be easily achieved by adapting the example above. Differentiating between local and global parameters is also simple. This example shows the definition of a full parametric model of two signals with same underlying distance distribution and background but different modulation depths::
    
        K1 = dipolarkernel(t1,r)
        K1 = dipolarkernel(t2,r)
        def Vmodel(par):
            # Extract parameters
            rmean,fwhm,lam1,lam,2conc = par
            B = bg_hom3d(t,conc)
            P = dd_gauss(r,[rmean,fwhm])
            V1 = (1-lam1 + lam1*K1@P)*B[:,np,newaxis]
            V2 = (1-lam2 + lam2*K2@P)*B[:,np,newaxis]
            return V
        
        # par: rmean fwhm lam1 lam2 conc
        par0 = [3.5, 0.5, 0.2, 0.3, 250] # start values
        lb   = [ 1,  0.1,  0,   0,  50 ] # lower bounds
        ub   = [20,   5,   1,   1,  800] # upper bounds
        parfit,paruq,stats = fitparamodel([V1,V2],Vmodel,par0,lb,ub)

    """

    lb,ub = np.atleast_1d(lb,ub)
    V, model, weights, Vsubsets = parse_multidatasets(V, model, weights)

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
    if MultiStart>1 and unboundedparams:
        raise ValueError('Multistart optimization cannot be used with unconstrained parameters.')
    multistarts_par0 = multistarts(MultiStart,par0,lb,ub)

    def lsqresiduals(p):
    # -------------------------------------------------------------------------------    
        """
        Residual vector function
        ------------------------------------------------------------------
        Function that provides vector of residuals, which is the objective
        function for the least-squares solvers
        """
        Vsim = model(p)

        # Check if there are invalid values...
        if any(np.isnan(Vsim)) or any(np.isinf(Vsim)):
            res = np.zeros_like(Vsim) # ...can happen when jacobianest() evaluates outside of bounds
            return res

        # Otherwise if requested, compute the scale of the signal via linear LSQ   
        if rescale:
            scale = np.squeeze(np.linalg.lstsq(np.atleast_2d(Vsim).T,np.atleast_2d(V).T,rcond=None)[0])
            Vsim = scale*Vsim
        
        # Compute the residual
        res = weights*(V-Vsim)
        
        return res
    # -------------------------------------------------------------------------------    
    fvals = []
    parfits = []
    for par0 in multistarts_par0:
        # Solve the non-linear least squares (NLLS) problem 
        sol = least_squares(lsqresiduals ,par0, bounds=(lb,ub), max_nfev=int(maxIter), ftol=tolFun, method='dogbox')
        parfits.append(sol.x)
        fvals.append(sol.cost)        

    # Find global minimum from multiple runs
    globmin = np.argmin(fvals)
    parfit = parfits[globmin]

    # Issue warnings if fitted parameter values are at search range boundaries
    tol = 1e-5
    atLower = abs(parfit-lb)<tol
    atUpper = abs(parfit-ub)<tol
    fixedPar = lb==ub
    for p in range(len(parfit)):
        # Skip if parameter is fixed
        if fixedPar[p]:
            continue
        if atLower[p]:
            warnings.warn('The fitted value of parameter #{}, is at the lower bound of the range.'.format(p))
        if atUpper[p]:
            warnings.warn('The fitted value of parameter #{},  is at the upper bound of the range.'.format(p))

    # Calculate parameter confidence intervals
    if uqanalysis:
        
        # Compute residual vector and estimate variance from that
        residuals = lsqresiduals(parfit)
        
        # Calculate numerical estimates of the Jacobian and Hessian
        # of the negative log-likelihood
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            J,_ = jacobianest(lsqresiduals,parfit)
        
        # Estimate the heteroscedasticity-consistent covariance matrix
        if isempty(covmatrix):
            # Use estimated data covariance matrix
            covmatrix = hccm(J,residuals,'HC1')
        else:
            # Use user-given data covariance matrix
            covmatrix = hccm(J,covmatrix)
        # Construct confidence interval structure
        paruq = uqst('covariance',parfit,covmatrix,lb,ub)
    else:
        paruq = []


    # Calculate goodness of fit
    stats = []
    for subset in Vsubsets: 
        Ndof = len(V[subset]) - len(par0)
        stats.append(goodness_of_fit(V[subset],model(parfit)[subset],Ndof))
    if len(stats)==1: 
        stats = stats[0]

    return parfit, paruq, stats


# solvers.py - Collection of least-squares solvers
# --------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

# External dependencies
import numpy as np
import cvxopt as cvx
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, lsq_linear
# DeerLab dependencies
import deerlab as dl
from deerlab.classes import UQResult, FitResult
from deerlab.utils import multistarts, hccm, parse_multidatasets, goodness_of_fit, Jacobian, isempty


def _plot(ys,yfits,show):
# ===========================================================================================
    """
    Plot method for the FitResult object
    ====================================

    Plots the input dataset(s), their fits, and uncertainty bands.
    """
    nSignals = len(ys)
    fig,axs = plt.subplots(nSignals,figsize=[7,3*nSignals])
    axs = np.atleast_1d(axs)
    for i,(y,yfit) in enumerate(zip(ys,yfits)): 
        # Plot the experimental signal and fit
        axs[i].plot(y,'.',color='grey',alpha=0.5)
        axs[i].plot(yfit,'tab:blue')
        axs[i].grid(alpha=0.3)
        axs[i].set_xlabel('Array elements')
        axs[i].set_ylabel(f'Data #{i}')
        axs[i].legend(('Data','Fit'))
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.close()
    return fig
# ===========================================================================================

# ===========================================================================================
def _check_bounds(lb,ub,par0=None,N=None):
    """
    Boundary conditions parsing and validation
    ==========================================

    Parses the input boundary conditions and determines whether they are correct.
    """
    if par0 is not None:
        N = len(par0)

    if lb is None or isempty(lb):
        lb = np.full(N, -np.inf)
    if ub is None or isempty(ub):
        ub = np.full(N, np.inf)
    lb, ub = np.atleast_1d(lb, ub)

    # Check that the boundaries are valid
    if np.any(ub < lb):
        raise ValueError('The upper bounds cannot be larger than the lower bounds.')
    if par0 is not None:
        if len(lb) != len(par0) or len(ub) != len(par0):
            raise TypeError('The lower/upper bounds and start values must have the same number of elements')
                # Check that the non-linear start values are inside the box constraint
        if np.any(par0 > ub) or np.any(par0 < lb):
            raise ValueError('The start values are outside of the specified bounds.')
    if N is not None:
        if len(lb) != N or len(ub) != N:
            raise TypeError(f'The lower/upper bounds and start values must have {N} elements')

    # Check if the nonlinear problem is bounded
    isbounded = (not np.all(np.isinf(lb))) or (not np.all(np.isinf(ub)))

    return lb,ub,isbounded
# ===========================================================================================


# ==============================================================================================
def _lsqcomponents(V, K, L=None, alpha=0, weights=None):
    """
    Linear least-squares components
    ===============================

    Calculate the components needed for the linear least-squares (LSQ) solvers. 
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
    
    # Compute the regularization term
    regterm = L.T@L
    
    KtKreg = KtK + alpha**2*regterm
    
    return KtKreg, KtV
# ==============================================================================================


# ===========================================================================================
def _prepare_linear_lsq(A,lb,ub,reg,L,tol,maxiter,nnlsSolver):
    """
    Preparation of linear least-squares 
    ===================================

    Evaluates the conditions of the linear least-squares problem and determines:
    - whether regularization is required
    - the type of boundary constraints
    - optimal linear least-squares solver to use
    """
    Nlin = np.shape(A)[1]

    # Determine whether to use regularization penalty
    if reg == 'auto':
        illConditioned = np.linalg.cond(A) > 10
        includeRegularization  = illConditioned
    else:
        includeRegularization  = reg

    # Check if the nonlinear and linear problems are constrained
    linearConstrained = (not np.all(np.isinf(lb))) or (not np.all(np.isinf(ub)))
    # Check for non-negativity constraints on the linear solution
    nonNegativeOnly = (np.all(lb == 0)) and (np.all(np.isinf(ub)))

    # Use an arbitrary axis
    axis = np.arange(Nlin)
    if L is None and includeRegularization: 
        L = dl.regoperator(axis,2)


    # Prepare the linear solver
    # ----------------------------------------------------------
    if not linearConstrained:
        # Unconstrained linear LSQ
        linSolver = np.linalg.solve
        parseResult = lambda result: result

    elif linearConstrained and not nonNegativeOnly:
        # Constrained linear LSQ
        linSolver = lambda AtA, Aty: lsq_linear(AtA, Aty, bounds=(lb, ub), method='bvls')
        parseResult = lambda result: result.x

    elif linearConstrained and nonNegativeOnly:
        # Non-negative linear LSQ
        if nnlsSolver == 'fnnls':
            linSolver = lambda AtA, Aty: fnnls(AtA, Aty, tol=tol, maxiter=maxiter)
        elif nnlsSolver == 'nnlsbpp':
            linSolver = lambda AtA, Aty: nnlsbpp(AtA, Aty, np.linalg.solve(AtA, Aty))
        elif nnlsSolver == 'cvx':
            linSolver = lambda AtA, Aty: cvxnnls(AtA, Aty, tol=tol, maxiter=maxiter)
        parseResult = lambda result: result
    
    # Ensure correct formatting and shield against float-point errors
    validateResult = lambda result: np.maximum(lb,np.minimum(ub,np.atleast_1d(result)))
    # ----------------------------------------------------------
    return axis, L, linSolver, parseResult, validateResult, includeRegularization
# ===========================================================================================

# ===========================================================================================
def _model_evaluation(ymodels,parfit,paruq,uq):
    """
    Model evaluation
    ================

    Evaluates the model(s) response(s) at the fitted solution and progates the uncertainty in the fit
    parameters to the model response. 
    """
    modelfit, modelfituq = [],[]
    for ymodel in ymodels: 
        modelfit.append(ymodel(parfit))
        if uq: 
            modelfituq.append(paruq.propagate(ymodel))
        else:
            modelfituq.append(UQResult('void'))
    return modelfit, modelfituq
# ===========================================================================================

# ===========================================================================================
def _goodness_of_fit_stats(ys,yfits,noiselvl,nParam):
    """
    Evaluation of goodness-of-fit statistics
    ========================================

    Evaluates a series of statistical criteria for the goodness-of-fit of a list of dataset
    and returns a list of dictionaries with the statistics for each dataset.
    """
    stats = []
    for y,yfit,sigma in zip(ys,yfits,noiselvl):
        Ndof = len(y) - nParam
        stats.append(goodness_of_fit(y, yfit, Ndof, sigma))
    return stats 
# ===========================================================================================

def _penalty_augmentation(alpha,L,P,type):
# ===========================================================================================
    """ 
    LSQ residual and Jacobian augmentation
    =======================================

    Augments the residual and the Jacobian of a LSQ problem to include the
    regularization penalty. The residual and Jacobian contributions of the 
    specific regularization methods are analytically introduced. 
    """
    # Compute the regularization penalty augmentation for the residual and the Jacobian
    resreg = L@P
    Jreg = L

    # Include regularization parameter
    resreg = alpha*resreg
    Jreg = alpha*Jreg
    if type=='residual':
        return resreg
    if type=='Jacobian':
        return Jreg
# ===========================================================================================

# ===========================================================================================
def _optimize_scale(y,yfit,subsets):
    """
    Dataset scale optimization
    ==========================

    For two datasets, find the least-squares scaling factor that fits one to each other.
    """
    scales, scales_vec = [], np.zeros_like(yfit) 
    for subset in subsets:
        yfit_,y_ = (np.atleast_2d(y[subset]) for y in [yfit, y]) # Rescale the subsets corresponding to each signal
        scale = np.squeeze(np.linalg.lstsq(yfit_.T,y_.T,rcond=None)[0])
        scales.append(scale) # Store the optimized scales of each signal
        scales_vec[subset] = scale 
    return scales, scales_vec
# ===========================================================================================



# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================



def nlls(y, model, par0, lb=None, ub=None, weights=None, noiselvl=None,
                 multistart=1, tol=1e-10, maxiter=3000, extrapenalty=None,
                 fitscale=True, uq=True, covmatrix=None):
    r""" 
    Non-linear least squares (NLLS) solver

    Fits a non-linear model `y(\theta)` to the data by solving the following (penalized) non-linear least squares problem: 

    .. math::

        & \min_{\theta} \left\{ \Vert y - y(\theta) \Vert^2  + \mathcal{P}(\theta) \right\} \\
        & \text{subject to} \\
        & \theta_\mathrm{lb} \leq \theta \leq \theta_\mathrm{ub}

    Penalty terms `\mathcal{P}(\theta)` can be included if specified.

    Parameters
    ----------
    y : array_like or list of array_like
        Input dataset(s) to be fitted.
    
    model : callable
        Function taking an array of parameters and returning the dipolar signal(s) as an array or a list of arrays.
    
    par0  : array_like
        Start values of the model parameters.
    
    lb : array_like, optional      
        Lower bounds for the model parameters. If not specified, it is left unbounded.
    
    ub : array_like, optional     
        Upper bounds for the model parameters. If not specified, it is left unbounded.
    
    noiselvl : array_like, optional
        Noise standard deviation of the input signal(s), if not specified it is estimated automatically. 

    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting. 
        If not specified all datasets are weighted inversely proportional to their noise levels.
    
    extrapenalty: callable or list thereof
        Custom penalty function(s) to impose upon the solution. A single penalty must be specified as a callable function. 
        Multiple penalties can be specified as a list of callable functons. Each function must take a vector of parameters
        and return a vector to be added to the residual vector (``pen = fcn(param)``).  
        The square of the penalty is computed internally.

    fitscale : boolean, optional
        Enable/disable fitting of the signal scale, by default it is enabled.
    
    multistart : scalar, optional
        Number of starting points for global optimization, the default is 1.
    
    uq : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.    
    
    tol : scalar, optional
        Optimizer function tolerance, the default is 1e-10.
    
    maxiter : scalar, optional
        Maximum number of optimizer iterations, the default is 3000.
    
    covmatrix : array_like with shape(n,n), optional
        Covariance matrix of the noise in the dataset(s). If not specified it is automatically computed.


    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    param : ndarray
        Fitted model parameters.
    
    model : ndarray or list of ndarrays
        Fitted model.

    paramUncert : :ref:`UQResult`
        Covariance-based uncertainty quantification of the fitted parameters.
 
    modelUncert : :ref:`UQResult`
        Covariance-based uncertainty quantification of the fitted model.

    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    
    plot : callable
        Function to display the results. It will 
        display the fitted signals. The function returns the figure object 
        (``matplotlib.figure.Figure``) object as output, which can be 
        modified. Using ``fig = plot(show=False)`` will not render
        the figure unless ``display(fig)`` is called. 
    
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
        fit = dl.nlls(V,Vmodel,par0,lb,ub)


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
        fit = dl.nlls([V1,V2],Vmodel,par0,lb,ub)

    """
    y, model, weights, subsets, noiselvl = parse_multidatasets(y, model, weights, noiselvl)
    Nsignals = len(subsets)
    scales = [1]*Nsignals

    if not np.all(np.isreal(y)):
        raise ValueError('The input signal(s) cannot be complex-valued.')

    lb,ub,boundedParams = _check_bounds(lb,ub,par0=par0)

    includeExtrapenalty = extrapenalty is not None
    if includeExtrapenalty:
        extrapenalty = np.atleast_1d(extrapenalty)
        for penalty in extrapenalty:
            if not callable(penalty):
                raise TypeError("The keyword argument 'extrapenalty' must be a callable function or a list thereof.")

    # Preprare multiple start global optimization if requested
    if multistart>1 and boundedParams:
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

        # Evaluate the model
        ysim = model(p)

        # Check if there are invalid values...
        if any(np.isnan(ysim)) or any(np.isinf(ysim)):
            res = np.zeros_like(ysim) # ...can happen when Jacobian is evaluated outside of bounds
            return res

        # If requested, fit the scale of the signal via linear LSQ
        scales_vec = np.ones_like(y)
        if fitscale:
            scales,scales_vec = _optimize_scale(y,ysim,subsets)
        
        # Compute the residual
        res = weights*(y-scales_vec*ysim)

        # Compute residual from user-defined penalties
        if includeExtrapenalty:
            for penalty in extrapenalty:
                penres = penalty(p)
                penres = np.atleast_1d(penres)
                res = np.concatenate((res,penres))
        
        return res
    # -------------------------------------------------------------------------------    
    fvals = []
    parfits = []
    sols = []
    for par0 in multistarts_par0:
        # Solve the non-linear least squares (NLLS) problem 
        sol = least_squares(lsqresiduals,par0, bounds=(lb,ub), max_nfev=int(maxiter), ftol=tol, method='dogbox')
        sols.append(sol)
        parfits.append(sol.x)
        fvals.append(2*sol.cost) # least_squares uses 0.5*sum(residual**2)          

    # Find global minimum from multiple runs
    globmin = np.argmin(fvals)
    parfit = parfits[globmin]
    sol = sols[globmin]

    # Compute residual vector
    residuals = lsqresiduals(parfit)

    # Calculate parameter confidence intervals
    if uq:
                
        J = Jacobian(lsqresiduals,parfit,lb,ub)

        # Estimate the heteroscedasticity-consistent covariance matrix
        if covmatrix is None:
            # Use estimated data covariance matrix
            covmatrix = hccm(J,residuals,'HC1')
        else:
            # Use user-given data covariance matrix
            covmatrix = hccm(J,covmatrix)
        # Construct confidence interval structure
        paramuq = UQResult('covariance',parfit,covmatrix,lb,ub)
    else:
        paramuq = UQResult('void')

    # Evaluate fitted model and propagate the uncertainty
    def ymodel(n): return lambda p: scales[n]*model(p)[subsets[n]]
    ymodels = [ymodel(n) for n in range(Nsignals)]
    modelfit,modelfituq = _model_evaluation(ymodels,parfit,paramuq,uq)

    # Make lists of data and fits
    ys = [y[subset] for subset in zip(subsets)]
    yfits = modelfit.copy()

    # Calculate goodness of fit
    stats = _goodness_of_fit_stats(ys,yfits,noiselvl,len(par0))

    # Get plot function
    def plotfcn(show=False):
        fig = _plot(ys,yfits,show)
        return fig

    if Nsignals==1: 
        stats = stats[0]
        scales = scales[0]
        fvals = fvals[0]
        modelfit = modelfit[0]
        modelfituq = modelfituq[0]

    return FitResult(
            param=parfit, model=modelfit, paramUncert=paramuq, modelUncert=modelfituq, scale=scales, stats=stats, cost=fvals,
            plot=plotfcn, residuals=residuals, success=sol.success)


# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================




def rlls(y, A, lb=None, ub=None, regparam='aic', regop=None, solver='cvx', reg='auto',
                weights=None, noiselvl=None, uq=True, tol=None, maxiter=None):
    r"""
    Regularized linear least-squares (RLLS) solver

    Fits a set of linear parameters `x` to the data by solving the following regularized linear least squares problem: 

    .. math::

        & \min_{\theta} \left\{ \Vert y - Ax \Vert^2  + \mathcal{R}(x) \right\} \\
        & \text{subject to} \\
        & \theta_\mathrm{lb} \leq x \leq \theta_\mathrm{ub}

    Different regularization terms `\mathcal{R}(x)` can be specified (by default uses Tikhonov regularization).


    Parameters 
    ----------
    y : array_like or list of array_like  
        Input dataset(s) to be fitted.
    
    A : 2D-array_like or list of 2D-array_like
        Design or regressor matrix. 
    
    lb : array_like, optional
        Lower bounds for the linear parameters, assumed unconstrained if not specified.

    ub : array_like, optional
        Upper bounds for the linear parameters, assumed unconstrained if not specified.
    
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
    
    noiselvl : array_like, optional
        Noise standard deviation of the input signal(s), if not specified it is estimated automatically. 

    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting. 
        If not specified all datasets are weighted inversely proportional to their noise levels.
    
    regop : 2D array_like, optional
        Regularization operator matrix, the default is the second-order differential operator.

    solver : string, optional
        Optimizer used to solve the non-negative least-squares problem: 

        * ``'cvx'`` - Optimization of the NNLS problem using cvxopt
        * ``'fnnls'`` - Optimization using the fast NNLS algorithm.
        * ``'nnlsbpp'`` - Optimization using the block principal pivoting NNLS algorithm.
        
        The default is ``'cvx'``.

    uq : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled. 

    tol : scalar, optional 
        Tolerance value for convergence of the NNLS algorithm. If not specified (``None``), the value is set automatically to ``tol = max(K.T@K.shape)*norm(K.T@K,1)*eps``.
        It is only valid for the ``'cvx'`` and ``'fnnls'`` solvers, it has no effect on the ``'nnlsbpp'`` solver.
        
    maxiter: scalar, optional  
        Maximum number of iterations before termination. If not specified (``None``), the value is set automatically to ``maxiter = 5*shape(K.T@K)[1]``.
        It is only valid for the ``'cvx'`` and ``'fnnls'`` solvers, it has no effect on the ``'nnlsbpp'`` solver.

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    param : ndarray
        Fitted linear parameters.
    
    model : ndarray 
        Fitted model.

    paramUncert : :ref:`UQResult`
        Uncertainty quantification of the linear parameters.

    modelUncert : :ref:`UQResult`
        Uncertainty quantification of the fitted model.
     
    regparam : float int
        Regularization parameter used in the optimization
    
    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    
    plot : callable
        Function to display the results. It will 
        display the fitted signals and distance distributions with
        confidence intervals. The function returns the figure object 
        (``matplotlib.figure.Figure``) object as output, which can be 
        modified. Using ``fig = plot(show=False)`` will not render
        the figure unless ``display(fig)`` is called. 
    
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
    y, A, weights, subsets, noiselvl, prescales = dl.utils.parse_multidatasets(y, A, weights, noiselvl, precondition=True)

    Ndatasets = len(subsets)
    N = np.shape(A)[1]

    # Check and validate boundaries
    lb,ub,_ = _check_bounds(lb,ub,N=N)

    # Determine optimal choice of solver and settings
    ax, L, linSolver, parseResult, validateResult, includeRegularization = _prepare_linear_lsq(A, lb, ub, reg, regop, tol, maxiter, solver)

    # Determine an optimal value of the regularization parameter if requested
    if includeRegularization:
        if type(regparam) is str:
            alpha = dl.selregparam(y, A, ax, regparam, L=L, weights=weights, noiselvl=noiselvl)
        else:
            alpha = regparam
    else: alpha=None

    # Prepare components of the LSQ problem
    AtAreg, Aty = _lsqcomponents(y, A, L, alpha, weights)

    # Parameter estimation
    # ----------------------------------------------------------------
    # Solve the linear least-squares problem
    result = linSolver(AtAreg, Aty)
    result = parseResult(result)
    xfit = validateResult(result)

    # Get fit final status
    success = ~np.all(xfit==0)

    # Construct residual parts for for the residual and regularization terms
    res = weights*(y - A@xfit)

    # Construct Jacobians for the residual and penalty terms
    J = A*weights[:,np.newaxis]
    if  includeRegularization:
        res = np.concatenate((res,_penalty_augmentation(alpha,L,xfit,'residual')))
        J = np.concatenate((J,_penalty_augmentation(alpha,L,xfit,'Jacobian')))

    # Get objective function value
    fval = np.linalg.norm(res)**2

    # Uncertainty quantification
    # ----------------------------------------------------------------
    if uq:
        covmat = hccm(J,res,'HC1')
        paramuq = UQResult('covariance',xfit,covmat,lb,ub)
    else:
        paramuq = UQResult('void')

    # Get fits and their uncertainty
    # --------------------------------------
    ymodels = [lambda x: prescales[n]*(A@x)[subsets[n]] for n in range(Ndatasets)]
    def ymodel(n): return lambda x: prescales[n]*(A@x)[subsets[n]]
    ymodels = [ymodel(n) for n in range(Ndatasets)]
    modelfit,modelfituq = _model_evaluation(ymodels,xfit,paramuq,uq)
    
    # Make lists
    ys = [prescale*y[subset] for prescale,subset in zip(prescales,subsets)]
    yfits = modelfit.copy()

    # Goodness-of-fit
    # --------------------------------------
    H = A@(np.linalg.pinv(AtAreg)@A.T)
    nParam = np.trace(H)
    stats = _goodness_of_fit_stats(ys,yfits,noiselvl,nParam)

    # Results display
    # ---------------------------------------

    def plotfcn(show=False):
        fig = _plot(ys,yfits,show)
        return fig 

    # Do not return lists if there is a single dataset
    scales,stats,modelfit,modelfituq = [var[0] if len(subsets)==1 else var for var in [prescales,stats,modelfit,modelfituq]]

    return FitResult(param=xfit, model=modelfit, paramUncert=paramuq, modelUncert=modelfituq, regparam=alpha, scale=scales, stats=stats, 
                     plot=plotfcn, cost=fval, residuals=res, success=success)
# ===========================================================================================




# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================

def snlls(y, Amodel, par0, lb=None, ub=None, lbl=None, ubl=None, nnlsSolver='cvx', reg='auto', weights=None,
          regparam='aic', multistart=1, regop=None, alphareopt=1e-3, extrapenalty=None,
          nonlin_tol=1e-9, nonlin_maxiter=1e8, lin_tol=1e-15, lin_maxiter=1e4, noiselvl=None,
          uq=True, fitscale=True):
    r""" Separable non-linear least squares (SNLLS) solver

    Fits a linear set of parameters `\theta_\mathrm{lin}` and non-linear parameters `\theta_\mathrm{nonlin}`
    by solving the following general penalized separable non-linear least-squares problem: 

    .. math::

        & \quad\quad \min_{\theta_\mathrm{nonlin}} \left\{ \Vert y - A(\theta_\mathrm{nonlin})\theta_\mathrm{lin}(\theta_\mathrm{nonlin}) \Vert^2  + \mathcal{P}(\theta_\mathrm{nonlin},\theta_\mathrm{lin}) \right\} \\
        & \text{with} \\
        & \quad\quad \theta_\mathrm{lin}(\theta_\mathrm{nonlin}) = {\arg\!\min}_{\theta} \left\{ \Vert y - A(\theta_\mathrm{nonlin})\theta \Vert^2 +  \mathcal{R}(\theta) \right\} \\
        & \text{subject to} \\
        & \quad\quad \theta_\mathrm{lb} \leq \theta_\mathrm{nonlin} \leq \theta_\mathrm{ub} \quad\quad \theta_\mathrm{lbl} \leq \theta_\mathrm{lin} \leq \theta_\mathrm{ubl} 

    If the non-linear function `A(\theta_\mathrm{nonlin})` yields an ill-conditioned problem, the solver will include a regularization penalty `\mathcal{R}(\theta)` (Tikhonov by default) 
    with automated selection of the regularization parameter. Additional penalty terms `\mathcal{P}(\theta_\mathrm{nonlin},\theta_\mathrm{lin})` can be included if specified.

    Parameters
    ----------
    y : array_like or list of array_like
        Input dataset(s) to be fitted.
        
    Amodel : callable
        Function taking an array of non-linear parameters and
        returning a matrix array or a list thereof.

    par0 : array_like
        Start values of the non-linear parameters.

    lb : array_like, optional
        Lower bounds for the non-linear parameters, assumed unconstrained if not specified.

    ub : array_like, optional
        Upper bounds for the non-linear parameters, assumed unconstrained if not specified.

    lbl : array_like, optional
        Lower bounds for the linear parameters, assumed unconstrained if not specified.

    ubl : array_like, optional
        Upper bounds for the linear parameters, assumed unconstrained if not specified.

    reg : boolean or string, optional
        Determines the use of regularization on the solution of the linear problem.
        
        * ``'auto'`` - Automatic decision based con the condition number of the non-linear model ``Amodel``.
        * ``True`` - Forces regularization regardless of the condition number
        * ``False`` - Disables regularization regardless of the condition number
        
        The default is ``'auto'``.

    regparam : string or float scalar, optional
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
        
        The regularization parameter can be manually specified by passing a scalar value
        instead of a string. The default ``'aic'``.

    regop : 2D array_like, optional
        Regularization operator matrix, the default is the second-order differential operator.

    extrapenalty: callable or list thereof
        Custom penalty function(s) to impose upon the solution. A single penalty must be specified as a callable function. 
        Multiple penalties can be specified as a list of callable functons. Each function must take two inputs, a vector of non-linear parameters
        and a vector of linear parameters, and return a vector to be added to the residual vector (``pen = fcn(pnonlin,plin)``).  
        The square of the penalty is computed internally.

    alphareopt : float scalar, optional
        Relative parameter change threshold for reoptimizing the regularization parameter
        when using a selection method, the default is 1e-3.

    fitscale : boolean, optional
        Enable/disable fitting of the signal scale, by default it is enabled.

    nnlsSolver : string, optional
        Solver used to solve a non-negative least-squares problem (if applicable):

        * ``'cvx'`` - Optimization of the NNLS problem using the cvxopt package.
        * ``'fnnls'`` - Optimization using the fast NNLS algorithm.
        * ``'nnlsbpp'`` - Optimization using the block principal pivoting NNLS algorithm.
        
        The default is ``'cvx'``.

    noiselvl : array_like, optional
        Noise standard deviation of the input signal(s), if not specified it is estimated automatically. 

    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting.
        If not specified all datasets are weighted inversely proportional to their noise levels.

    multistart : int scalar, optional
        Number of starting points for global optimization, the default is 1.
        
    nonlin_maxiter : float scalar, optional
        Non-linear solver maximal number of iterations, the default is 1e8.

    nonlin_tol : float scalar, optional
        Non-linear solver function tolerance, the default is 1e-9.

    lin_maxiter : float scalar, optional
        Linear solver maximal number of iterations, the default is 1e4.

    lin_tol : float scalar, optional
        Linear solver function tolerance, the default is 1e-15.

    uq : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.


    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    nonlin : ndarray
        Fitted non-linear parameters.
    lin : ndarray
        Fitted linear parameters.
    modelfit : ndarray
        Fitted model.
    nonlinUncert : :ref:`UQResult`
        Uncertainty quantification of the non-linear parameter set.
    linUncert : :ref:`UQResult`
        Uncertainty quantification of the linear parameter set.
    modelUncert : :ref:`UQResult`
        Uncertainty quantification of the fitted model.
    regparam : scalar
        Regularization parameter value used for the regularization of the linear parameters.
    plot : callable
        Function to display the results. It will display the fitted data.
        The function returns the figure object (``matplotlib.figure.Figure``)
        object as output, which can be modified. Using ``fig = plot(show=False)`` 
        will not render the figure unless ``display(fig)`` is called. 
    stats : dict
        Goodness of fit statistical estimators

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
    # Ensure that all arrays are numpy.nparray
    par0 = np.atleast_1d(par0)

    # Parse multiple datsets and non-linear operators into a single concatenated vector/matrix
    y, Amodel, weights, subsets, noiselvl, prescales = dl.utils.parse_multidatasets(y, Amodel, weights, noiselvl, precondition=True)
    Ndatasets = len(subsets)    
    
    # Get info on the problem parameters and non-linear operator
    A0 = Amodel(par0)
    Nnonlin = len(par0)
    Nlin = np.shape(A0)[1]

    includeExtrapenalty = extrapenalty is not None
    if includeExtrapenalty:
        extrapenalty = np.atleast_1d(extrapenalty)
        for penalty in extrapenalty:
            if not callable(penalty):
                raise TypeError("The keyword argument 'extrapenalty' must be a callable function or a list thereof.")

    # Checks for bounds constraints
    # ----------------------------------------------------------
    lb,ub,nonLinearBounded = _check_bounds(lb,ub,par0=par0)
    lbl,ubl,linearBounded = _check_bounds(lbl,ubl,N=Nlin)

    # Prepare the optimal solver setup for the linear problem
    ax, L, linSolver, parseResult, validateResult, includeRegularization = _prepare_linear_lsq(A0,lbl,ubl,reg,regop,lin_tol,lin_maxiter,nnlsSolver)

    # Pre-allocate nonlocal variables
    check = False
    regparam_prev = 0
    par_prev = [0]*len(par0)
    alpha = None
    linfit = np.zeros(Nlin)
    scales = [1 for _ in subsets]

    def linear_problem(A,optimize_alpha,alpha):
    #===========================================================================
        """
        Linear problem
        ------------------
        Solves the linear subproblem of the SNLLS objective function via linear LSQ 
        constrained, unconstrained, or regularized.
        """
        
        # Optimiza the regularization parameter only if needed
        if optimize_alpha:
            alpha = dl.selregparam(y, A, ax, regparam, weights=weights, L=L)

        # Components for linear least-squares
        AtA, Aty = _lsqcomponents(y, A, L, alpha, weights=weights)
         
        # Solve the linear least-squares problem
        result = linSolver(AtA, Aty)
        result = parseResult(result)
        linfit = validateResult(result)

        return linfit, alpha
    #===========================================================================

    def ResidualsFcn(p):
    #===========================================================================
        """
        Residuals function
        ------------------
        Provides vector of residuals, which is the objective function for the
        non-linear least-squares solver. 
        """

        nonlocal par_prev, check, regparam_prev, linfit, alpha, scales

        # Non-linear model evaluation
        A = Amodel(p)

        # Check whether optimization of the regularization parameter is needed
        if includeRegularization :
            if type(regparam) is str:
                # If the parameter vector has not changed by much...
                if check and all(abs(par_prev-p)/p < alphareopt):
                    # ...use the alpha optimized in the previous iteration
                    optimize_alpha = False
                    alpha = regparam_prev
                else:
                    # ...otherwise optimize with current settings
                    alpha = regparam
                    optimize_alpha = True
                    check = True
            else:
                # Fixed regularization parameter
                alpha = regparam
                optimize_alpha = False
            # Store current iteration data for next one
            par_prev = p
            regparam_prev = alpha
        else:
            # Non-linear operator without penalty
            optimize_alpha = False
            alpha = 0

        linfit,alpha = linear_problem(A,optimize_alpha,alpha)
        regparam_prev = alpha

        # Evaluate full model residual
        yfit = A@linfit

        # Optimize the scale yfit
        if fitscale:
            scales, scales_vec = _optimize_scale(y,yfit,subsets)

        # Compute residual vector
        res = weights*(scales_vec*(Amodel(p)@linfit) - y)

        # Compute residual from user-defined penalties
        if includeExtrapenalty:
            for penalty in extrapenalty:
                penres = penalty(p,linfit)
                penres = np.atleast_1d(penres)
                res = np.concatenate((res,penres))

        if includeRegularization:
            # Augmented residual
            res_reg = _penalty_augmentation(alpha, L, linfit,'residual')
            res = np.concatenate((res,res_reg))
        
        return res
    #===========================================================================

    # Preprare multiple start global optimization if requested
    if multistart > 1 and not nonLinearBounded:
        raise TypeError('Multistart optimization cannot be used with unconstrained non-linear parameters.')
    multiStartPar0 = dl.utils.multistarts(multistart, par0, lb, ub)

    # Pre-allocate containers for multi-start run
    fvals, nonlinfits, linfits, sols = ([] for _ in range(4))

    # Multi-start global optimization
    for par0 in multiStartPar0:
        # Run the non-linear solver
        sol = least_squares(ResidualsFcn, par0, bounds=(lb, ub), max_nfev=int(nonlin_maxiter), ftol=nonlin_tol)
        nonlinfits.append(sol.x)
        linfits.append(linfit)
        fvals.append(2*sol.cost) # least_squares uses 0.5*sum(residual**2)          
        sols.append(sol)

    # Find global minimum from multiple runs
    globmin = np.argmin(fvals)
    linfit = linfits[globmin]
    nonlinfit = nonlinfits[globmin]
    sol = sols[globmin]
    Afit = Amodel(nonlinfit)

    scales_vec = np.zeros_like(y) 
    for subset,scale in zip(subsets,scales): 
        scales_vec[subset] = scale 
    yfit = scales_vec*(Afit@linfit)

    # Uncertainty analysis
    #---------------------
    if uq:

        def uq_subset(uq_full,subset,subset_lb,subset_ub):
        #-----------------------------------------------------------------------------
            "Get the uncertainty quantification for a subset of parameters"

            subset_model = lambda x: x[subset]
            uq_subset = uq_full.propagate(subset_model,lbm=subset_lb, ubm=subset_ub)

            return uq_subset
        #-----------------------------------------------------------------------------
        
        # Compute the fit residual
        res = ResidualsFcn(nonlinfit)
        
        # Jacobian (non-linear part)
        Jnonlin = Jacobian(ResidualsFcn,nonlinfit,lb,ub)
        # Jacobian (linear part)
        Jlin = scales_vec[:,np.newaxis]*Amodel(nonlinfit)
        if includeExtrapenalty:
            for penalty in extrapenalty:
                Jlin = np.concatenate((Jlin, Jacobian(lambda plin: penalty(nonlinfit,plin),linfit,lbl,ubl)))
        if includeRegularization:
            Jlin = np.concatenate((Jlin, _penalty_augmentation(alpha, L, linfit,'Jacobian')))
        # Full Jacobian
        J = np.concatenate((Jnonlin,Jlin),axis=1)

        # Calculate the heteroscedasticity consistent covariance matrix
        covmatrix = hccm(J, res, 'HC1')

        # Get combined parameter sets and boundaries
        parfit = np.concatenate((nonlinfit, linfit))
        lbs = np.concatenate((lb, lbl))
        ubs = np.concatenate((ub, ubl))

        # Construct the uncertainty quantification object
        paramuq = UQResult('covariance', parfit, covmatrix, lbs, ubs)

        # Split the uncertainty quantification of nonlinear/linear parts
        nonlin_subset = np.arange(0,Nnonlin)
        lin_subset = np.arange(Nnonlin,Nnonlin+Nlin)
        paramuq_nonlin = uq_subset(paramuq,nonlin_subset,lb,ub)
        paramuq_lin = uq_subset(paramuq,lin_subset,lbl,ubl)

    else:
        paramuq_nonlin = UQResult('void')
        paramuq_lin = UQResult('void')
        paramuq = UQResult('void')

    for i in range(len(subsets)):
        scales[i] *= prescales[i]

    # Get fitted signals and their uncertainty
    parfit = np.concatenate((nonlinfit, linfit))
    nonlin_idx = np.arange(len(nonlinfit))
    lin_idx = np.arange(len(nonlinfit),len(parfit))
    def ymodel(n): return lambda p: scales[n]*((Amodel(p[nonlin_idx])@p[lin_idx])[subsets[n]]) 
    ymodels = [ymodel(n) for n in range(len(subsets))]
    modelfit,modelfituq = _model_evaluation(ymodels,parfit,paramuq,uq)

    # Make lists of data and fits
    ys = [prescales[n]*y[subset] for n,subset in enumerate(subsets)]
    yfits = modelfit.copy()

    # Goodness-of-fit
    # ---------------
    stats = _goodness_of_fit_stats(ys,yfits,noiselvl,Nnonlin)
    
    # Display function
    def plotfcn(show=False):
        fig = _plot(ys,yfits,show)
        return fig

    if len(stats) == 1: 
        stats = stats[0]
        fvals = fvals[0]
        modelfit = modelfit[0]
        modelfituq = modelfituq[0]

    return FitResult(nonlin=nonlinfit, lin=linfit, model=modelfit, nonlinUncert=paramuq_nonlin,
                     linUncert=paramuq_lin, modelUncert=modelfituq, regparam=alpha, plot=plotfcn,
                     stats=stats, cost=fvals, residuals=sol.fun, success=sol.success, scale=scales)
# ===========================================================================================



# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================
# ===========================================================================================





def fnnls(AtA, Atb, tol=None, maxiter=None, verbose=False):
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
    if tol is None:
        eps = np.finfo(float).eps
        tol = 10*eps*np.linalg.norm(AtA,1)*max(np.shape(AtA))
    if maxiter is None:
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
        w[passive] = -np.inf
        if verbose:
            print(f'{outIteration:10.0f}{iIteration:15.0f}{max(w):20.4e}\n')

    if verbose:
        if unsolvable:
            print('Optimization stopped because the solution cannot be further changed. \n')
        elif any(~passive):
            print('Optimization stopped because the active set has been completely emptied. \n')
        elif w>tol:
            print(f'Optimization stopped because the gradient (w) is inferior than the tolerance value TolFun = {tol:.6e}. \n')
        else:
            print('Solution found. \n')
    
    return x

#=====================================================================================




def cvxnnls(AtA, Atb, tol=None, maxiter=None):
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
    if tol is None:
        eps = np.finfo(float).eps
        tol = 10*eps*np.linalg.norm(AtA,1)*max(np.shape(AtA))
    if maxiter is None:
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


def nnlsbpp(AtA, AtB, x0=None):
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

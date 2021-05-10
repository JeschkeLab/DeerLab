# snlls.py - Separable non-linear least-squares solver
# ---------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import copy
import numpy as np
from scipy.optimize import least_squares, lsq_linear
import matplotlib.pyplot as plt
from numpy.linalg import solve

# Import DeerLab depencies
import deerlab as dl
from deerlab.utils import goodness_of_fit, hccm, isempty, Jacobian
from deerlab.nnls import cvxnnls, fnnls, nnlsbpp
from deerlab.classes import UQResult, FitResult

def snlls(y, Amodel, par0, lb=None, ub=None, lbl=None, ubl=None, nnlsSolver='cvx', reg='auto', weights=None,
          regtype='tikhonov', regparam='aic', multistart=1, regorder=2, alphareopt=1e-3, extrapenalty=None,
          nonlin_tol=1e-9, nonlin_maxiter=1e8, lin_tol=1e-15, lin_maxiter=1e4, huberparam=1.35,
          uq=True):
    r""" Separable Non-linear Least Squares Solver

    Parameters
    ----------
    y : array_like or list of array_like
        Input data to be fitted.
        
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

    regType : string, optional
        Regularization penalty type:

        * ``'tikhonov'`` - Tikhonov regularizaton
        * ``'tv'``  - Total variation regularization
        * ``'huber'`` - Huber regularization
        
        The default is ``'tikhonov'``.

    regorder : int scalar, optional
        Order of the regularization operator

    regParam : string or float scalar, optional
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

    extrapenalty: callable 
        Custom penalty function to impose upon the solution. Must take two inputs, a vector of non-linear parameters
        and a vector of linear parameters, and return a vector to be added to the residual vector (``pen = fcn(pnonlin,plin)``).  
        The square of the penalty is computed internally.

    alphareopt : float scalar, optional
        Relative parameter change threshold for reoptimizing the regularization parameter
        when using a selection method, the default is 1e-3.

    nnlsSolver : string, optional
        Solver used to solve a non-negative least-squares problem (if applicable):

        * ``'cvx'`` - Optimization of the NNLS problem using the cvxopt package.
        * ``'fnnls'`` - Optimization using the fast NNLS algorithm.
        * ``'nnlsbpp'`` - Optimization using the block principal pivoting NNLS algorithm.
        
        The default is ``'cvx'``.

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
        Fitted non-linear parameters
    lin : ndarray
        Fitted linear parameters
    modelfit : ndarray
        Fitted model
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


    Notes
    -----
    Fits a linear set of parameters ``plin`` and non-linear parameters ``pnlin``
    by solving the following non-linear least squares problem

    .. code-block:: text

        [pnlin,plin] = argmin ||Amodel(pnlin)*plin - y||^2
                       subject to  pnlin in [lb,ub]
                                   plin  in [lbl,ubl]


    where the parameter space is composed of a set of non-linear parameters pnlin
    and linear parameters ``plin``. If the non-linear function Amodel yields an
    ill-conditioned problem, the solver will include a regularization penalty
    (Tikhonov by default) and solve in the Tikhonov case the following problem

    .. code-block:: text

        [pnlin,plin] = argmin ||Amodel(pnlin)*plin - y||^2 + alpha^2*||L*plin||^2
                       subject to  pnlin in [lb,ub]
                                   plin  in [lbl,ubl]


    where ``alpha`` and ``L`` are the regularization parameter and operator, respectively.

    When solving the linear problem th function will
    identify and adapt automatically to the following scenarios:

        * Well-conditioned + unconstrained       ``plin = solve(A,y)``
        * Well-conditioned + constrained         ``plin = lsqlin(A,y,lb,ub)``
        * Ill-conditioned  + unconstrained       ``plin = solve(AtA + alpha^2*LtL, Aty)``
        * Ill-conditioned  + constrained         ``plin = lsqlin(AtA + alpha^2*LtL,Kty,lb,ub)``
        * Ill-conditioned  + non-negativity      ``plin = fnnls((AtA + alpha^2*LtL),Aty)``

    By default, for poorly conditioned cases, Tikhonov regularization with
    automatic AIC regularization parameter selection is used.

    """
    # Ensure that all arrays are numpy.nparray
    par0 = np.atleast_1d(par0)

    # Parse multiple datsets and non-linear operators into a single concatenated vector/matrix
    y, Amodel, weights, subsets, prescales = dl.utils.parse_multidatasets(y, Amodel, weights, precondition=True)

    # Get info on the problem parameters and non-linear operator
    A0 = Amodel(par0)
    Nnonlin = len(par0)
    Nlin = np.shape(A0)[1]
    linfit = np.zeros(Nlin)
    scales = [1 for _ in subsets]

    # Determine whether to use regularization penalty
    illConditioned = np.linalg.cond(A0) > 10
    if reg == 'auto':
        includeRegularization  = illConditioned
    else:
        includeRegularization  = reg

    # Checks for bounds constraints
    # ----------------------------------------------------------
    if lb is None or isempty(lb):
        lb = np.full(Nnonlin, -np.inf)

    if ub is None or isempty(ub):
        ub = np.full(Nnonlin, np.inf)

    if lbl is None or isempty(lbl):
        lbl = np.full(Nlin, -np.inf)

    if ubl is None or isempty(ubl):
        ubl = np.full(Nlin, np.inf)

    lb, ub, lbl, ubl = np.atleast_1d(lb, ub, lbl, ubl)

    # Check that the correct number of boundaries are given
    if len(lb) != Nnonlin or len(ub) != Nnonlin:
        raise TypeError('The lower/upper bounds of the non-linear problem must have ', Nnonlin, ' elements')
    if len(lbl) != Nlin or len(ubl) != Nlin:
        raise TypeError('The lower/upper bounds of the linear problem must have ', Nlin, ' elements')

    # Check that the boundaries are valid
    if np.any(ub < lb) or np.any(ubl < lbl):
        raise ValueError('The upper bounds cannot be larger than the lower bounds.')
    # Check that the non-linear start values are inside the box constraint
    if np.any(par0 > ub) or np.any(par0 < lb):
        raise ValueError('The start values are outside of the specified bounds.')
    # ----------------------------------------------------------


    # Check if the nonlinear and linear problems are constrained
    nonLinearConstrained = (not np.all(np.isinf(lb))) or (not np.all(np.isinf(ub)))
    linearConstrained = (not np.all(np.isinf(lbl))) or (not np.all(np.isinf(ubl)))
    # Check for non-negativity constraints on the linear solution
    nonNegativeOnly = (np.all(lbl == 0)) and (np.all(np.isinf(ubl)))

    # Use an arbitrary axis
    ax = np.arange(1, Nlin+1)
    if includeRegularization :
        # Get regularization operator
        regorder = np.minimum(Nlin-1, regorder)
        L = dl.regoperator(ax, regorder)
    else:
        L = None

    # Prepare the linear solver
    # ----------------------------------------------------------
    if not linearConstrained:
        # Unconstrained linear LSQ
        linSolver = solve
        parseResult = lambda result: result

    elif linearConstrained and not nonNegativeOnly:
        # Constrained linear LSQ
        linSolver = lambda AtA, Aty: lsq_linear(AtA, Aty, bounds=(lbl, ubl), method='bvls')
        parseResult = lambda result: result.x

    elif linearConstrained and nonNegativeOnly:
        # Non-negative linear LSQ
        if nnlsSolver == 'fnnls':
            linSolver = lambda AtA, Aty: fnnls(AtA, Aty, tol=lin_tol, maxiter=lin_maxiter)
        elif nnlsSolver == 'nnlsbpp':
            linSolver = lambda AtA, Aty: nnlsbpp(AtA, Aty, np.linalg.solve(AtA, Aty))
        elif nnlsSolver == 'cvx':
            linSolver = lambda AtA, Aty: cvxnnls(AtA, Aty, tol=lin_tol, maxiter=lin_maxiter)
        parseResult = lambda result: result
    # ----------------------------------------------------------

    # Containers for alpha-update checks
    check = False
    regparam_prev = 0
    par_prev = [0]*len(par0)
    alpha = None

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
            alpha = dl.selregparam(y, A, ax, regtype, regparam, weights=weights, regorder=regorder)

        # Components for linear least-squares
        AtA, Aty = dl.lsqcomponents(y, A, L, alpha, weights=weights, regtype=regtype)
         
        # Solve the linear least-squares problem
        result = linSolver(AtA, Aty)
        linfit = parseResult(result)
        linfit = np.atleast_1d(linfit)
        
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

        nonlocal par_prev, check, regparam_prev, scales, linfit, alpha

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
        scales, scales_vec = [], np.zeros_like(yfit) 
        for subset in subsets:
            yfit_,y_ = (np.atleast_2d(y[subset]) for y in [yfit, y]) # Rescale the subsets corresponding to each signal
            scale = np.squeeze(np.linalg.lstsq(yfit_.T,y_.T,rcond=None)[0])
            scales.append(scale) # Store the optimized scales of each signal
            scales_vec[subset] = scale 

        # Compute residual vector
        res = weights*(scales_vec*(Amodel(p)@linfit) - y)

        # Compute residual from custom penalty
        if callable(extrapenalty):
            penres = extrapenalty(p,linfit)
            penres = np.atleast_1d(penres)
            res = np.concatenate((res,penres))

        if includeRegularization:
            # Augmented residual
            res_reg, _ = reg_penalty(regtype, alpha, L, linfit, huberparam, Nnonlin)
            res = np.concatenate((res,res_reg))
        
        return res
    #===========================================================================

    # Preprare multiple start global optimization if requested
    if multistart > 1 and not nonLinearConstrained:
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
        # Compute the fit residual
        res = ResidualsFcn(nonlinfit)
        
        # Jacobian (non-linear part)
        Jnonlin = Jacobian(ResidualsFcn,nonlinfit,lb,ub)

        # Jacobian (linear part)
        Jlin = scales_vec[:,np.newaxis]*Amodel(nonlinfit)
        if callable(extrapenalty):
            Jlin = np.concatenate((Jlin, Jacobian(lambda plin: extrapenalty(nonlinfit,plin),linfit,lbl,ubl)))
        if includeRegularization:
            Jlin = np.concatenate((Jlin, reg_penalty(regtype, alpha, L, linfit, huberparam, Nnonlin)[1]))

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


    for i in range(len(subsets)):
        scales[i] *= prescales[i]

    # Get fitted signals and their uncertainty
    parfit = np.concatenate((nonlinfit, linfit))
    nonlin_idx = np.arange(len(nonlinfit))
    lin_idx = np.arange(len(nonlinfit),len(parfit))
    modelfit, modelfituq = [],[]
    for i,subset in enumerate(subsets): 
        subset_model = lambda p: scales[i]*((Amodel(p[nonlin_idx])@p[lin_idx])[subset])
        modelfit.append(subset_model(parfit))
        if uq: 
            modelfituq.append(paramuq.propagate(subset_model))
        else:
            modelfituq.append(UQResult('void'))

    # Goodness-of-fit
    # ---------------
    stats = []
    for subset in subsets:
        Ndof = len(y[subset]) - Nnonlin
        stats.append(goodness_of_fit(y[subset], yfit[subset], Ndof))
    if len(stats) == 1: 
        stats = stats[0]
        fvals = fvals[0]
        modelfit = modelfit[0]
        modelfituq = modelfituq[0]

    # Display function
    def plotfcn(show=False):
        fig = _plot(subsets,y,yfit,show)
        return fig

    return FitResult(nonlin=nonlinfit, lin=linfit, model=modelfit, nonlinUncert=paramuq_nonlin,
                     linUncert=paramuq_lin, modelUncert=modelfituq, regparam=alpha, plot=plotfcn,
                     stats=stats, cost=fvals, residuals=sol.fun, success=sol.success, scale=scales)
# ===========================================================================================


def uq_subset(uq_full,subset,subset_lb,subset_ub):
#===========================================================================
    "Get the uncertainty quantification for a subset of parameters"

    subset_model = lambda x: x[subset]
    uq_subset = uq_full.propagate(subset_model,lbm=subset_lb, ubm=subset_ub)

    return uq_subset
#===========================================================================

def reg_penalty(regtype, alpha, L, x, eta, Nnonlin):
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
        resreg = L@x
        Jreg = L
    elif regtype == 'tv':
        resreg =((L@x)**2 + eps)**(1/4)
        Jreg = 2/4*((( ( (L@x)**2 + eps)**(-3/4) )*(L@x))[:, np.newaxis]*L)
    elif regtype == 'huber':
        resreg = np.sqrt(np.sqrt((L@x/eta)**2 + 1) - 1)
        Jreg = 0.5/(eta**2)*((((np.sqrt((L@x/eta)**2 + 1) - 1 + eps)**(-1/2)*(((L@x/eta)**2 + 1+ eps)**(-1/2)))*(L@x))[:, np.newaxis]*L)

    # Include regularization parameter
    resreg = alpha*resreg
    Jreg = alpha*Jreg

    return resreg, Jreg
# ===========================================================================================


def _plot(subsets,y,yfit,show):
# ===========================================================================================
    nSignals = len(subsets)
    fig,axs = plt.subplots(nSignals,figsize=[7,3*nSignals])
    axs = np.atleast_1d(axs)
    for i in range(nSignals): 
        subset = subsets[i]
        # Plot the experimental signal and fit
        axs[i].plot(y[subset],'.',color='grey',alpha=0.5)
        axs[i].plot(yfit[subset],'tab:blue')
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


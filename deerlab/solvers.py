# solvers.py - Collection of least-squares solvers
# --------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

# External dependencies
import numpy as np
import cvxopt as cvx
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, lsq_linear
# DeerLab dependencies
import deerlab as dl
from deerlab.classes import UQResult, FitResult
from deerlab.utils import multistarts, hccm, parse_multidatasets, goodness_of_fit, Jacobian, isempty
import time 
from functools import partial

def timestamp():
# ===========================================================================================
    # Get the seconds since epoch
    secondsSinceEpoch = time.time()
    # Convert seconds since epoch to struct_time
    timeObj = time.localtime(secondsSinceEpoch)
    # get the current timestamp elements from struct_time object i.e.
    return '[%d-%d-%d %d:%d:%d]' % ( timeObj.tm_mday, timeObj.tm_mon, timeObj.tm_year, timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
# ===========================================================================================

def _plot(ys,yfits,yuqs,axis=None,xlabel=None):
# ===========================================================================================
    """
    Plot method for the FitResult object
    ====================================

    Plots the input dataset(s), their fits, and uncertainty bands.
    """

    complexy = np.any([np.iscomplex(y).any() for y in ys])

    nSignals = len(ys)
    ncols = 2 if complexy else 1
    fig,axs = plt.subplots(nSignals,ncols,figsize=[7*ncols,4*nSignals])
    axs = np.atleast_1d(axs)

    if axis is None: 
        axis = [np.arange(len(y)) for y in ys]
    if not isinstance(axis,list): 
        axis = [axis]
    axis = [np.real(ax) for ax in axis]
    if xlabel is None: 
        xlabel = 'Array elements'


    n = 0
    for i,(y,yfit,yuq) in enumerate(zip(ys,yfits,yuqs)): 
        # Plot the experimental signal and fit
        axs[n].plot(axis[i],y.real,'.',color='grey')
        axs[n].plot(axis[i],yfit.real,color='#4550e6')
        if yuq.type!='void': 
            axs[i].fill_between(axis[i],yuq.ci(95)[:,0].real,yuq.ci(95)[:,1].real,alpha=0.4,linewidth=0,color='#4550e6')
        axs[n].set_xlabel(xlabel)
        axs[n].set_ylabel(f'Dataset #{i+1}')
        axs[n].legend(('Data (real)','Fit','95%-CI'),loc='best',frameon=False)
        n += 1
        
        if complexy: 
            axs[n].plot(axis[i],y.imag,'.',color='grey')
            axs[n].plot(axis[i],yfit.imag,color='tab:orange')
            if yuq.type!='void': 
                axs[n].fill_between(axis[i],yuq.ci(95)[:,0].imag,yuq.ci(95)[:,1].imag,alpha=0.4,color='tab:orange',linewidth=0)
            axs[n].set_xlabel(xlabel)
            axs[n].set_ylabel(f'Dataset #{i+1}')
            axs[n].legend(('Data (imag)','Fit','95%-CI'),loc='best',frameon=False)
            n += 1

    plt.tight_layout()
    plt.autoscale(enable=True, axis='both', tight=True)

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
    if par0 is not None:
        if len(lb) != len(par0) or len(ub) != len(par0):
            raise TypeError('The lower/upper bounds and start values must have the same number of elements')
                # Check that the non-linear start values are inside the box constraint
        if np.any(par0 > ub) or np.any(par0 < lb):
            raise ValueError('The start values are outside of the specified bounds.')
    if N is not None:
        if len(lb) != N or len(ub) != N:
            raise TypeError(f'The lower/upper bounds and start values must have {N} elements')
    if np.any(ub < lb):
        raise ValueError('The upper bounds cannot be larger than the lower bounds.')

    # Check if the nonlinear problem is bounded
    isbounded = (not np.all(np.isinf(lb))) or (not np.all(np.isinf(ub)))

    return lb,ub,isbounded
# ===========================================================================================

# ===========================================================================================
def _check_frozen(frozen,N=None,par0=None):
    """
    Frozen parameters parsing and validation 
    ========================================

    Parses the input frozen conditions and determines whether they are correct.
    """
    if par0 is not None:
        N = len(par0)

    if frozen is None: 
        frozen = np.full(N,False)
        xfrozen = np.full(N,None)
    else:
        frozen = np.atleast_1d(frozen.copy())
        xfrozen = frozen.copy()
        frozen[np.where(frozen!=None)[0]] = True
        frozen[np.where(frozen==None)[0]] = False
        frozen = frozen.astype(bool)
    return frozen,xfrozen
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
    Kw = weights[:,np.newaxis]*K
    Vw = weights*V
    KtK = Kw.T@Kw
    KtV = Kw.T@Vw
    
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
        d = np.minimum(2,len(axis))
        L = dl.regoperator(axis,d)


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
def _unfrozen_subset(param,frozen,parfrozen):
    param,frozen,parfrozen = np.atleast_1d(param,frozen,parfrozen)
    param = param.tolist()
    # Account for frozen parameters
    param = np.atleast_1d([parfrozen[n] if frozen[n] else param.pop(0) for n in range(len(frozen))])
    return param
# ===========================================================================================

# ===========================================================================================
def _unfrozen_subset_inv(param,frozen):
    param,frozen = np.atleast_1d(param,frozen)
    # Account for frozen parameters
    return np.atleast_1d([param[n] for n in range(len(frozen)) if not frozen[n] ])
# ===========================================================================================


# ===========================================================================================
def _insertfrozen(parfit,parfrozen,frozen):
    _parfit = parfrozen.copy()
    _parfit[~frozen] = parfit
    return _parfit.astype(float)
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


def snlls(y, Amodel, par0=None, lb=None, ub=None, lbl=None, ubl=None, nnlsSolver='cvx', reg='auto', weights=None, verbose=0,
          regparam='aic', regparamrange=None, multistart=1, regop=None, alphareopt=1e-3, extrapenalty=None, subsets=None,
          ftol=1e-8, xtol=1e-8, max_nfev=1e8, lin_tol=1e-15, lin_maxiter=1e4, noiselvl=None, lin_frozen=None, mask=None,
          nonlin_frozen=None, uq=True):
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

    regparamrange : array_like, optional 
        Search range for the optimization of the regularization parameter. Must be specified as a list ``[regparam_lb,regparam_ub]`` 
        with the lower/upper boundaries of the regularization parameter. 

    regop : 2D array_like, optional
        Regularization operator matrix, the default is the second-order differential operator.

    extrapenalty: callable or list thereof, optional
        Custom penalty function(s) to impose upon the solution. A single penalty must be specified as a callable function. 
        Multiple penalties can be specified as a list of callable functons. Each function must take two inputs, a vector of non-linear parameters
        and a vector of linear parameters, and return a vector to be added to the residual vector (``pen = fcn(pnonlin,plin)``).  
        The square of the penalty is computed internally.

    alphareopt : float scalar, optional
        Relative parameter change threshold for reoptimizing the regularization parameter
        when using a selection method, the default is ``1e-3``.

    nnlsSolver : string, optional
        Solver used to solve a non-negative least-squares problem (if applicable):

        * ``'cvx'`` - Optimization of the NNLS problem using the cvxopt package.
        * ``'fnnls'`` - Optimization using the fast NNLS algorithm.
        
        The default is ``'cvx'``.

    noiselvl : array_like, optional
        Noise standard deviation of the input signal(s), if not specified it is estimated automatically. 

    mask : arrau_like or list thereof, optional
        Array (or list of arrays) containing boolean (True/False) values defining a mask over one or multiple datasets.
        All dataset point with enabled mask values (True) will be acoounted for during the fitting procedure. All disabled
        values (False) will not be inlcuded. This does not affect model evaluation and uncertainty propagation.   

    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting.
        If not specified all datasets are weighted inversely proportional to their noise levels.

    multistart : int scalar, optional
        Number of starting points for global optimization, the default is ``1``.

    xtol : float scalar, optional
        Tolerance for termination by the change of the independent variables. Default is 1e-8. The optimization process is stopped when ``norm(dx) < xtol * (xtol + norm(x))``.
        If set to ``None``, the termination by this condition is disabled.

    ftol : float scalar, optional
        Tolerance for termination by the change of the cost function. Default is 1e-8. The optimization process is stopped when ``dF < ftol*F``, 
        and there was an adequate agreement between a local quadratic model and the true model in the last step.
        If set to ``None``, the termination by this condition is disabled.
        
    max_nfev : float scalar, optional
        Maximum number of function evaluations before the termination. the default is ``1e8``.

    lin_maxiter : float scalar, optional
        Linear solver maximal number of iterations, the default is ``1e4``.

    lin_tol : float scalar, optional
        Linear solver tolerance to decide when a value is defined as a zero, the default is ``1e-15``.

    nonlin_frozen : array_like 
        Values for non-linear parameters to be frozen during the optimization. If set to ``None`` a non-linear parameter 
        will be optimized, if set to a numerical value it will be fixed to it and ignored during the optimization.

    lin_frozen : array_like 
        Values for linear parameters to be frozen during the optimization. If set to ``None`` a linear parameter 
        will be optimized, if set to a numerical value it will be fixed to it and ignored during the optimization.

    uq : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.

    subsets : array_like, optional 
        Vector of dataset subset indices, if the datasets are passed concatenated instead of a list.

    verbose : scalar integer, optional
        Level of verbosity during the analysis:

            * ``0`` : Work silently (default).
            * ``1`` : Display progress including the non-linear least-squares' solver termination report.
            * ``2`` : Display progress including the non-linear least-squares' solver iteration details.

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    nonlin : ndarray
        Fitted non-linear parameters.
    lin : ndarray
        Fitted linear parameters.
    param : ndarray
        Fitted full parameter set. 
    model : ndarray
        Fitted model.
    nonlinUncert : :ref:`UQResult`
        Uncertainty quantification of the non-linear parameter set.
    linUncert : :ref:`UQResult`
        Uncertainty quantification of the linear parameter set.
    paramUncert : :ref:`UQResult`
        Uncertainty quantification of the full parameter set.
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
        
    if verbose>0: 
        print(f'{timestamp()} Preparing the SNLLS analysis...')

    
    # Ensure that all arrays are numpy.nparray
    par0 = np.atleast_1d(par0)
    
    # Parse multiple datsets and non-linear operators into a single concatenated vector/matrix
    y, Amodel, weights, mask, subsets, noiselvl = parse_multidatasets(y, Amodel, weights, noiselvl, masks=mask, subsets=subsets)    
    
    if not callable(Amodel):
        Amatrix = Amodel.copy()
        Amodel = lambda _: Amatrix
        par0 = np.array([])

    if regparamrange is None: 
        regparamrange = [1e-8,1e3]

    # Get info on the problem parameters and non-linear operator
    A0 = Amodel(par0)
    if len(np.shape(A0))!=2:
        Amodel_ = Amodel
        Amodel = lambda p: np.atleast_2d(Amodel_(p)).T
        A0 = Amodel(par0)

    if np.shape(A0)[0]!=np.shape(y)[0]:
        raise RuntimeError(f"The number of elements ({np.shape(A0)[0]}) in the model's output does not match the number of elements ({np.shape(y)[0]}) in the data.")

    complexy = np.iscomplexobj(y)
    imagsubsets = []
    if complexy: 
        # If the data is complex-valued
        for subset in subsets:
            imagsubsets.append(subset + len(y))
        y = np.concatenate([y.real,y.imag])
        weights = np.concatenate([weights,weights])   
        mask = np.concatenate([mask,mask])   
    Amodel__ = Amodel
    if np.iscomplexobj(A0):
       # If the design matrix is complex-valued
        Amodel = lambda p: np.concatenate([Amodel__(p).real,Amodel__(p).imag]) 
        A0 = np.concatenate([A0.real,A0.imag]) 
        A0 = Amodel(par0)
        if not complexy: 
            # If the data is not complex-valued
            y = np.concatenate([y,np.zeros_like(y)])  
            weights = np.concatenate([weights,weights])
            mask = np.concatenate([mask,mask])      
    elif complexy:
        # If the design matrix is not complex-valued, but the data is
        Amodel = lambda p: np.concatenate([Amodel__(p),np.zeros_like(A0)]) 


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

    nonlin_frozen,nonlin_parfrozen = _check_frozen(nonlin_frozen,par0=par0)
    lin_frozen,lin_parfrozen = _check_frozen(lin_frozen,N=Nlin)

    Nlin_notfrozen = Nlin - np.sum(lin_frozen)
    Nnonlin_notfrozen = Nnonlin - np.sum(nonlin_frozen)

    lb_red,ub_red,par0_red = [var[~nonlin_frozen] for var in [lb,ub,par0]]
    lbl_red,ubl_red = [var[~lin_frozen] for var in [lbl,ubl]]
    A0red = A0[:,~lin_frozen]
    # Redefine model to take just the unfrozen parameter subset

    _Amodel = Amodel
    Amodel = lambda param: _Amodel(_unfrozen_subset(param,nonlin_frozen,nonlin_parfrozen))

    if includeExtrapenalty:
        extrapenalty_ = [penalty for penalty in extrapenalty]
        extrapenalty = [lambda pnonlin, plin: penalty(_unfrozen_subset(pnonlin,nonlin_frozen,nonlin_parfrozen), plin) for penalty in extrapenalty_]

    # Prepare the optimal solver setup for the linear problem

    if Nlin_notfrozen>0:
        _, L, linSolver, parseResult, validateResult, includeRegularization = _prepare_linear_lsq(A0red,lbl_red,ubl_red,reg,regop,lin_tol,lin_maxiter,nnlsSolver)
    else: 
        includeRegularization = False

    # Pre-allocate nonlocal variables
    check = False
    regparam_prev = 0
    par_prev = [0]*len(par0_red)
    alpha = None
    xfit = np.zeros(Nlin)
    Ndof_lin = 0

    if verbose>0: 
        print(f'{timestamp()} Preparations completed.')


    def linear_problem(y,A,optimize_alpha,alpha):
    #===========================================================================
        """
        Linear problem
        ------------------
        Solves the linear subproblem of the SNLLS objective function via linear LSQ 
        constrained, unconstrained, or regularized.
        """

        # If all linear parameters are frozen, do not optimize 
        if Nlin_notfrozen==0:
            return lin_parfrozen.astype(float), None, 0

        # Remove columns corresponding to frozen linear parameters 
        Ared = A[:,~lin_frozen]
        # Frozen component of the model response
        yfrozen = (A[:,lin_frozen]@lin_parfrozen[lin_frozen]).astype(float)

        # Optimiza the regularization parameter only if needed
        if optimize_alpha:
            alpha = dl.selregparam((y-yfrozen)[mask], Ared[mask,:], linSolver, regparam, 
                                        weights=weights[mask], regop=L, candidates=regparamrange, 
                                        noiselvl=noiselvl,searchrange=regparamrange)

        # Components for linear least-squares
        AtA, Aty = _lsqcomponents((y-yfrozen)[mask], Ared[mask,:], L, alpha, weights=weights[mask])

        Ndof = np.maximum(0,np.trace(Ared@np.linalg.pinv(AtA)))

        # Solve the linear least-squares problem
        result = linSolver(AtA, Aty)
        result = parseResult(result)
        xfit = validateResult(result)

        # Insert back the frozen linear parameters
        xfit = _insertfrozen(xfit,lin_parfrozen,lin_frozen)


        return xfit, alpha, Ndof
    #===========================================================================

    def ResidualsFcn(p):
    #===========================================================================
        """
        Residuals function
        ------------------
        Provides vector of residuals, which is the objective function for the
        non-linear least-squares solver. 
        """

        nonlocal par_prev, check, regparam_prev, xfit, alpha, Ndof, Ndof_lin

        # Non-linear model evaluation
        A = Amodel(p)

        # Check whether optimization of the regularization parameter is needed
        if includeRegularization :
            if type(regparam) is str:
                if Nnonlin_notfrozen>0:
                    # If the parameter vector has not changed by much...
                    if check and all(abs(par_prev-p)/(p+np.finfo(np.float64).eps) < alphareopt):
                        # ...use the alpha optimized in the previous iteration
                        optimize_alpha = False
                        alpha = regparam_prev
                    else:
                        # ...otherwise optimize with current settings
                        alpha = regparam
                        optimize_alpha = True
                        check = True
                else:
                    alpha = regparam
                    optimize_alpha = True
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

        xfit,alpha,Ndof_lin = linear_problem(y,A,optimize_alpha,alpha)
        regparam_prev = alpha

        # Compute residual vector
        res = weights*(Amodel(p)@xfit - y)

        # Apply mask to residual
        res = res[mask]

        # Compute residual from user-defined penalties
        if includeExtrapenalty:
            for penalty in extrapenalty:
                penres = penalty(p,xfit)
                penres = np.atleast_1d(penres)
                res = np.concatenate((res,penres))

        if includeRegularization:
            # Augmented residual
            res_reg = _penalty_augmentation(alpha, L, xfit,'residual')
            res = np.concatenate((res,res_reg))
        
        return res
    #===========================================================================

    # -------------------------------------------------------------------
    #  Only linear parameters
    # -------------------------------------------------------------------
    if Nnonlin_notfrozen==0 and Nlin_notfrozen>0:

        if verbose>0: 
            print(f'{timestamp()} Linear least-squares routine in progress...')
    
        if type(regparam) is str and includeRegularization:
            # Optimized regularization parameter
            alpha = regparam
            optimize_alpha = True
        else:
            # Fixed regularization parameter
            alpha = regparam
            optimize_alpha = False

        # Solve the linear LSQ problem
        linfit, alpha, Ndof_lin = linear_problem(y,Amodel(None),optimize_alpha,alpha)
        
        # Compute the residual
        res = ResidualsFcn(None)

        fvals = np.sum(res**2)
        nonlinfit = None

    # -------------------------------------------------------------------
    #  With non-linear parameters
    # -------------------------------------------------------------------
    elif Nnonlin_notfrozen>0:

        if verbose>0: 
            print(f'{timestamp()} Non-linear least-squares routine in progress...')

        # Preprare multiple start global optimization if requested
        if multistart > 1 and not nonLinearBounded:
            raise TypeError('Multistart optimization cannot be used with unconstrained non-linear parameters.')
        multiStartPar0 = multistarts(multistart, par0_red, lb_red, ub_red)

        # Pre-allocate containers for multi-start run
        fvals, nonlinfits, linfits, sols = [],[],[],[]

        # Multi-start global optimization
        for par0 in multiStartPar0:

            # Run the non-linear solver
            sol = least_squares(ResidualsFcn, par0, bounds=(lb_red, ub_red), max_nfev=int(max_nfev), xtol=xtol, ftol=ftol, verbose=verbose)
            nonlinfits.append(sol.x)
            linfits.append(xfit)
            fvals.append(2*sol.cost) # least_squares uses 0.5*sum(residual**2)          
            sols.append(sol)

        # Find global minimum from multiple runs
        globmin = np.argmin(fvals)
        linfit = linfits[globmin]
        nonlinfit = nonlinfits[globmin]
        fvals = np.min(fvals)
            
    # Insert frozen parameters back into the nonlinear parameter vector  
    if nonlinfit is not None: 
        par_prev = nonlinfit.copy()
    nonlinfit = _insertfrozen(nonlinfit,nonlin_parfrozen,nonlin_frozen)

    # Compute the fit residual
    _ResidualsFcn = lambda nonlinfit: ResidualsFcn(_unfrozen_subset_inv(nonlinfit,nonlin_frozen))
    res = _ResidualsFcn(nonlinfit)
    fvals = np.sum(res**2) 

    if verbose>0: 
        print(f'{timestamp()} Least-squares routine finished.')

    # Uncertainty analysis
    #---------------------
    if uq:
        
        if verbose>0: 
            print(f'{timestamp()} Uncertainty analysis in progress...')

        # Define the indices of the parameter subsets
        nonlin_subset = np.arange(0,Nnonlin)
        lin_subset = np.arange(Nnonlin,Nnonlin+Nlin)

        def uq_subset(uq_full,subset,subset_lb,subset_ub):
        #-----------------------------------------------------------------------------
            "Get the uncertainty quantification for a subset of parameters"

            subset_model = lambda x: x[subset]
            uq_subset = uq_full.propagate(subset_model,lb=subset_lb, ub=subset_ub)

            return uq_subset
        #-----------------------------------------------------------------------------

        # Jacobian (non-linear part)
        Jnonlin = Jacobian(_ResidualsFcn,nonlinfit,lb,ub)
        # Jacobian (linear part)
        scale = np.trapz(linfit,np.arange(Nlin))
        Jlin = (weights[:,np.newaxis]*Amodel(nonlinfit))[mask,:]
        if includeExtrapenalty:
            for penalty in extrapenalty:
                Jlin = np.concatenate((Jlin, Jacobian(lambda plin: penalty(nonlinfit,plin),linfit,lbl,ubl)))
        if includeRegularization:
            Jlin = np.concatenate((Jlin, _penalty_augmentation(alpha, L, linfit,'Jacobian')))
        Jlin *= scale
        # Full Jacobian
        J = np.concatenate((Jnonlin,Jlin),axis=1)

        # Calculate the heteroscedasticity consistent covariance matrix
        covmatrix = hccm(J, res)

        # Get combined parameter sets and boundaries
        parfit = np.concatenate((nonlinfit, linfit))
        lbs = np.concatenate((lb, lbl))
        ubs = np.concatenate((ub, ubl))

        # Account for scale of linear parameters since covariance matrix is scale invariant
        covmatrix[np.ix_(lin_subset,lin_subset)] *= scale**2
        
        # Set rows/columns corresponding to frozen parameters to zero
        frozen = np.concatenate([nonlin_frozen,lin_frozen])
        covmatrix[frozen,:] = 0
        covmatrix[:,frozen] = 0

        # Construct the uncertainty quantification object
        paramuq = UQResult('covariance', parfit, covmatrix, lbs, ubs)

        # Split the uncertainty quantification of nonlinear/linear parts
        paramuq_nonlin = uq_subset(paramuq,nonlin_subset,lb,ub)
        paramuq_lin = uq_subset(paramuq,lin_subset,lbl,ubl)

    else:
        paramuq_nonlin = UQResult('void')
        paramuq_lin = UQResult('void')
        paramuq = UQResult('void')

    if verbose>0: 
        print(f'{timestamp()} Uncertainty analysis completed.')
        print(f'{timestamp()} Model evaluation in progress...')

    # Get fitted signals and their uncertainty
    parfit = np.concatenate((nonlinfit, linfit))
    nonlin_idx = np.arange(len(nonlinfit))
    lin_idx = np.arange(len(nonlinfit),len(parfit))
    Amodel = _Amodel # Use the model with the full parameter set
    def ymodel(n):
        return lambda p: (Amodel(p[nonlin_idx])@p[lin_idx])[subsets[n]]
    if complexy: 
        ymodel_ = ymodel
        def ymodel(n):
            return lambda p: ymodel_(n)(p) + 1j*(Amodel(p[nonlin_idx])@p[lin_idx])[imagsubsets[n]]
    ymodels = [ymodel(n) for n in range(len(subsets))]
    modelfit,modelfituq = _model_evaluation(ymodels,parfit,paramuq,uq)

    # Make lists of data and fits
    ys = [y[subset] for subset in subsets]
    if complexy: 
        ys = [ys[n] + 1j*y[imagsubset] for n,imagsubset in enumerate(imagsubsets)]
    yfits = modelfit.copy()
    yuq = modelfituq.copy()

    if verbose>0: 
        print(f'{timestamp()} Model evaluation completed.')

    # Goodness-of-fit
    # ---------------
    Ndof = Nnonlin + Ndof_lin 
    stats = _goodness_of_fit_stats(ys,yfits,noiselvl,Ndof)

    # Display function
    plotfcn = partial(_plot,ys,yfits,yuq)

    if len(stats) == 1: 
        stats = stats[0]
        modelfit = modelfit[0]
        modelfituq = modelfituq[0]

    return FitResult(nonlin=nonlinfit, lin=linfit, param=parfit, model=modelfit, nonlinUncert=paramuq_nonlin,
                     linUncert=paramuq_lin, paramUncert=paramuq, modelUncert=modelfituq, regparam=alpha, plot=plotfcn,
                     stats=stats, cost=fvals, residuals=res)
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
    Fast non-negative least-squares (NNLS) solver.

    Solves the problem 
    
    .. math:: \min \Vert b - Ax \Vert^2 
    
    where ``AtA`` `= A^TA` and ``Atb`` `= A^Tb` using the fast 
    non-negative least-squares (FNNLS) algorithm [1]_, [2]_.

    Parameters
    ----------
    AtA : matrix_like 
        Matrix `A^TA`
    Atb : array_like 
        Vector `A^Tb`
    tol : float scalar 
        Tolerance  used for deciding when elements of x are less than zero.
        The default is ``tol = max(shape(AtA))*norm(AtA,1)*eps``. 
    maxiter : integer scalar 
        Maximum number of iterations. The default is ``maxiter = 5*max(shape(AtA))``
    verbose : boolean
        If set to ``True``, display the iteration details during the optimization.

    Returns
    -------
    x : ndarray 
        Solution vector.

    References
    ----------
    .. [1] R. Bro, S. De Jong,
        A Fast Non-Negativity-Constrained Least Squares Algorithm
        Journal of Chemometrics 11 (1997) 393-401
        
    .. [2] Lawson and Hanson,
        Solving Least Squares Problems
        Prentice-Hall, 1974.

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
    r"""
    Non-negative least-squares (NNLS) via the CVXOPT package.

    Solves the problem 
    
    .. math:: \min \Vert b - Ax \Vert^2 
    
    where ``AtA`` `= A^TA` and ``Atb`` `= A^Tb` using the non-linear convex programming 
    tools from CVXOPT [1]_, [2]_.

    Parameters
    ----------
    AtA : matrix_like 
        Matrix `A^TA`
    Atb : array_like 
        Vector `A^Tb`
    tol : float scalar 
        Tolerance  used for deciding when elements of x are less than zero.
        The default is ``tol = max(shape(AtA))*norm(AtA,1)*eps``. 
    maxiter : integer scalar 
        Maximum number of iterations. The default is ``maxiter = 5*max(shape(AtA))``
    verbose : boolean
        If set to ``True``, display the iteration details during the optimization.
        
    Returns
    -------
    x : ndarray 
        Solution vector.

    References
    ----------
    .. [1] M. Andersen, J. Dahl, and L. Vandenberghen,
       CVXOPT, Python Software for Convex Optimization 
       http://cvxopt.org/
        
    .. [2] Rein, Lewe, Andrade, Kacprzak and Weber
        J. Magn. Reson., 295 (2018) 17–26.
        Global analysis of complex PELDOR time traces

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
    try:
        P = cvx.solvers.qp(cAtA, cAtb, I, lb, initvals=cvx.matrix(x0))['x']
        P = np.squeeze(np.asarray(P))
    except: 
        P = x0
    return P
#=====================================================================================
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2023: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.solvers import snlls
from deerlab.fitresult import FitResult
from deerlab.utils import formatted_table, parse_multidatasets
from deerlab.bootstrap_analysis import bootstrap_analysis
from deerlab.classes import UQResult
from sys import stdout

from deerlab.model import Model
from copy import copy
#--------------------------------------------------------------------------

def _outerOptimization(fitfcn,penalty_objects,sigma):
    """
    (Private function)

    A method to optimize the fit of a model with penalties.

    This method returns a function that can be used to evaluate the fit of a model with penalties. It takes in the following arguments:

    fitfcn : callable
    The function to be optimized, which takes in a set of parameters and returns a scalar value representing the fit of the model.

    penalty_objects : list of Penalty objects
    A list of penalty objects that define the penalty functions to be applied to the fit function. The list can have up to three penalty objects.

    sigma : numpy.ndarray
    The vector of observation uncertainties to be used in the penalty functions.
    Returns

    fitfcn_ : callable
    A function that can be used to evaluate the fit of a model with penalties. This function takes in a set of parameters and returns a scalar value representing the fit of the model.
    """
    # If there are no penalties
    if len(penalty_objects)==0:
        fitfcn_ = lambda y: fitfcn(y,[None])

    # Otherwise, prepare to solve multiobjective problem 
    elif len(penalty_objects)==3:
        thirdfcn = lambda y,*param: penalty_objects[2].optimize(lambda weight: fitfcn(y,[*param,weight]),y,sigma)[1]
        secondfcn = lambda y,*param: penalty_objects[1].optimize(lambda weight: fitfcn(y,[*param,weight,thirdfcn(y,*param,weight)]),y,sigma)[1]
        fitfcn_ = lambda y: penalty_objects[0].optimize(lambda weight: fitfcn(y,[weight,secondfcn(y,weight),thirdfcn(y,weight,secondfcn(weight))]),y,sigma)[0]

    elif len(penalty_objects)==2:
        secondfcn = lambda y,*param: penalty_objects[1].optimize(lambda weight: fitfcn(y,[*param,weight]),y,sigma)[1]
        fitfcn_ = lambda y: penalty_objects[0].optimize(lambda weight: fitfcn(y,[weight,secondfcn(y,weight)]),y,sigma)[0]

    elif len(penalty_objects)==1:
        fitfcn_ = lambda y: penalty_objects[0].optimize(lambda weight: fitfcn(y,[weight]),y,sigma)[0]
    else: 
        raise RuntimeError('The fit() function can only handle up to three penalties.')

    return fitfcn_

#--------------------------------------------------------------------------
def _print_fitresults(fitresult,model):
    """
    (Private function)

    Construct summary table of fit results to print
    
    This helper method takes the output of a model fit and constructs a
    summary table of the goodness-of-fit statistics, along with the 
    hyperparameters of the model (if any). The method is intended to be used
    to print the results of the model fit.
    
    Parameters
    ----------
    fitresult : ``FitResult`` object
        Result of a model fit.
        
    model : ``Model`` object
        Model object instance used for the fit.
        
    Returns
    -------
    table_string : str
        A string containing the summary table of the fit results.
    """

    #-----------------------------------------------------
    def colortxt(str, color, spaces, is_tty=stdout.isatty()):
        """
        (Private function)

        Helper method that applies ANSI codes to add colors to the text
        in a terminal. 

        Parameters
        ----------
        str : str 
            The string to be colored.
        color : str
            The color to be applied to the text. It can be 'red', 'yellow', or 'white'.
        spaces : int
            The number of spaces to add to the beginning and end of the string.
        is_tty : bool
            A boolean value indicating whether the output is a terminal or not.
            If it is not a terminal, the method simply returns the string without
            applying any coloring.

        Returns
        -------
        colored_str : str
            Colored string.
        """
        if color=='red': color = '\033[91m'
        if color=='yellow': color = '\033[93m'
        if color=='white': color = ''
        if is_tty:
            return str
        else: 
            return f"{color}" +" "*spaces + f"{str}"+" "*spaces + "\033[00m" 
    #-----------------------------------------------------

    # Start printout string
    string = ''
    # Get number of models in the fit
    modelfits = fitresult.model
    if not isinstance(modelfits,list):
        modelfits = [modelfits]
    Ndatasets = len(modelfits)

    # Construct table of goodness-of-fit statistics
    table = []
    table.append([f'Dataset','Noise level','Reduced 𝛘2','Residual autocorr.','RMSD']) # Header
    alignment = ['^','^','^','^','^'] # Tab alignment
    stats = np.atleast_1d(fitresult.stats)
    noiselevels = np.atleast_1d(fitresult.noiselvl)
    for n in range(Ndatasets):
        noiselvl = noiselevels[n]
        chi2red = stats[n]['chi2red']
        rmsd = stats[n]['rmsd']
        autocorr = stats[n]['autocorr']
        # Use colored text to warn of very poor fits
        autocorrcolor = lambda str: colortxt(str,'white',7)
        if autocorr>0.5 and autocorr<1:
            # Relatively acceptable autocorrelations (yellow)
            autocorrcolor = lambda str: colortxt(str,'yellow',7)
        elif autocorr>1:
            # Worrisome autocorrelations (red)
            autocorrcolor = lambda str: colortxt(str,'red',7)
        chicolor = lambda str: colortxt(str,'white',3)
        # Standard deviation of reduced 𝛘2 statistic's uncertainty (Gaussian limit)
        chi2red_sigma = np.sqrt(2/len(modelfits[n]))*3 
        if abs(1-chi2red)>3*chi2red_sigma and abs(1-chi2red)<6*chi2red_sigma:
            # Poor case (yellow), 𝛘2 exceeds thrice the expected uncertainty 
            chicolor = lambda str: colortxt(str,'yellow',3)
        elif abs(1-chi2red)>6*chi2red_sigma:
            # Horrible case (red), 𝛘2 exceeds six times the expected uncertainty 
            chicolor = lambda str: colortxt(str,'red',3)
        # Convert numbers to well-formatted strings
        noiselvl,chi2red,autocorr,rmsd = [f'{var:.3f}' if var<1e3 or var>1e-3 else f'{var:.2e}' for var in [noiselvl,chi2red,autocorr,rmsd]] 
        table.append([f'#{1+n}',noiselvl,chicolor(chi2red),autocorrcolor(autocorr),rmsd])
    # Add auto-formatted table string
    string += 'Goodness-of-fit: \n'
    string += formatted_table(table,alignment) + '\n'

    # Construct table of model hyperparameters
    hasregularization = fitresult.regparam!=0
    haspenalties = fitresult.penweights
    if hasregularization or haspenalties:
        string += 'Model hyperparameters: \n'
        tags,values,alignment = [],[],[]
        # If regularization was used, add regularization parameter
        if hasregularization:
            alignment.append('^')
            tags.append('Regularization parameter')      
            regparam = fitresult.regparam
            if regparam is None: regparam = 0  
            values.append(regparam) 
        # If additional penalties were used, add their weights
        if haspenalties:
            for n,penweight in enumerate(fitresult.penweights):
                alignment.append('^')
                tags.append(f'Penalty weight #{n+1}')        
                values.append(penweight) 
        # Format the values
        values = [f'{var:.3f}' if var<1e3 and var>1e-3 else f'{var:.2e}' for var in values] 
        table = [tags,values] 
        # Add to the table
        string += formatted_table(table,alignment) + '\n'

    # Construct table of model parameters fits
    table = []
    table.append([f'Parameter','Value','95%-Confidence interval','Unit','Description']) # Header
    alignment = ['<','<','<','^','<'] # Alignment
    for param in model._parameter_list('vector'):
        if len(np.atleast_1d(getattr(model,param).idx))==1:
            if np.any(getattr(model,param).frozen): 
                # If parameter is frozen, print just the value
                value = getattr(model,param).value
                try:
                    if isinstance(value, (list, tuple, np.ndarray)): value = value[0]
                except: pass
                value = f'{value:.3f}' if abs(value)<1e3 and abs(value)>1e-3 else f'{value:.2e}'
                ci = '(frozen)'
            else:
                # If parameter is scalar, report values and CIs
                value = getattr(fitresult,param)
                if getattr(fitresult,param+'Uncert').type == 'void':
                    ci = ''
                else:
                    ci_lower,ci_upper = getattr(fitresult,param+'Uncert').ci(95)
                    value,ci_lower,ci_upper = [f'{var:.3f}' if abs(var)<1e3 and abs(var)>1e-3 else f'{var:.2e}' for var in [value,ci_lower,ci_upper]]
                    ci = f'({ci_lower},{ci_upper})'
        else:
            # If parameter is vectorial, print just dots
            value = '...'
            ci = '(...,...)'
        unit = str(getattr(model,param).unit)
        description = str(getattr(model,param).description)
        table.append([f'{param}',value,ci,unit,description])
    # Add auto-formatted table string
    string += 'Model parameters: \n'
    string += formatted_table(table,alignment)
    string += '\n'
    return string
#--------------------------------------------------------------------------

def _insert_snlls_optionals_docstrings():
    """
    (Private decorator)
    
    A decorator that takes a function ``func` as input and replaces the 
    string ``'snlls_keyargs_docstrings'`` in the function's docstring with
    the optional keyword arguments documentation for the ``snlls.py`` 
    function. This is done by splitting the ``snlls.py`` docstring into 
    paragraphs, filtering out paragraphs that are already included in the
    outer function's docstring, and concatenating the remaining paragraphs.
    The resulting string is then inserted into ``func``'s docstring and 
    the modified function is returned. This allows for the optional keyword
    arguments documentation to be easily updated in the docstring of any 
    function that uses the ``snlls.py`` function.
    """
    # Get the documentation for the optional keyword arguments in snlls.py also used by fit()
    text = snlls.__doc__
    text = text.split('\n\n')
    # Exclude arguments already set by the outer function
    exclude = ['lb','ub','lbl','ubl','subsets','lin_frozen','nonlin_frozen','regparam','reg','regparamrange', 'extrapenalty']
    paragraphs = [s for s in text if not any(e in s for e in exclude)]
    # Concatenate the arguments
    snlls_keyargs_docs = ''
    for paragraph in paragraphs: 
        # Only keep optional keyword arguments
        if 'optional' in paragraph:
            snlls_keyargs_docs += paragraph + '\n' 

    def decorator(func):
        func.__doc__ = func.__doc__.replace('snlls_keyargs_docstrings',snlls_keyargs_docs)
        return func
    return decorator

#==============================================================================================
@_insert_snlls_optionals_docstrings()
def fit(model_, y, *constants, par0=None, penalties=None, bootstrap=0, noiselvl=None, mask=None, weights=None,
                regparam='aic',reg='auto',regparamrange=None, bootcores=1,**kwargs):
    r"""
    Fit the model(s) to the dataset(s)

    Fit the input model to the data ``y`` via one of the three following approaches: 
    
    - Non-linear least-squares 
    - Regularized linear-least-squares 
    - Separable non-linear least-squares 

    The most appropiate solver is chosen automatically based on the model structure. 

    Parameters
    ----------
    model : :ref:`Model`
        Model object. 

    y : array_like 
        Data to be fitted. 

    par0 : array_like, optional 
        Value at which to initialize the parameter at the start of a fit routine. 
        Must be specified if not defined in the model. Otherwise, it overrides the definition in the model. 

    penalties: callable or list thereof, optional
        Custom penalty function(s) to impose upon the solution. A single penalty must be specified as a callable function. 
        Multiple penalties can be specified as a list of callable functons. Each function must take two inputs, a vector of non-linear parameters
        and a vector of linear parameters, and return a vector to be added to the residual vector (``pen = fcn(pnonlin,plin)``).  
        The square of the penalty is computed internally.

    bootstrap : scalar, optional,
        Bootstrap samples for uncertainty quantification. If ``bootstrap>0``, the uncertainty quantification will be 
        performed via the boostrapping method with based on the number of samples specified as the argument.
    
    bootcores : scalar, optional
        Number of CPU cores/processes for parallelization of the bootstrap uncertainty quantification. If ``cores=1`` no parallel 
        computing is used. If ``cores=-1`` all available CPUs are used. The default is one core (no parallelization).

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
        Search range for the optimization of the regularization parameter. Must be specified as a list ``[regparam_lb, regparam_ub]`` 
        with the lower/upper boundaries of the regularization parameter. The default range is ``[1e-8, 1e3]``. 

    regop : 2D array_like, optional
        Regularization operator matrix, the default is the second-order differential operator.

    alphareopt : float scalar, optional
        Relative parameter change threshold for reoptimizing the regularization parameter
        when using a selection method, the default is ``1e-3``.

    nnlsSolver : string, optional
        Solver used to solve a non-negative least-squares problem (if applicable):

        * ``'qp'`` - Optimization of the NNLS problem using the ``quadprog`` package. Only Python <= 3.10.
        * ``'cvx'`` - Optimization of the NNLS problem using the ``cvxopt`` package.
        * ``'fnnls'`` - Optimization using the fast NNLS algorithm.
        
        The default is ``'cvx'``.

    snlls_keyargs_docstrings

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    <parameter_name> : :ref:`Parameter`
        Fitted value of the <parameter_name> model parameter.
    <parameter_name>Uncert : :ref:`UQResult`
        Uncertainty quantification of the <parameter_name> model parameter.
    param : ndarray
        Fitted parameter vector ordered according to the model parameter indices.
    paramUncert : :ref:`UQResult`
        Uncertainty quantification of the parameter vector ordered according to the model parameter indices.
    paramlist : list
        List of the fitted parameter names ordered according to the model parameter indices.
    model : ndarray
        Fitted model response.     
    regparam : scalar
        Regularization parameter value used for the regularization of the linear parameters.
    penweights : scalar or list thereof 
        Penalty weight value(s) used for the penalties specified through ``penalties``.
    stats : dict
        Goodness of fit statistical estimators

        * ``stats['chi2red']`` - Reduced \chi^2 test
        * ``stats['r2']`` - R^2 test
        * ``stats['rmsd']`` - Root-mean squared deviation (RMSD)
        * ``stats['aic']`` - Akaike information criterion
        * ``stats['aicc']`` - Corrected Akaike information criterion
        * ``stats['bic']`` - Bayesian information criterion
    cost : float
        Value of the cost function at the solution.
    noiselvl : ndarray
        Estimated or user-given noise standard deviations of the individual datasets.
    """

    if not isinstance(model_,Model):
        raise TypeError('The input model must be a valid deerlab.Model object.')
    else:
        model = copy(model_)

    required = len(model._constantsInfo)
    if len(constants)!=required: 
        raise SyntaxError(f'The input model requires {required} constant(s) to be specified. Specify them via fit(model,y,*constants).')
    elif len(constants)>0:
        constants = np.atleast_1d(constants)

    if model.Nlin==0:
        model.addlinear('scale',lb=-np.inf,ub=np.inf,description='Scaling factor')

    normalization = False
    normfactor_keys = []
    for key in model._parameter_list():
        param = getattr(model,key)
        if np.all(param.linear):
            if param.normalization is not None:
                normfactor_key = f'{key}_scale'
                normfactor_keys.append(normfactor_key)
                model.addnonlinear(normfactor_key,lb=-np.inf,ub=np.inf,par0=1,description=f'Normalization factor of {key}')
                getattr(model,normfactor_key).freeze(1)
                normalization = True

    # Get boundaries and conditions for the linear and nonlinear parameters
    ubl,ub = model._split_linear(model._vecsort(model._getvector('ub')))
    lbl,lb = model._split_linear(model._vecsort(model._getvector('lb')))
    frozenl,frozen = model._split_linear(model._vecsort(model._getvector('frozen')))
    valuesl,values = model._split_linear(model._vecsort(model._getvector('value')))

    # Check the initial conditions and whether they are defined
    if par0 is None:
        _,par0 = model._split_linear(model._vecsort(model._getvector('par0')))
    if np.any(par0==None):
        raise RuntimeError(f"It appears some start values (par0) have not been specified. Either specify them in the model definition or using the keyword.")

    linfrozen = np.full(model.Nlin,None)
    linfrozen[frozenl] = valuesl[frozenl]
    nonlinfrozen = np.full(model.Nnonlin,None)
    nonlinfrozen[frozen] = values[frozen]

    if len(linfrozen)==0: linfrozen = [1]

    if type(y) is not list: y = [y]
    ysplit = y.copy()
    y, _, weights, mask, ysubsets, noiselvl = parse_multidatasets(y, None, weights, noiselvl, precondition=False, masks=mask)
    sigmas = np.concatenate([np.full_like(yset,sigma) for sigma,yset in zip(noiselvl,ysplit)])


    if model.Nlin==0 and model.Nnonlin==0:
        raise AssertionError(f'The model has no parameters to fit.')    

    # Get parameter indices in the order spitted out by the solver
    param_idx = [[] for _ in model._parameter_list('vector')]
    idxprev = 0
    for islinear in [False,True]:
        for n,param in enumerate(model._parameter_list('vector')):
            if np.all(getattr(model,param).linear == islinear):
                N = len(np.atleast_1d(getattr(model,param).idx))
                param_idx[n] = np.arange(idxprev,idxprev + N)
                idxprev += N  

    # If there are penalties in the model
    if penalties is not None:
        if not hasattr(penalties, '__iter__'): 
            penalties = [penalties]
        # Get the parameter names of the model
        modelparam = model._parameter_list('vector')
        penaltyfcns = []
        for penalty in penalties:
            # Determine the indices of the subset of parameters the model depends on
            subsets = [getattr(model,modelparam[np.where(np.asarray(modelparam)==param)[0][0]]).idx for param in penalty.signature]
            # Adapt the signature of penaltyfcn for snlls
            penaltyfcns.append(lambda pnonlin,plin,weight: penalty.penaltyfcn(weight,*[np.concatenate([pnonlin,plin])[subset] for subset in subsets]))

        # Prepare the penalties to input to snlls
        extrapenalties = lambda weights: [lambda nonlin,lin: penaltyfcn(nonlin,lin,weight) for penaltyfcn,weight in zip(penaltyfcns,weights)]
    else: 
        # If there are no penalties in the model
        penalties = []
        extrapenalties = lambda *_: None

    # Prepare the separable non-linear least-squares solver
    Amodel_fcn = lambda param: model.nonlinmodel(*constants,*param)
    fitfcn = lambda y,penweights: snlls(y, Amodel_fcn, par0, lb=lb, ub=ub, lbl=lbl, ubl=ubl, mask=mask, weights=weights, 
                                                subsets=ysubsets, lin_frozen=linfrozen, nonlin_frozen=nonlinfrozen,
                                                regparam=regparam, reg=reg, regparamrange=regparamrange, noiselvl=noiselvl,
                                                extrapenalty=extrapenalties(penweights), **kwargs)        

    # Prepare outer optimization of the penalty weights, if necessary
    fitfcn = _outerOptimization(fitfcn,penalties,sigmas)

    # Run the fitting algorithm 
    fitresults = fitfcn(y)
    penweights = [penalty._weight_value for penalty in penalties]

    # If requested, perform a bootstrap analysis
    if bootstrap>0: 
        def bootstrap_fcn(ysim): 
            fit = fitfcn(np.concatenate(ysim))
            if not isinstance(fit.model,list): fit.model = [fit.model]
            return (fit.param,*fit.model)
        # Bootstrapped uncertainty quantification
        if not ('verbose' in kwargs): # Only show boostrap progress bar when verbose is set in kwargs.
            bootstrap_verbose = True
        else:
            bootstrap_verbose = False
            
        param_uq = bootstrap_analysis(bootstrap_fcn,ysplit,fitresults.model,samples=bootstrap,noiselvl=noiselvl,cores=bootcores, verbose=bootstrap_verbose)
        # Include information on the boundaries for better uncertainty estimates
        paramlb = model._vecsort(model._getvector('lb'))[np.concatenate(param_idx)] 
        paramub = model._vecsort(model._getvector('ub'))[np.concatenate(param_idx)] 
        fitresults.paramUncert = UQResult('bootstrap',data=param_uq[0].samples,lb=paramlb,ub=paramub)
        fitresults.param = fitresults.paramUncert.median
        # Get the uncertainty estimates for the model response
        fitresults.model = [param_uq[n].median for n in range(1,len(param_uq))]
        if len(fitresults.model)==1: 
            fitresults.model = fitresults.model[0]
    # Get some basic information on the parameter vector
    keys = model._parameter_list(order='vector')

    # Dictionary of parameter names and fitted values
    FitResult_param = {key : fitvalue if len(fitvalue)>1 else fitvalue[0] for key,fitvalue in zip(keys,[fitresults.param[idx] for idx in param_idx])}
    # Dictionary of parameter names and fit uncertainties
    FitResult_paramuq = {f'{key}Uncert': model._getparamuq(fitresults.paramUncert,idx) for key,idx in zip(keys,param_idx)}
    # Dictionary of other fit quantities of interest
    FitResult_dict = {key: getattr(fitresults,key) for key in ['y','mask','param','paramUncert','model','cost','plot','residuals','stats','regparam','regparam_stats','__plot_inputs']}
    _paramlist = model._parameter_list('vector')

    param_idx = [[] for _ in _paramlist]
    idxprev = 0
    for islinear in [False,True]:
        for n,param in enumerate(_paramlist):
            if np.all(getattr(model,param).linear == islinear):
                N = len(np.atleast_1d(getattr(model,param).idx))
                param_idx[n] = np.arange(idxprev,idxprev + N)
                idxprev += N  

    # Enforce normalization of the linear parameters (if needed) for the final output
    FitResult_param_,FitResult_paramuq_ = FitResult_param.copy(),FitResult_paramuq.copy()
    if normalization:
        for key in keys:
            param = getattr(model,key)
            if key in normfactor_keys:
                param.unfreeze() 
            if np.all(param.linear):
                if param.normalization is not None:
                    def _scale(x):
                        x = x + np.finfo(float).eps
                        return np.mean(x/param.normalization(x))
                    FitResult_param_[f'{key}_scale'] = _scale(FitResult_param_[key]) # Normalization factor
                    FitResult_param_[key] = param.normalization(FitResult_param_[key]) # Normalized value

                    FitResult_paramuq_[f'{key}_scaleUncert'] = FitResult_paramuq_[f'{key}Uncert'].propagate(_scale)
                    FitResult_paramuq_[f'{key}Uncert'] = FitResult_paramuq_[f'{key}Uncert'].propagate(lambda x: x/FitResult_param_[f'{key}_scale'], lb=param.lb, ub=param.ub) # Normalization of the uncertainty
    if len(noiselvl)==1: 
        noiselvl = noiselvl[0]
    
    # Generate FitResult object from all the dictionaries
    fitresult = FitResult({**FitResult_param_,**FitResult_paramuq_, **FitResult_dict,'penweights':penweights,'noiselvl':noiselvl,'paramlist':_paramlist, '_param_idx':param_idx}) 

    fitresult._summary = _print_fitresults(fitresult,model)

    return fitresult
#==============================================================================================



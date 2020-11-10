# fitmultimodel.py - Multi-component distributions SNNLS fit function
# --------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import copy
import numpy as np
import matplotlib.pyplot as plt
from types import FunctionType
import deerlab as dl
from deerlab.utils import hccm, goodness_of_fit, Jacobian
from deerlab.classes import FitResult

def fitmultimodel(V, Kmodel, r, model, maxModels, method='aic', lb=None, ub=None, lbK=None, ubK=None,
                 weights=1, renormalize = True, uqanalysis=True, **kwargs):
    r""" 
    Fits a multi-model parametric distance distribution model to a dipolar signal using separable 
    non-linear least-squares (SNLLS).

    Parameters
    ----------
    V : array_like or list of array_like
        Dipolar signal(s) to be fitted.
    
    Kmodel : callable or 2D-array_like
        Dipolar kernel model. If no kernel parameters must be fitted, it can be specified as a matrix array 
        (or a list thereof if multiple signals are globally fitted). 
        Otherwise, it is a callable function that accepts an array of kernel parameters and returns
        a kernel matrix array or a list thereof.
    
    r : array_like 
        Distance axis, in nanometers.
    
    model : callable 
        Basis component of the multi-component distance distribution. 
        Must be a callable DeerLab model function (e.g. ``dd_gauss`` or ``dd_rice``).
    
    maxModels : scalar 
        Maximal number of components in the multi-component distance distribution.
    
    method : string, optional
        Functional metric used for the selection of the optimal number of components:

        * ``'aic'``  Akaike information criterion
        * ``'aicc'`` corrected Akaike information criterion
        * ``'bic'``  Bayesian information criterion
        * ``'rmsd'`` Root-mean squared deviation
        The default is ``'aic'``.
        
    lb : array_like, optional
        Lower bounds for the distribution basis model parameters. If not specified, parameters are unbounded.
    
    ub : array_like, optional
        Upper bounds for the distribution basis model parameters. If not specified, parameters are unbounded.
    
    ubK : array_like, optional
        Lower bounds for the kernel model parameters. If not specified, parameters are unbounded.
    
    ubK : array_like, optional
        Upper bounds for the kernel model parameters. If not specified, parameters are unbounded.
    
    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting, the default is all weighted equally.
    
    renormalize : boolean, optional
        Enable/disable renormalization of the fitted distribution, by default it is enabled.
    
    uqanalysis : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.  

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    P : ndarray
        Fitted distance distribution with optimal number of components.
    
    Pparam : ndarray
        Fitted distance distribution components parameters
    
    amps : ndarray
        Fitted components amplitudes
    
    Kparam : ndarray
        Fitted kernel parameters.
    
    Puncert : :ref:`UncertQuant`
        Covariance-based uncertainty quantification of the fitted distance distribution
    
    paramUncert : :ref:`UncertQuant`
        Covariance-based uncertainty quantification of the fitted parameters
    
    Nopt : int scalar
        Optimized number of components in model.
    
    Pn : list of ndarrays
        List of all fitted multi-component distance distributions. 
    
    selfun : ndarray
        Selection functional values (as specified as ``method``) for the all fitted multi-component models.
    
    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    
    plot : callable
        Function to display the results. It will display the 
        fitted signals, the distance distribution with confidence intervals, 
        and the values of the selection functional. If requested, the function
        returns the `matplotlib.axes` object as output. 
    
    stats : dict
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


    Notes
    -----
    This function takes advantage of the special structure of a multi-component model, i.e. the separability
    of the component amplitudes as linear parameters from the rest of the non-linear parameters. This makes it
    suitable to be solved as a SNLLS problem.

    Examples
    --------
    A classical example involves the fit of a multi-Gauss distance distribution to a 4-pulse DEER dipolar signal. 
    Since the signal requires additional parameters  (e.g. modulation depth, background parameters,…) a kernel model
    can be defined to account for these::
    
        def Kmodel(Kpar):
            # Unpack parameters
            lam,conc = Kpar
            # Construct kernel
            K = dl.dipolarkernel(t,r,lam,dl.bg_hom3d(t,conc,lam))
            return K
        
        fit = dl.fitmultimodel(V,Kmodel,r,dd_model,Nmax,'aicc')


    If multiple signals are to be fitted globally the example abova can be easily adapted by passing multiple
    signals to the fit function and by returning multiple kernels with the kernel model function::
    
        def Kmodel(Kpar):
            # Unpack parameters
            lam,conc = Kpar
            # Construct kernels for both signals
            K1 = dl.dipolarkernel(t1,r,lam,bg_hom3d(t1,conc,lam))
            K2 = dl.dipolarkernel(t2,r,lam,bg_hom3d(t2,conc,lam))
            return K1,K2
        
        fit = dl.fitmultimodel([V1,V2],Kmodel,r,dd_model,Nmax,'aicc')


    """
    # Ensure that all arrays are numpy.nparray
    r = np.atleast_1d(r)
    
    # Parse multiple datsets and non-linear operators into a single concatenated vector/matrix
    V, Kmodel, weights, Vsubsets, prescales = dl.utils.parse_multidatasets(V, Kmodel, weights,precondition=True)

    # Check kernel model
    if type(Kmodel) is FunctionType:
        # If callable, determine how many parameters the model requires
        nKparam = 0
        notEnoughParam = True
        while notEnoughParam:
            nKparam = nKparam + 1
            try:
                Kmodel(np.random.uniform(size=nKparam))
                notEnoughParam = False
            except ValueError:
                notEnoughParam = True
    else:
        # If the kernel is just a matrix make it a callable without parameters
        nKparam = 0
        K = copy.deepcopy(Kmodel) # need a copy to avoid infite recursion on next step
        Kmodel = lambda _: K
    
    # Extract information about the model
    info = model()
    nparam = len(info['Start'])
    if lb is None:
        lb = info['Lower']
    if ub is None:
        ub = info['Upper']
    paramnames = info['Parameters']

    if lbK is None:
        lbK = []
    if ubK is None:
        ubK = []

    # Ensure that all arrays are numpy.nparray
    lb,ub,lbK,ubK = np.atleast_1d(lb,ub,lbK,ubK)

    if len(lbK) is not nKparam or len(ubK) is not nKparam:
        raise ValueError('The upper/lower bounds of the kernel parameters must be ',nKparam,'-element arrays')

    areCenterDistances = [str in ['Mean','Location'] for str in paramnames]
    if any(areCenterDistances):
        # If the center of the basis function is a parameter limit it 
        # to the distance axis range (stabilizes parameter search)
        ub[areCenterDistances] = max(r)
        lb[areCenterDistances] = min(r)



    def nonlinmodel(par,Nmodels):
    #===============================================================================
        """
        Non-linear augmented kernel model
        ----------------------------------
        This function constructs the actual non-linear function which is 
        passed to the SNLLS problem. The full signal is obtained by multiplication
        of this matrix by a vector of amplitudes. 
        """
        # Get kernel with current non-linear parameters
        K = Kmodel(par[np.arange(0,nKparam)])
        paramsused = nKparam
        Knonlin = np.zeros((K.shape[0],Nmodels))
        for iModel in range(Nmodels):
            subset = np.arange(paramsused, paramsused+nparam)
            paramsused = paramsused + nparam
            # Get basis functions
            Pbasis = model(r,par[subset])
            # Combine all non-linear functions into one
            Knonlin[:,iModel] = K@Pbasis
        return Knonlin
    #===============================================================================

    def Pmodel(nlinpar,linpar):
    #===============================================================================
        """
        Multi-component distribution model
        ----------------------------------
        This function constructs the distance distribution from a set of
        non-linear and linear parameters given certain number of components and
        their basis function.
        """
        paramsused = nKparam
        Pfit = 0
        for amp in linpar:
            subset = np.arange(paramsused, paramsused+nparam)
            paramsused = paramsused + nparam
            # Add basis functions
            Pfit = Pfit + amp*model(r,nlinpar[subset])
        return Pfit
    #===============================================================================

    def logestimators(V,Vfit,plin,pnonlin,functionals):
    #===============================================================================
        """
        Log-Likelihood Estimators
        ---------------------------
        Computes the estimated likelihood of a multi-component model being the
        optimal choice.
        """
        if type(functionals) is not dict:
            functionals = {'rmsd':[],'aic':[],'aicc':[],'bic':[]}
        nParams = len(pnonlin) + len(plin)
        Q = nParams + 1
        N = len(V)
        SquaredSumRes = np.sum((weights*(V - Vfit))**2)
        logprob = N*np.log(SquaredSumRes/N)
        # Compute the estimators
        rmsd = np.sqrt(1/N*SquaredSumRes)
        aic = logprob + 2*Q
        aicc = logprob + 2*Q + 2*Q*(Q+1)/(N-Q-1)
        bic =  logprob + Q*np.log(N)

        # Append results to existing dictionary
        functionals['rmsd'].append(rmsd)
        functionals['aic'].append(aic)
        functionals['aicc'].append(aicc)
        functionals['bic'].append(bic)
        return functionals
    #===============================================================================

    def spread_within_box(lb,ub,Nmodels):
    #===============================================================================
        """
        Start values within boundary box
        ---------------------------------
        Computes the start values of the parameters of the multiple basis functions
        spread equally within the box constraints.
        """
        par = []
        for lbi,ubi in zip(lb,ub):
            par.append(np.linspace(lbi,ubi,Nmodels+2)[1:-1])
        par0 = []        
        for i in range(Nmodels):
            for pari in par:
                par0.append(pari[i])
        return par0
    #===============================================================================

    # Pre-allocate containers
    fits,Vfit,Pfit,plin_,pnonlin_,nlin_ub_,nlin_lb_,lin_ub_,lin_lb_ = ([] for _ in range(9))
    logest = []

    # Loop over number of components in model
    # =======================================
    for Nmodels in np.arange(1,maxModels+1):
        
        # Prepare non-linear model with N-components
        # ===========================================
        Knonlin = lambda par: nonlinmodel(par,Nmodels)
        
        # Start values of non-linear parameters
        par0_K = (ubK - lbK)/2 # start in the middle
        par0_P = spread_within_box(lb,ub,Nmodels) # start spread within boundaries
        par0 = np.concatenate((par0_K, par0_P),axis=None)

        # Box constraints for the model parameters (non-linear parameters)
        nlin_lb = np.tile(lb,(1,Nmodels))
        nlin_ub = np.tile(ub,(1,Nmodels))

        # Add the box constraints on the non-linear kernel parameters
        nlin_lb = np.concatenate((lbK, nlin_lb),axis=None)
        nlin_ub = np.concatenate((ubK, nlin_ub),axis=None)
        
        # Box constraints for the components amplitudes (linear parameters)
        lin_lb = np.ones(Nmodels)        
        lin_ub = np.full(Nmodels,np.inf)
        
        # Separable non-linear least-squares (SNLLS) fit
        scale = 1e2
        fit = dl.snlls(V*scale,Knonlin,par0,nlin_lb,nlin_ub,lin_lb,lin_ub, 
                        weights=weights, reg=False, uqanalysis=False, lin_tol=[], lin_maxiter=[],**kwargs)
        pnonlin = fit.nonlin
        plin = fit.lin

        plin = plin/scale

        # Store the fitted parameters
        pnonlin_.append(pnonlin)
        plin_.append(plin)

        # Get fitted kernel
        Kfit = nonlinmodel(pnonlin,Nmodels)
        
        # Get fitted signal
        Vfit.append(Kfit@plin)
        
        # Get fitted distribution
        Pfit.append(Pmodel(pnonlin,plin))

        # Likelihood estimators
        # =====================
        logest = logestimators(V,Vfit[Nmodels-1],plin,pnonlin,logest)

        # Store other parameters for later
        nlin_ub_.append(nlin_ub)
        nlin_lb_.append(nlin_lb)   
        lin_ub_.append(lin_ub)
        lin_lb_.append(lin_lb-1) # rescale to zero
        fits.append(fit)
    # Select the optimal model
    # ========================
    Peval = Pfit
    fcnals = logest[method]
    idx = np.argmin(fcnals)
    Nopt = idx+1
    fit,Pfit,Vfit,pnonlin,plin = (var[idx] for var in [fits,Pfit,Vfit,pnonlin_,plin_])
    nlin_lb,nlin_ub,lin_lb,lin_ub = (var[idx] for var in [nlin_lb_,nlin_ub_,lin_lb_,lin_ub_])

    # Package the fitted parameters
    # =============================
    fitparam_K = pnonlin[0:nKparam] # Kernel parameters
    fitparam_P = pnonlin[nKparam:nKparam+nparam] # Components parameters
    fitparam_amp = plin # Components amplitudes

    # Uncertainty quantification analysis (if requested)
    # ==================================================
    if uqanalysis:
        # Compute the residual vector
        Knonlin = lambda par: nonlinmodel(par,Nopt)
        res = weights*(Vfit - V)

        lb_full = np.concatenate((nlin_lb, lin_lb))
        ub_full = np.concatenate((nlin_ub, lin_ub))

        # Compute the Jacobian
        Jnonlin = Jacobian(lambda p: weights*(Knonlin(p)@plin),pnonlin,nlin_lb,nlin_ub)
        Jlin = weights[:,np.newaxis]*Knonlin(pnonlin)
        J = np.concatenate((Jnonlin, Jlin),1)

        # Estimate the heteroscedasticity-consistent covariance matrix
        covmatrix = hccm(J,res,'HC1')
        
        # Construct uncertainty quantification structure for fitted parameters
        paramuq = dl.UncertQuant('covariance',np.concatenate((pnonlin, plin)),covmatrix,lb_full,ub_full)
        
        P_subset = np.arange(0, nKparam+nparam*Nopt)
        amps_subset = np.arange(P_subset[-1]+1, P_subset[-1]+1+Nopt)
        Puq = paramuq.propagate(lambda p: Pmodel(p[P_subset],p[amps_subset]), np.zeros(len(r)))
    else: 
        Puq = dl.UncertQuant('void')
        paramuq = dl.UncertQuant('void')

    # Goodness of fit
    stats = []
    for subset in Vsubsets: 
        Ndof = len(V[subset]) - (nKparam + nparam + Nopt)
        stats.append(goodness_of_fit(V[subset],Vfit[subset],Ndof))
    if len(stats)==1: 
        stats = stats[0]

    # If requested re-normalize the distribution
    postscale = np.trapz(Pfit,r)
    if renormalize:
        Pfit = Pfit/postscale
        fitparam_amp = fitparam_amp/sum(fitparam_amp)
        if uqanalysis:
            Puq_ = copy.deepcopy(Puq) # need a copy to avoid infite recursion on next step
            Puq.ci = lambda p: Puq_.ci(p)/postscale

    # Dataset scales
    scales = []
    for i in range(len(Vsubsets)): 
        scales.append(prescales[i]*postscale)
    if len(scales)==1:
        scales = scales[0]

    # Results display function
    plotfcn = lambda: _plot(Vsubsets,V,Vfit,r,Pfit,Puq,fcnals,maxModels,method)

    return FitResult(P=Pfit, Pparam=fitparam_P, Kparam=fitparam_K, amps=fitparam_amp, Puncert=Puq, 
                    paramUncert=paramuq, selfun=fcnals, Nopt=Nopt, Pn=Peval, scale=scales, plot=plotfcn,
                    stats=stats, cost=fit.cost, residuals=fit.residuals, success=fit.success)
# =========================================================================


def _plot(Vsubsets,V,Vfit,r,Pfit,Puq,fcnals,maxModels,method):
# =========================================================================
    nSignals = len(Vsubsets)
    _,axs = plt.subplots(nSignals+1,figsize=[7,3+3*nSignals])
    for i in range(nSignals): 
        subset = Vsubsets[i]
        # Plot the experimental signal and fit
        axs[i].plot(V[subset],'.',color='grey',alpha=0.5)
        axs[i].plot(Vfit[subset],'tab:blue')
        axs[i].grid(alpha=0.3)
        axs[i].set_xlabel('Array Elements')
        axs[i].set_ylabel('V[{}]'.format(i))
        axs[i].legend(('Data','Fit'))

    # Confidence intervals of the fitted distance distribution
    Pci95 = Puq.ci(95) # 95#-confidence interval
    Pci50 = Puq.ci(50) # 50#-confidence interval

    ax = plt.subplot(nSignals+1,2,2*(nSignals+1)-1)
    ax.plot(r,Pfit,color='tab:blue',linewidth=1.5)
    ax.fill_between(r,Pci50[:,0],Pci50[:,1],color='tab:blue',linestyle='None',alpha=0.45)
    ax.fill_between(r,Pci95[:,0],Pci95[:,1],color='tab:blue',linestyle='None',alpha=0.25)
    ax.grid(alpha=0.3)
    ax.legend(['truth','optimal fit','95%-CI'])
    ax.set_xlabel('Distance [nm]')
    ax.set_ylabel('P [nm⁻¹]')
    axs = np.append(axs,ax)

    # Compute the Akaike weights
    dfcnals = fcnals - min(fcnals)
    ax = plt.subplot(nSignals+1,2,2*(nSignals+1))
    ax.bar(np.arange(maxModels)+1,dfcnals + abs(min(dfcnals)),facecolor='tab:blue',alpha=0.6)
    ax.grid(alpha=0.3)
    ax.set_ylabel('Δ{}'.format(method.upper()))
    ax.set_xlabel('Number of components')
    axs = np.append(axs,ax)

    plt.tight_layout()
    plt.show()
    return axs
# =========================================================================


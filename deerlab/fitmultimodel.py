# fitmultimodel.py - Multi-component distributions SNNLS fit function
# --------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import copy
import numpy as np
import matplotlib.pyplot as plt
from types import FunctionType
import deerlab as dl
from deerlab.utils import hccm, goodness_of_fit, Jacobian
from deerlab.classes import FitResult

def fitmultimodel(V, Kmodel, r, model, maxModels, method='aic', lb=None, ub=None, lbK=None, ubK=None,
                 strategy='split', weights=None, renormalize = True, uq=True, tol=1e-9, maxiter=1e8):
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
    
    strategy: string, optional
        Strategy for the initialization of the multi-component non-linear parameters: 

        * ``'spread'`` For each N-component model, the non-linear parameters start values are spread equidistantly over the box constraints. 
                       The number of components are changed in a forward matter, i.e. 1,2,...,N.
        * ``'split'`` For each N-component model, the non-linear parameters start values are selected by splitting the location and spread of 
                      the components obtained from the N-1 component fit. The number of components are changed in a forward matter, i.e. 1,2,...,N.
        * ``'merge'`` For each N-component model, the non-linear parameters start values are selected by merging the location and spread of 
                      the components obtained from the N+1 component fit. The number of components are changed in a backward matter, i.e. N,N-1,...,1.
        
        The default is ``'split'``. 

    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting. 
        If not specified all datasets are weighted inversely proportional to their noise levels.
    
    renormalize : boolean, optional
        Enable/disable renormalization of the fitted distribution, by default it is enabled.
    
    uq : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.  

    tol : scalar, optional 
        Tolerance value for convergence of the NNLS algorithm. If not specified, the value is set to ``tol = 1e-9``.
        
    maxiter: scalar, optional  
        Maximum number of iterations before termination. If not specified, the value is set to ``maxiter = 1e8``.

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
    
    V : ndarray or list thereof
        Fitted dipolar signal(s).

    Puncert : :ref:`UQResult`
        Covariance-based uncertainty quantification of the fitted distance distribution
    
    paramUncert : :ref:`UQResult`
        Covariance-based uncertainty quantification of the fitted parameters

    Vuncert : ndarray or list thereof
        Covariance-based uncertainty quantification of the fitted dipolar signal(s).

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
        and the values of the selection functional. The function returns the figure object 
        (``matplotlib.figure.Figure``) object as output, which can be 
        modified. Using ``fig = plot(show=False)`` will not render
        the figure unless ``display(fig)`` is called. 
    
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
            K = dl.dipolarkernel(t,r,mod=lam,bg=dl.bg_hom3d(t,conc,lam))
            return K
        
        fit = dl.fitmultimodel(V,Kmodel,r,dd_model,Nmax,'aicc')


    If multiple signals are to be fitted globally the example abova can be easily adapted by passing multiple
    signals to the fit function and by returning multiple kernels with the kernel model function::
    
        def Kmodel(Kpar):
            # Unpack parameters
            lam,conc = Kpar
            # Construct kernels for both signals
            K1 = dl.dipolarkernel(t1,r,mod=lam,bg=bg_hom3d(t1,conc,lam))
            K2 = dl.dipolarkernel(t2,r,mod=lam,bg=bg_hom3d(t2,conc,lam))
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
            except:
                notEnoughParam = True
    else:
        # If the kernel is just a matrix make it a callable without parameters
        nKparam = 0
        K = copy.deepcopy(Kmodel) # need a copy to avoid infite recursion on next step
        Kmodel = lambda _: K

    # Extract information about the model
    nparam = len(model.start)
    if lb is None:
        lb = model.lower
    if ub is None:
        ub = model.upper
    paramnames = model.parameters

    if lbK is None:
        lbK = []
    if ubK is None:
        ubK = []

    # Ensure that all arrays are numpy.nparray
    lb,ub,lbK,ubK = np.atleast_1d(lb,ub,lbK,ubK)

    if len(lbK)!=nKparam or len(ubK)!=nKparam:
        raise ValueError('The upper/lower bounds of the kernel parameters must be ',nKparam,'-element arrays')

    areLocations = [str in ['Mean','Location'] for str in paramnames]
    areSpreads   = [str in ['Spread','Width','Standard deviation'] for str in paramnames]
    if any(areLocations):
        # If the center of the basis function is a parameter limit it 
        # to the distance axis range (stabilizes parameter search)
        ub[areLocations] = max(r)
        lb[areLocations] = min(r)



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
        Knonlin = weights[:,np.newaxis]*Knonlin
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

    def initialize_nonlin_param(par,Ncomp,lb,ub,strategy):
    #===============================================================================
        """
        Initialize start values of nonlinear parameters
        ------------------------------------------------
        Computes the start values of the parameters of the multiple basis functions
        according to different strategies.
        """
        # On the first run, use a spread strategy regardless of choice
        if par is None:
            # Spread all components parameters equidistantly within box boundaries
            par = []
            for lbi,ubi in zip(lb,ub):
                par.append(np.linspace(lbi,ubi,Ncomp+2)[1:-1])
            par0 = []        
            for i in range(Ncomp):
                for pari in par:
                    par0.append(pari[i])
            return par0
        else: 
            par = np.asarray(par)

        # Extract information on location and spread of components of previous solution
        Nold = int(len(par)/len(lb))
        areLocations_old = np.tile(areLocations,Nold)
        areSpreads_old = np.tile(areSpreads,Nold)
        locations_old = par[areLocations_old]
        spreads_old = par[areSpreads_old]

        # Pre-allocate containers
        locations_new = np.zeros((Ncomp))
        spreads_new = np.zeros((Ncomp))
        
        # Strategy #1: Merge N+1 components into N components
        if strategy=='merge':
            for n in range(Ncomp):
                locations_new[n] = (locations_old[n] + locations_old[n+1])/2
                spreads_new[n] = (spreads_old[n] + spreads_old[n+1])/2

        # Strategy #2: Split N-1 components into N components
        elif strategy=='split':
            locations_new[0] = locations_old[0] - spreads_old[0]
            spreads_new[0] = spreads_old[0]/2
            for n in range(1,Ncomp-1,1):
                locations_new[n] = (locations_old[n-1] + spreads_old[n-1] + locations_old[n] - spreads_old[n])/2
                spreads_new[n] = spreads_old[n]/2
            locations_new[-1] = locations_old[-1] + spreads_old[-1]
            spreads_new[-1] = spreads_old[-1]/2

        # Strategy #3: Spread N components equidistantly within box boundaries
        elif strategy=='spread': 
            locations_new = np.squeeze(np.linspace(lb[areLocations],ub[areLocations],Ncomp+2)[1:-1])
            spreads_new = np.squeeze(np.linspace(lb[areSpreads],ub[areSpreads],Ncomp+2)[1:-1])

        # Ensure box constraints of new parameter start values
        locations_new = np.maximum(locations_new,lb[areLocations])
        locations_new = np.minimum(locations_new,ub[areLocations])
        spreads_new = np.maximum(spreads_new,lb[areSpreads])
        spreads_new = np.minimum(spreads_new,ub[areSpreads])

        # Mount the new par0 vector for the N components
        areLocations_new = np.tile(areLocations,Ncomp)
        areSpreads_new = np.tile(areSpreads,Ncomp)
        par0_new = np.tile(par[0:len(lb)],Ncomp)
        par0_new[areLocations_new] = locations_new
        par0_new[areSpreads_new] = spreads_new

        return par0_new
    #===============================================================================

    # Pre-allocate containers
    fits,Vfit,Pfit,plin_,pnonlin_,nlin_ub_,nlin_lb_,lin_ub_,lin_lb_ = ([] for _ in range(9))
    logest = []
    par_prev = None

    if strategy=='spread' or strategy=='split':
        # Start with single components, add sequentially (forward)
        Ncomponents = np.arange(1,maxModels+1)
    elif strategy=='merge':
        # Start with maximum components, remove sequentially (backward)
        Ncomponents = np.arange(maxModels,0,-1)

    # Loop over number of components in model
    # =======================================
    for Ncomp in Ncomponents:
        
        # Prepare non-linear model with N-components
        # ===========================================
        Knonlin = lambda par: nonlinmodel(par,Ncomp)
        
        # Start values of non-linear parameters
        par0_K = (ubK - lbK)/2 # start in the middle
        par0_P = initialize_nonlin_param(par_prev,Ncomp,lb,ub,strategy) # start spread within boundaries
        par0 = np.concatenate((par0_K, par0_P),axis=None)

        # Box constraints for the model parameters (non-linear parameters)
        nlin_lb = np.tile(lb,(1,Ncomp))
        nlin_ub = np.tile(ub,(1,Ncomp))

        # Add the box constraints on the non-linear kernel parameters
        nlin_lb = np.concatenate((lbK, nlin_lb),axis=None)
        nlin_ub = np.concatenate((ubK, nlin_ub),axis=None)
        
        # Box constraints for the components amplitudes (linear parameters)
        lin_lb = np.ones(Ncomp)        
        lin_ub = np.full(Ncomp,np.inf)
        
        # Separable non-linear least-squares (SNLLS) fit
        upscale = 1e2
        fit = dl.snlls(V*upscale,Knonlin,par0,nlin_lb,nlin_ub,lin_lb,lin_ub, 
                        weights=weights, reg=False, nonlin_tol=tol, nonlin_maxiter=maxiter)
        pnonlin = fit.nonlin
        plin = fit.lin
        par_prev = pnonlin

        plin /= upscale
        fit.model /= upscale
        fit.modelUncert = fit.modelUncert.propagate(lambda x: x/upscale)

        # Store the fitted parameters
        pnonlin_.append(pnonlin)
        plin_.append(plin)

        # Get fitted kernel
        Kfit = nonlinmodel(pnonlin,Ncomp)
        
        # Get fitted signal
        Vfit.append(Kfit@plin)
        
        # Get fitted distribution
        Pfit.append(Pmodel(pnonlin,plin))

        # Likelihood estimators
        # =====================
        logest = logestimators(V,Vfit[-1],plin,pnonlin,logest)

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
    Nopt = Ncomponents[idx]
    fit,Pfit,Vfit,pnonlin,plin = (var[idx] for var in [fits,Pfit,Vfit,pnonlin_,plin_])
    nlin_lb,nlin_ub,lin_lb,lin_ub = (var[idx] for var in [nlin_lb_,nlin_ub_,lin_lb_,lin_ub_])

    # Package the fitted parameters
    # =============================
    fitparam_K = pnonlin[0:nKparam] # Kernel parameters
    fitparam_P = pnonlin[nKparam:nKparam+nparam] # Components parameters
    fitparam_amp = plin # Components amplitudes

    # Uncertainty quantification analysis (if requested)
    # ==================================================
    if uq:
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
        paramuq = dl.UQResult('covariance',np.concatenate((pnonlin, plin)),covmatrix,lb_full,ub_full)
        
        P_subset = np.arange(0, nKparam+nparam*Nopt)
        amps_subset = np.arange(P_subset[-1]+1, P_subset[-1]+1+Nopt)
        Puq = paramuq.propagate(lambda p: Pmodel(p[P_subset],p[amps_subset]), np.zeros(len(r)))
    else: 
        Puq = dl.UQResult('void')
        paramuq = dl.UQResult('void')

    # If requested re-normalize the distribution
    postscale = np.trapz(Pfit,r)
    if renormalize:
        Pfit = Pfit/postscale
        if uq:
            Puq_ = copy.deepcopy(Puq) # need a copy to avoid infite recursion on next step
            Puq = Puq_.propagate(lambda P: P/postscale, lbm=np.zeros_like(Pfit))
        postscale /= sum(fitparam_amp)
        fitparam_amp = fitparam_amp/sum(fitparam_amp)

    # Dataset scales
    scales = []
    for i in range(len(Vsubsets)): 
        scales.append(prescales[i]*postscale)

    # Get fitted signals and their uncertainty
    modelfit, modelfituq = [],[]
    for i,subset in enumerate(Vsubsets): 
        V[subset] = V[subset]*prescales[i]
        modelfit.append(  scales[i]*fit.model[subset] )
        if uq: 
            modelfituq.append( fit.modelUncert.propagate(lambda V: scales[i]*V[subset]) )
        else:
            modelfituq.append(dl.UQResult('void'))

    # Goodness of fit
    stats = []
    for i,subset in enumerate(Vsubsets): 
        Ndof = len(V[subset]) - (nKparam + nparam + Nopt)
        stats.append(goodness_of_fit(V[subset],modelfit[i],Ndof))

    # Results display function
    def plotfcn(show=True):
        fig = _plot(Vsubsets,V,modelfit,modelfituq,r,Pfit,Puq,fcnals,maxModels,method,uq,show)
        return fig

    # If just one dataset, return vector instead of list
    if len(Vsubsets)==1: 
        stats = stats[0]
        scales = scales[0]
        modelfit = modelfit[0]
        modelfituq = modelfituq[0]


    return FitResult(P=Pfit, Pparam=fitparam_P, Kparam=fitparam_K, amps=fitparam_amp, V=modelfit, Puncert=Puq, 
                    paramUncert=paramuq, Vuncert=modelfituq, selfun=fcnals, Nopt=Nopt, Pn=Peval, scale=scales, plot=plotfcn,
                    stats=stats, cost=fit.cost, residuals=fit.residuals, success=fit.success)
# =========================================================================


def _plot(Vsubsets,V,Vfit,Vuq,r,Pfit,Puq,fcnals,maxModels,method,uq,show):
# =========================================================================
    if not isinstance(Vfit, list): 
        Vuq = [Vuq]
        Vfit = [Vfit]

    nSignals = len(Vsubsets)
    fig,axs = plt.subplots(nSignals+1,figsize=[7,3+3*nSignals])
    for i in range(nSignals): 
        subset = Vsubsets[i]
        # Plot the experimental signal and fit
        axs[i].plot(V[subset],'.',color='grey',alpha=0.5)
        axs[i].plot(Vfit[i],'tab:blue')
        if uq:
            # Confidence intervals of the fitted datasets
            Vci95 = Vuq[i].ci(95) # 95#-confidence interval
            Vci50 = Vuq[i].ci(50) # 50#-confidence interval
            tax = np.arange(len(subset))
            axs[i].fill_between(tax,Vci50[:,0],Vci50[:,1],color='tab:blue',linestyle='None',alpha=0.45)
            axs[i].fill_between(tax,Vci95[:,0],Vci95[:,1],color='tab:blue',linestyle='None',alpha=0.25)
        axs[i].grid(alpha=0.3)
        axs[i].set_xlabel('Array Elements')
        axs[i].set_ylabel(f'V[{i}]')
        axs[i].legend(('Data','Fit'))

    if uq:
        # Confidence intervals of the fitted distance distribution
        Pci95 = Puq.ci(95) # 95#-confidence interval
        Pci50 = Puq.ci(50) # 50#-confidence interval

    ax = plt.subplot(nSignals+1,2,2*(nSignals+1)-1)
    ax.plot(r,Pfit,color='tab:blue',linewidth=1.5)
    if uq:
        ax.fill_between(r,Pci50[:,0],Pci50[:,1],color='tab:blue',linestyle='None',alpha=0.45)
        ax.fill_between(r,Pci95[:,0],Pci95[:,1],color='tab:blue',linestyle='None',alpha=0.25)
    ax.grid(alpha=0.3)
    if uq:
        ax.legend(['truth','optimal fit','95%-CI'])
    else:
        ax.legend(['truth','optimal fit'])
    ax.set_xlabel('Distance (nm)')
    ax.set_ylabel('P (nm⁻¹)')
    axs = np.append(axs,ax)

    # Compute the Akaike weights
    dfcnals = fcnals - min(fcnals)
    ax = plt.subplot(nSignals+1,2,2*(nSignals+1))
    ax.bar(np.arange(maxModels)+1,np.log10(1 + dfcnals + abs(min(dfcnals))),facecolor='tab:blue',alpha=0.6)
    ax.grid(alpha=0.3)
    ax.set_ylabel(f'log10 Δ{method.upper()}')
    ax.set_xlabel('Number of components')
    axs = np.append(axs,ax)

    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.close()
    return fig
# =========================================================================


# fitmodel.py - Dipolar signal fit function
# ---------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
import copy
import matplotlib.pyplot as plt
import deerlab as dl
from collections.abc import Sequence
from deerlab.classes import UQResult, FitResult
from deerlab.bg_models import bg_hom3d
from deerlab.ex_models import ex_4pdeer
from deerlab.utils import isempty, goodness_of_fit, Jacobian

def fitmodel(Vexp, t, r, dd_model='P', bg_model=bg_hom3d, ex_model=ex_4pdeer,
              dd_par0=None, bg_par0=None, ex_par0=None, verbose=False, 
              dd_lb=None, bg_lb=None, ex_lb=None, dd_ub=None, bg_ub=None, ex_ub=None,
              weights=None, uq='covariance', regparam='aic', tol=1e-10,maxiter=1e8):
    r"""
    Fits a dipolar model to the experimental signal ``V`` with time axis ``t``, using
    distance axis ``r``. The model is specified by the distance distribution (dd),
    the background (bg), and the experiment (ex).

    If multiple signals (``V1``, ``V2``,...) and their corresponding time axes (``t1``, ``t2``,...)
    are given, they will be fitted globally with a single distance distribution (``dd``).
    For each signal, a specific background (``bg1``, ``bg2``,...) and experiment (``ex1``, ``ex2``)
    models can be assigned.

    This function can handle both parametric and non-parametric distance distribution models. It 
    accepts only built-in parametric models, no user-defined models.

    Parameters
    ----------
    V : array_like or list of array_like
        Time-domain signal(s) to fit.
    
    t : array_like or list of array_like
        Time axis, in microseconds. 
    
    r : array_like
        Distance axis, in nanometers.
    
    dd : callable or string
        Distance distribution model, the following modes are allowed:

        * A string ``'P'`` to indicate a non-parametric distribution
        * A callable function of a DeerLab parametric distribution model (e.g. ``dd_gauss``).
        * A string ``None`` to indicate no distribution, i.e. only background.

        The default is ``'P'``.

    bg : callable or string
        Background model, the following modes are allowed:

        * A callable function of a DeerLab parametric background model (e.g. ``bg_hom3d``).
        * ``None`` to indicate no background decay.

        The default is ``bg_hom3d``.

    ex : callable or string
        Experiment model, the following modes are allowed:

        * Function handle to experiment model (e.g. ``ex_4pdeer``)
        * ``None`` to indicate simple dipolar oscillation (modulation depth equal 1)

        The default is ``ex_4pdeer``.

    dd_par0, bg_par0, ex_par0 : array_like, optional
        Initial parameter values of the distance distribution/background/experimental model parameters.
        If a model does not require parameters or are to be determined automatically it can be omitted or specified 
        as ``None`` (default).
    
    dd_lb, bg_lb, ex_lb : array_like, optional    
        Lower boundary values of the distance distribution/background/experimental model parameters.
        If a model does not require parameters or are to be determined automatically it can be omitted or specified 
        as ``None`` (default).
    
    dd_ub, bg_ub, ex_ub : array_like, optional    
        Upper boundary values of the distance distribution/background/experimental model parameters.
        If a model does not require parameters or are to be determined automatically it can be omitted or specified 
        as ``None`` (default).
    
    uq : string or list, optional
        Type of uncertainty quantification analysis. Any ``UQResult`` output returned by this function will
        be adjusted accordingly. The options are:

        * ``'covariance'`` - Covariance-based uncertainty quantification. Fast, but approximate.   
        * ``'bootstrap'`` - Bootstrapped uncertainty quantification. Slow, but accurate. By default, 1000 bootstrap
        samples are used. Alternatively, a different number can be specified as follows ``uq=['bootstrap',Nsamples]``.
        * ``None`` - Disable the uncertainty quantification analysis. 

        The default is ``'covariance'``.

    verbose : boolean, optional
        Enable/disable printing a table of fit results, by default is disabled

    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting,
        the default is all weighted equally. If not specified all datasets are weighted inversely proportional to their noise levels.
    
    regparam : str or scalar, optional
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

    tol : scalar, optional 
        Tolerance value for convergence of the NNLS algorithm. If not specified, the value is set to ``tol = 1e-10``.
        
    maxiter: scalar, optional  
        Maximum number of iterations before termination. If not specified, the value is set to ``maxiter = 1e8``.

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    V :  ndarray or list thereof
        Fitted time-domain signal(s).
    P : ndarray
        Fitted distance distribution.
    B : ndarray or list thereof
        Fitted background decay(s).
    Vmod : ndarray or list thereof
        Fitted modulated contribution (foreground and background).
    Vunmod : ndarray or list thereof
        Fitted unmodulated contribution (foreground and background). 
    ddparam : ndarray
        Fitted parameters for distance distribution model.
    bgparam : ndarray or list thereof
        Fitted parameters for background model(s).
    exparam : ndarray or list thereof
        Fitted parameters for experiment model(s).
    Vuncert : :ref:`UQResult` or list thereof
        Uncertainty quanfitication for fitted dipolar signal(s).
    Puncert : :ref:`UQResult`
        Uncertainty quanfitication for fitted distance distribution.
    Buncert : :ref:`UQResult` or list thereof
        Uncertainty quanfitication for fitted background(s).
    VmodUncert : ndarray or list thereof
        Uncertainty quanfitication for fitted modulated contribution(s). 
    VunmodUncert : ndarray or list thereof
        Uncertainty quanfitication for fitted unmodulated contribution(s). 
    ddparamUncert : :ref:`UQResult` 
        Uncertainty quanfitication for distribution parameters
    bgparamUncert : :ref:`UQResult` or list thereof
        Uncertainty quanfitication for background(s) parameters
    exparamUncert : :ref:`UQResult` or list thereof
        Uncertainty quanfitication for experiment(s) parameters
    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    regparam : scalar
        Optimal regularization parameter value
    plot : callable
        Function to display the results. It will 
        display the fitted signals and distance distributions with
        confidence intervals. The function returns 
        the figure object (``matplotlib.figure.Figure``) object as output, 
        which can be modified. Using ``fig = plot(show=False)`` will not render
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
    Fit a 4pDEER signal with homogenous 3D background with Gaussian distribution::

        fit = dl.fitmodel(V,t,r,dl.dd_gauss,dl.bg_hom3d,dl.ex_4pdeer)  


    Fit a 5pDEER signal with exponential background and Tikhonov regularization::

        fit = dl.fitmodel(V,t,r,'P',dl.bg_hom3d,dl.ex_5pdeer)
    
    
    Fit a 4pDEER stretched exponential background (no foreground)::
    
        fit = dl.fitmodel(V,t,r,None,dl.bg_strexp,None)  
    
    
    Fit a dipolar evolution function with Rician distribution::
    
        fit = dl.fitmodel(V,t,r,dl.dd_rice,None,None)
    
    
    Fit a 4pDEER form factor (no background) with bimodal Gaussian distribution:: 

        fit = dl.fitmodel(V,t,r,dl.dd_gauss2,None,dl.ex_4pdeer)   
   
     
    Fit a 4pDEER signal and a 5pDEER signal with same type of background::
    
        fit = fitmodel([V1,V2],t,r,'P',dl.bg_hom3d,[dl.ex_4pdeer,dl.ex_5pdeer])    
    """

    # Default optional settings
    normP = True
    bootsamples = 1000

    # Make inputs into a list if just a singal is passed
    Vexp,t,bg_model,ex_model = ([var] if not isinstance(var, Sequence) else var for var in (Vexp,t,bg_model,ex_model))

    # Basic validation of experimental signals
    nSignals = len(Vexp)
    if len(t)!=nSignals:
        raise KeyError('The same number of signals V and time axes t must be provided.')
    for i in range(nSignals):
        if len(Vexp[i])!=len(t[i]):
            raise ValueError(f'V[{i}] and t[{i}] must have the same number of elements.')
        if not np.all(np.isreal(Vexp[i])):
            raise ValueError('Input signals cannot be complex-valued.')
    
    # If just one model is specified, assume same model for all signals
    if len(bg_model)!=nSignals:
        bg_model = bg_model*nSignals
    if len(ex_model)!=nSignals:
        ex_model = ex_model*nSignals

    # Process uncertainty quantification settings
    if isinstance(uq, str) or uq is None:
        uq = [uq]
    if uq[0]!='bootstrap' and uq[0]!='covariance'and uq[0] is not None:
        raise KeyError("Uncertainty quantification must be either 'covariance', 'bootstrap', or None.")
    if uq[0]=='bootstrap':
        # Override default if user has specified bootstraped samples
        if len(uq)>1:
            bootsamples = uq[1]
    
    uq = uq[0]
    uqanalysis = uq is not None

    # Combine input boundary and start conditions
    par0 = [[] if par0_i is None else par0_i for par0_i in [dd_par0,bg_par0,ex_par0]]
    lb = [[] if lb_i is None else lb_i for lb_i in [dd_lb,bg_lb,ex_lb]]
    ub = [[] if ub_i is None else ub_i for ub_i in [dd_ub,bg_ub,ex_ub]]
     
    # Get information about distance distribution parameters
    # ------------------------------------------------------
    parametricDistribution  = type(dd_model) is types.FunctionType # Check if P(r) is parametric model
    parfreeDistribution = dd_model == 'P' # Check if P(r) is non-parametric model
    includeForeground = np.array(dd_model is not None) # Check if there is a foreground or just background in the model

    par0_dd, lower_dd, upper_dd, N_dd = _getmodelparams(dd_model) # Get model built-in info

    # Catch nonsensical case
    if includeForeground and not parametricDistribution and not parfreeDistribution:
        raise TypeError('Distribution model (4th input) must either be a function handle, ''P'', or None.')

    # Get information about background parameters
    # ------------------------------------------------------
    isparamodel = np.array([type(model) is types.FunctionType for model in bg_model]) # Check if B(t) is a parametric model
    includeBackground = np.array([model is not None for model in bg_model]) # Check if model includes B(t) or not

    par0_bg, lower_bg, upper_bg, N_bg = _getmodelparams(bg_model) # Get model built-in info

    # Check whether the background model is a physical or phenomenological model
    isphenomenological = [_iserror(model,t,par0,1) for model,par0 in zip(bg_model,par0_bg)]

    # Catch nonsensical case
    if np.any(~isparamodel & includeBackground):
        raise TypeError('Background model (5th input) must either be a function handle, or None.')
    
    # Get information about experiment parameters
    # ------------------------------------------------------
    isparamodel = np.array([type(model) is types.FunctionType for model in ex_model]) # Check if experiment is a parametric model
    includeExperiment = np.array([model is not None for model in ex_model]) # Check whether to include experiment in the model

    par0_ex, lower_ex, upper_ex, N_ex = _getmodelparams(ex_model) # Get model built-in info

    # Catch nonsensical situations
    if np.any(~isparamodel & includeExperiment):
        raise TypeError('Experiment models must either be a function handle, or None.')

    # Prepare full-parameter set
    # -------------------------------------------------------
    # Combine all initial, lower and upper bounds of parameters into vectors
    par0,par0_defined = _parcombine(par0,[par0_dd,par0_bg,par0_ex])
    lb,_ = _parcombine(lb,[lower_dd,lower_bg,lower_ex])
    ub,_ = _parcombine(ub,[upper_dd,upper_bg,upper_ex])
    lb,ub,par0 = _checkbounds(lb,ub,par0,par0_defined)

    # Build index vectors for accessing parameter subsets
    ddidx = np.arange(N_dd[0])
    bgidx,exidx = ([],[])
    for ii in range(nSignals):
        bgidx.append(N_dd + sum(N_bg[np.arange(0,ii)]) + np.arange(0,N_bg[ii]) )
        exidx.append(N_dd + sum(N_bg) + sum(N_ex[np.arange(0,ii)]) + np.arange(0,N_ex[ii]) )


    def multiPathwayModel(par):
    # =========================================================================
        """
        Multi-pathway dipolar model
        ------------------------------------------------------------------
        This function represent the core model of fitmodel, it takes the
        parameters and computes the dipolar kernel and background according to
        the dipolar multi-pathway theory. 
        This is the non-linear part of the SNLLS problem when fitting a
        parameter-free distribution.
        """
        Bs,Ks =([],[])
        for iSignal in range(nSignals):            
            bg_par = par[bgidx[iSignal]]
            ex_par = par[exidx[iSignal]]

            # Prepare background basis function
            if includeBackground[iSignal]:
                if isphenomenological[iSignal]:
                    Bfcn = lambda t: bg_model[iSignal](t,bg_par)
                else:
                    Bfcn = lambda t,lam: bg_model[iSignal](t,bg_par,lam)
            else:
                Bfcn = lambda _,__: np.ones_like(Vexp[iSignal])
            
            # Get pathway information
            if includeExperiment[iSignal]:
                pathways = ex_model[iSignal](ex_par)
            else:
                pathways = [[1,0]]
            
            # Compute the multipathway-background
            B_ = dl.dipolarbackground(t[iSignal],pathways,Bfcn)
            # Compute the multipathway-kernel
            K_ = dl.dipolarkernel(t[iSignal],r,pathways=pathways,bg=Bfcn)
            Ks.append(K_)
            Bs.append(B_)

        return Ks, Bs     
    # =========================================================================

    def splituq(param_uq,Pfit_uq,Pfit,Vfit,Bfit,parfit_,Kfit,scales=1):
    # =========================================================================
        """ 
        Uncertainty quantification
        ---------------------------
        This function takes the full combined covariance matrix of the
        linear+nonlinear parameter sets and splits it into the different
        components of the model. 
        """
        # Pre-allocation
        paruq_bg,paruq_ex,Bfit_uq,Vmod_uq,Vunmod_uq,Vfit_uq = ([],[],[],[],[],[])

        # Background parameters uncertainty
        # ---------------------------------
        for jj in range(nSignals):
            if includeBackground[jj]:
                bgsubcovmat  = param_uq.covmat[np.ix_(bgidx[jj],bgidx[jj])]
                paruq_bg.append( UQResult('covariance',parfit_[bgidx[jj]],bgsubcovmat,lb[bgidx[jj]],ub[bgidx[jj]]))
            else:
                paruq_bg.append([None])
        
        # Experiment parameters uncertainty
        # ----------------------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                exsubcovmat  = param_uq.covmat[np.ix_(exidx[jj],exidx[jj])]
                paruq_ex.append( UQResult('covariance',parfit_[exidx[jj]],exsubcovmat,lb[exidx[jj]],ub[exidx[jj]]))
            else:
                paruq_ex.append([None])
            
        # Distribution parameters uncertainty
        # ------------------------------------
        if parametricDistribution:
            ddsubcovmat  = param_uq.covmat[np.ix_(ddidx,ddidx)]
            paruq_dd = UQResult('covariance',parfit_[ddidx],ddsubcovmat,lb[ddidx],ub[ddidx])
        else:
            paruq_dd = [None]
        
        # Distance distribution uncertainty
        # ----------------------------------
        nonneg = np.zeros_like(r)
        if parametricDistribution:
            # Prepare parametric model
            if includeForeground:
                Pfcn = lambda par: dd_model(r,par[ddidx])
            else:
                Pfcn = lambda _: np.ones_like(r)/np.trapz(np.ones_like(r),r)
            Pfit_uq = param_uq.propagate(Pfcn,nonneg)
        
        # Background uncertainty
        # -----------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                Bfit_uq.append( param_uq.propagate(lambda par:scales[jj]*multiPathwayModel(par)[1][jj]) )
            else:
                Bfit_uq.append([None])
        
        # Unmodulated contribution uncertainty
        # -----------------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                Lam0fcn = lambda par: ex_model[jj](par)[0]
                Bfcn = lambda par: scales[jj]*multiPathwayModel(par)[1][jj]
                Vunmod_fcn = lambda par: Lam0fcn(par[exidx[jj]])*Bfcn(par)
                Vunmod_uq.append( param_uq.propagate(lambda par:Vunmod_fcn(par)) )
            else:
                Vunmod_uq.append([None])
        
        # Dipolar signal uncertainty
        # --------------------------
        for jj in range(nSignals):
            if includeForeground and parametricDistribution:
                # Full parametric signal
                Vmodel = lambda par: scales[jj]*multiPathwayModel(par)[0][jj]@Pfcn(par[ddidx])
                Vfit_uq.append( param_uq.propagate(Vmodel))
            elif includeForeground and np.all(~includeExperiment & ~includeBackground):
                Vmodel = lambda _: Kfit[jj]@Pfit
                # Dipola evolution function
                J = Kfit[jj]
                Vcovmat = J@Pfit_uq.covmat@J.T
                Vfit_uq.append( UQResult('covariance',Vfit[jj],Vcovmat))
            elif includeForeground:
                # Parametric signal with parameter-free distribution
                Vmodel = lambda par: scales[jj]*multiPathwayModel(par)[0][jj]@Pfit
                Vfit_uq.append( param_uq.propagate(Vmodel) )
            else:
                Vfit_uq.append([None])

        # Modulated contribution uncertainty
        # -----------------------------
        for jj in range(nSignals):
            if includeForeground and np.all(~includeExperiment & ~includeBackground):
                Vmod_uq.append(Vfit_uq) 
            elif includeForeground:
                Vmod_fcn = lambda par: Vmodel(par) - Vunmod_fcn(par)
                Vmod_uq.append( param_uq.propagate(Vmod_fcn))
            else: 
                Vmod_uq.append([None]) 

        return Vfit_uq, Pfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd   
    # =========================================================================


    def calculate_Vmod_Vunmod(parfit,Vfit,Bfit,scales):
    # =========================================================================
        " Calculation of the (un)modulated components of the dipolar signal" 

        # Calculate the unmodulated contribution (Vunmod)
        # --------------------------------------------------------
        Vunmod = []
        for j in range(nSignals):
            if includeExperiment[j]:
                Lam0 = ex_model[j](parfit[exidx[j]])[0][0]
                if includeBackground[j]:
                    Vunmod.append(Lam0*np.array(Bfit[j]))
                else:
                    Vunmod.append(np.full_like(t[j],scales[j]*Lam0))
            else:
                Vunmod.append(np.zeros_like(t[j]))

        # Calculate the modulated contribution (Vmod)
        # --------------------------------------------------------
        Vmod = []
        for j in range(nSignals):
            Vmod.append(Vfit[i] - Vunmod[i])

        return Vmod, Vunmod
    # =========================================================================

    def regularization_analysis(Vexp):
    # =========================================================================
        " Analysis workflow for non-parametric models based on regularized least-squares" 
        
        # Use basic dipolar kernel
        Ks = [dl.dipolarkernel(ts,r) for ts in t]
        
        # Linear regularization fit
        fit = dl.fitregmodel(Vexp,Ks,r,'tikhonov',regparam, weights=weights,uq=uqanalysis,tol=tol,maxiter=maxiter)
        Pfit = fit.P
        Pfit_uq = fit.Puncert
        scales = np.atleast_1d(fit.scale)
        alphaopt = fit.regparam

        # Get fitted models
        Vfit = [scale*K@Pfit for K,scale in zip(Ks,scales)]
        Bfit = [scale*np.ones_like(V) for V,scale in zip(Vexp,scales)]
        Vmod, Vunmod = calculate_Vmod_Vunmod(None,Vfit,Bfit,scales)

        # No parameters
        parfit = np.asarray([None])

        if uqanalysis and uq=='covariance':
            Vfit_uq, Pfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd = splituq(None,Pfit_uq,Pfit,Vfit,Bfit,parfit,Ks,scales)
            return fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit, Pfit_uq, Vfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd, scales, alphaopt
        else:
            return fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit, scales, alphaopt
    # =========================================================================

    def nonlinear_lsq_analysis(Vexp):
    # =========================================================================
        " Analysis workflow for fully parametric models based on nonlinear least-squares" 

        # Prepare the full-parametric model
        if includeForeground:
            Pfcn = lambda par: dd_model(r,par[ddidx])
        else:
            Pfcn = lambda _: np.ones_like(r)/np.trapz(np.ones_like(r),r)
        Vmodel =lambda par: [K@Pfcn(par) for K in multiPathwayModel(par)[0]]

        # Non-linear parametric fit
        fit = dl.fitparamodel(Vexp,Vmodel,par0,lb,ub,weights=weights,uq=uqanalysis,tol=tol,maxiter=maxiter)
        parfit = fit.param
        param_uq = fit.paramUncert
        scales = np.atleast_1d(fit.scale)
        alphaopt = None

        # Get fitted models
        Vfit = Vmodel(parfit)
        _,Bfit = multiPathwayModel(parfit)
        if includeForeground:
            Pfit = Pfcn(parfit)
        else:
            Pfit = []
        if type(Vfit) is not list:
            Vfit = [Vfit]
        if type(Bfit) is not list:
            Bfit = [Bfit]
        Bfit = [scale*B for B,scale in zip(Bfit,scales)]
        Vfit = [scale*V for V,scale in zip(Vfit,scales) ]
        Vmod, Vunmod = calculate_Vmod_Vunmod(parfit,Vfit,Bfit,scales)

        if uqanalysis and uq=='covariance':
            Vfit_uq, Pfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd = splituq(param_uq,None,Pfit,Vfit,Bfit,parfit,None, scales)
            return fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit, Pfit_uq, Vfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd,scales,alphaopt
        else:
            return fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit, scales, alphaopt
    # =========================================================================

    
    def separable_nonlinear_lsq_analysis(Vexp):
    # =========================================================================
        " Analysis workflow for semiparametric models based on separable nonlinear least-squares" 

        # Non-negativity constraint on distributions
        lbl = np.zeros_like(r)
        
        #prescales = [max(V) for V in Vexp]
        prescales = [1 for V in Vexp]
        Vexp_ = [Vexp[i]/prescales[i] for i in range(nSignals)]

        def scale_constraint(nonlinpar,linpar):
        # --------------------------------------------------------
            penalty = np.zeros(nSignals)
            for i in range(nSignals):
                ex_par = nonlinpar[exidx[i]]
                pathways = ex_model[i](ex_par)
                lams = [pathway[0] for pathway in pathways]
                if np.sum(lams)<1 or ex_model[i].__name__=='ex_4pdeer': 
                    penalty[i] = 0
                else:
                    penalty[i] = max(Vexp[i])*(np.sum(lams) - 1)

            return penalty 
        # --------------------------------------------------------


        # Separable non-linear least squares (SNNLS) 
        fit = dl.snlls(Vexp_,lambda par: multiPathwayModel(par)[0],par0,lb,ub,lbl, reg=True,
                            regparam=regparam, uq=uqanalysis, weights=weights,extrapenalty=scale_constraint,
                            nonlin_tol=tol,nonlin_maxiter=maxiter)
        parfit = fit.nonlin
        Pfit = fit.lin
        param_uq = fit.nonlinUncert
        Pfit_uq = fit.linUncert
        alphaopt = fit.regparam
        scales = fit.scale

        # Normalize distribution
        # -----------------------
        Pscale = np.trapz(Pfit,r)
        Pfit /= Pscale
        scales = np.atleast_1d([scale*Pscale for scale in scales])
        if uqanalysis and uq=='covariance':
            # scale CIs accordingly
            Pfit_uq_ = copy.deepcopy(Pfit_uq) # need a copy to avoid infite recursion on next step
            Pfit_uq = Pfit_uq_.propagate(lambda P: P/Pscale,lbm=lbl)

        # Get the fitted models
        Kfit,Bfit = multiPathwayModel(parfit)
        Bfit = [scale*B for B,scale in zip(Bfit,scales)]
        Vfit = [scale*K@Pfit for K,scale in zip(Kfit,scales)]
        Vmod, Vunmod = calculate_Vmod_Vunmod(parfit,Vfit,Bfit,scales)

        if uqanalysis and uq=='covariance':
            Vfit_uq, _, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd = splituq(param_uq,Pfit_uq, Pfit, Vfit, Bfit, parfit, Kfit, scales)
            return fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit, Pfit_uq, Vfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd,scales,alphaopt
        else:
            return fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit, scales, alphaopt
    # =========================================================================

    # Analyze the data
    # ----------------------

    # Determine type of model
    nonparametric = np.all(~parametricDistribution & ~includeExperiment & ~includeBackground)
    fullparametric = not nonparametric and (parametricDistribution or not includeForeground)
    semiparametric = not nonparametric and not fullparametric

    # Choose appropiate analysis for type of model
    if nonparametric:
        analysis = regularization_analysis
    elif fullparametric:
        analysis = nonlinear_lsq_analysis
    elif semiparametric:
        analysis = separable_nonlinear_lsq_analysis

    # Run the analysis
    results = analysis(Vexp)

    # Unpack results
    if uqanalysis and uq=='covariance':
        fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit_, Pfit_uq, Vfit_uq, Bfit_uq, Vmod_uq, Vunmod_uq, paruq_bg, paruq_ex, paruq_dd, scales, alphaopt = results
    else:
        fit, Pfit, Vfit, Bfit, Vmod, Vunmod, parfit_, scales, alphaopt = results

    # Bootstrapping uncertainty quantification
    # -----------------------------------------
    if uqanalysis and uq=='bootstrap':

        def bootstrapfcn(Vexp):
            # ======================================================
            # Fit the data
            _, Pfit_, Vfit_, Bfit_, Vmod_, Vunmod_, parfit, _, _ = analysis(Vexp)
            # Extract the individual parameter subsets
            parfit_bg = [parfit[bgidx[n]] for n in range(nSignals)]
            parfit_ex = [parfit[exidx[n]] for n in range(nSignals)]
            parfit_dd = parfit[ddidx]

            return (Pfit_,*Vfit_,*Bfit_,*Vmod_,*Vunmod_,*parfit_bg,*parfit_ex,parfit_dd)
            # ======================================================

        # Run bootstrapping
        boot_uq = dl.bootan(bootstrapfcn,Vexp,Vfit,samples=bootsamples,verbose=verbose)

        # Unpack bootstrapping results
        Pfit_uq   =  boot_uq[0]
        Vfit_uq   = [boot_uq[1+n] for n in range(nSignals)]
        Bfit_uq   = [boot_uq[1+nSignals+n] for n in range(nSignals)]
        Vmod_uq   = [boot_uq[1+2*nSignals+n] for n in range(nSignals)]
        Vunmod_uq = [boot_uq[1+3*nSignals+n] for n in range(nSignals)]
        paruq_bg  = [boot_uq[1+4*nSignals+n] for n in range(nSignals)]
        paruq_ex  = [boot_uq[1+5*nSignals+n] for n in range(nSignals)]
        paruq_dd  =  boot_uq[-1]

    # Normalize distribution
    # -----------------------
    scale = 1
    if normP and includeForeground:
        scale = np.trapz(Pfit,r)
        Pfit /= scale
        if uqanalysis:
            # scale CIs accordingly
            Pfit_uq_ = copy.deepcopy(Pfit_uq) # need a copy to avoid infite recursion on next step
            Pfit_uq = Pfit_uq_.propagate(lambda P: P/scale,lbm=np.zeros_like(r))


    # Calculate goodness of fit
    # -------------------------
    stats = []
    for j in range(nSignals):
        Ndof = len(Vexp[j]) - len(parfit_)
        stats.append(goodness_of_fit(Vexp[j],Vfit[j],Ndof))


    # Return fitted parameters and confidence intervals in structures
    # ---------------------------------------------------------------
    parfit = dict()
    parfit['dd'] = parfit_[ddidx]
    parfit['bg'] = [parfit_[bgidx[j]] for j in range(nSignals)]
    parfit['ex'] = [parfit_[exidx[j]] for j in range(nSignals)]
    if uqanalysis:
        paruq = dict()
        paruq['dd'] = paruq_dd
        paruq['bg'] = [paruq_bg[j] for j in range(nSignals)]
        paruq['ex'] = [paruq_ex[j] for j in range(nSignals)]
        modfituq = dict()
        modfituq['Pfit'] = Pfit_uq
        modfituq['Vfit'] = [Vfit_uq[j] for j in range(nSignals)]
        modfituq['Vmod'] = [Vmod_uq[j] for j in range(nSignals)]
        modfituq['Vunmod'] = [Vunmod_uq[j] for j in range(nSignals)]
        modfituq['Bfit'] = [Bfit_uq[j] for j in range(nSignals)]
    else:
        paruq = dict()
        paruq['dd'] = UQResult('void')
        paruq['bg'] = UQResult('void')
        paruq['ex'] = UQResult('void')
        modfituq = dict()
        modfituq['Pfit'] = UQResult('void')
        modfituq['Vfit'] = UQResult('void')
        modfituq['Bfit'] = UQResult('void')
        modfituq['Vmod'] = UQResult('void')
        modfituq['Vunmod'] = UQResult('void')


    Vfit_ = Vfit.copy()
    Vunmod_ = Vunmod.copy()
    def _display_results(show=True):
    # =========================================================================
        fig,axs = plt.subplots(nSignals+1,figsize=[7,3+3*nSignals])
        for i in range(nSignals):
            # Plot the signal
            axs[i].plot(t[i],Vexp[i],'.',color='grey',alpha=0.5)
            axs[i].plot(t[i],Vfit_[i],'tab:blue')
            if includeBackground[i]:
                axs[i].plot(t[i],Vunmod_[i],'--',color='tab:orange')
            if uqanalysis:
                # Get confidence intervals for the signal
                Vci95 = Vfit_uq[i].ci(95)
                Vci50 = Vfit_uq[i].ci(50)
                axs[i].fill_between(t[i],Vci95[:,0], Vci95[:,1],facecolor='tab:blue',linestyle='None',alpha=0.2)
                axs[i].fill_between(t[i],Vci50[:,0], Vci50[:,1],facecolor='tab:blue',linestyle='None',alpha=0.4)
                if includeBackground[i]:
                    Vunmod95 = Vunmod_uq[i].ci(95)
                    Vunmod50 = Vunmod_uq[i].ci(50)
                    axs[i].fill_between(t[i],Vunmod95[:,0], Vunmod95[:,1],facecolor='tab:orange',linestyle='None',alpha=0.2)
                    axs[i].fill_between(t[i],Vunmod50[:,0], Vunmod50[:,1],facecolor='tab:orange',linestyle='None',alpha=0.4)
            axs[i].grid(alpha=0.3)
            axs[i].set_xlabel('Time (µs)')
            axs[i].set_ylabel(f'V[{i}]')
            axs[i].legend(('Data','Vfit','Bfit','95%-CI','50%-CI'))

        # Plot the distribution
        if uqanalysis and uq=='bootstrap':
            axs[nSignals].plot(r,Pfit_uq.median,'tab:blue',label='Median')
        else:
            axs[nSignals].plot(r,Pfit,'tab:blue',label='Fit')
        if uqanalysis:
            # Get confidence intervals for the distance distribution
            Pci95 = Pfit_uq.ci(95)
            Pci50 = Pfit_uq.ci(50)
            axs[nSignals].fill_between(r,Pci95[:,0], Pci95[:,1],facecolor='tab:blue',linestyle='None',alpha=0.2,label='95%-CI')
            axs[nSignals].fill_between(r,Pci50[:,0], Pci50[:,1],facecolor='tab:blue',linestyle='None',alpha=0.4,label='50%-CI')
        axs[nSignals].set_xlabel('Distance (nm)')
        axs[nSignals].set_ylabel('P (nm⁻¹)')
        axs[nSignals].legend()
        axs[nSignals].grid(alpha=0.3)
        plt.tight_layout()
        if show:
            plt.show()
        else:
            plt.close()
        return fig
    # =========================================================================


    if verbose:
        print('----------------------------------------------------------------------------')
        print('Goodness of fit')
        for i in range(nSignals):
            print(f'  Vexp[{i}]: 𝛘2 = {stats[i]["chi2red"]:4f}  RMSD  = {stats[i]["rmsd"]:4e}')
        print('----------------------------------------------------------------------------')
        print('Fitted parameters and 95%-confidence intervals')
        if len(parfit['dd'])>0:
            print()
            if uqanalysis:
                ci = paruq['dd'].ci(95)
            else:
                ci = np.full((len(parfit['dd']),2),np.nan)
            for j in range(len(parfit['dd'])):
                c = parfit['dd'][j]
                print(f'  ddparam[{j:d}]:   {c:5.7f}  ({ci[j,0]:.7f}, {ci[j,1]:.7f})  {dd_model.parameters[j]} ({dd_model.units[j]})')

        for i in range(nSignals):
            print(f'Vfit[{i}]:')
            if scales[i]>1e2:
                print(f'  V0:  {scales[i]:.3e}  Signal scale (arb.u.)')
            else:
                print(f'  V0:  {scales[i]:2.2f}  Signal scale (arb.u.)')
            if includeBackground[i]:
                if uqanalysis:
                    ci = paruq['bg'][i].ci(95)
                else:
                    ci = np.full((len(parfit['bg'][i]),2),np.nan)
                for j in range(len(parfit['bg'][i])):
                    c = parfit['bg'][i][j]
                    print(f'  bgparam[{j:d}]:   {c:5.7f}  ({ci[j,0]:.7f}, {ci[j,1]:.7f})  {bg_model[i].parameters[j]} ({bg_model[i].units[j]})')
            if includeExperiment[i]:
                if uqanalysis:
                    ci = paruq['ex'][i].ci(95)
                else:
                    ci = np.full((len(parfit['ex'][i]),2),np.nan)
                for j in range(len(parfit['ex'][i])):
                    c = parfit['ex'][i][j]
                    print(f'  exparam[{j:d}]:   {c:5.7f}  ({ci[j,0]:.7f}, {ci[j,1]:.7f})  {ex_model[i].parameters[j]} ({ex_model[i].units[j]})')
        print('----------------------------------------------------------------------------')


    # Return numeric arrays and not lists if there is only one signal
    if nSignals==1:
        Vfit,Bfit,Vmod,Vunmod = (var[0] if type(var) is list else var for var in [Vfit,Bfit,Vmod,Vunmod])
        for subset in ('ex','bg'):
            parfit[subset] = parfit[subset][0]
            if uqanalysis:
                paruq[subset] = paruq[subset][0]
        for subset in ('Vfit','Bfit','Vmod','Vunmod'):
            if uqanalysis:
                modfituq[subset] = modfituq[subset][0]
        if not isempty(stats):
            stats = stats[0]

    # Do not return array for a single scale
    if len(scales)==1:
        scales = scales[0]


    return FitResult(V=Vfit, P=Pfit, B=Bfit, Vmod=Vmod, Vunmod=Vunmod, exparam=parfit['ex'], bgparam=parfit['bg'],
                      ddparam=parfit['dd'], Vuncert=modfituq['Vfit'], Puncert=modfituq['Pfit'], 
                      VmodUncert=modfituq['Vmod'], VunmodUncert=modfituq['Vunmod'],
                      Buncert=modfituq['Bfit'], exparamUncert=paruq['ex'], bgparamUncert=paruq['bg'],
                      ddparamUncert=paruq['dd'], regparam=alphaopt, plot=_display_results, scale=scales, stats=stats, cost=fit.cost,
                      residuals=fit.residuals, success=fit.success)
# ==============================================================================
# ==============================================================================
# ==============================================================================

            
def _getmodelparams(models):
# ==============================================================================
    "Extract parameter info from parametric models"
    if type(models) is not list:
        models = [models]

    par0,lo,up,N = ([],[],[],[])
    for model in models:
        if type(model) is types.FunctionType:
            par0.append(model.start)
            lo.append(model.lower)
            up.append(model.upper)
        else:
            par0.append([])
            lo.append([])
            up.append([])
        N.append(len(par0[-1]))

    par0 = np.concatenate(par0)
    lo = np.concatenate(lo)
    up = np.concatenate(up)
    N = np.atleast_1d(N)
    return par0,lo,up,N
# ==============================================================================

def _parcombine(p,d):
# ==============================================================================
    """
    Combine parameter subsets
    p: User-defined parameter values
    d: Built-in parameter values
    """
    defined = p.copy()
    # Loop over all three parameter subsets (dd,bg,ex)
    for k in range(3):
        p[k] = np.atleast_1d(p[k])
        d[k] = np.atleast_1d(d[k])
        # If user has not defined the values, use built-in defaults
        if isempty(p[k]):
            p[k] = d[k].flatten()
            defined[k] = np.full_like(d[k],False)
        else:
            defined[k] = np.full_like(d[k],True)
        # Otherwise use the user-defined values
        if type(p[k]) is list:
            p[k] = p[k].flatten()
            
    # Concatenate all parameter subsets into one vector
    pvec = np.concatenate(p,axis=0)
    defined = np.concatenate(defined,axis=0).astype(bool)

    return pvec,defined
# ==============================================================================

def _checkbounds(lb,ub,par0,par0_passed):
# ==============================================================================
    "Boundary checking"
    nParams = len(par0)
    if len(lb)!=nParams:
        raise ValueError('Lower bounds and initial values of parameters must have the same number of elements.')
    if len(ub)!=nParams:
        raise ValueError('Upper bounds and initial values of parameters must have the same number of elements.')
    if np.any(lb>ub):
        raise ValueError('Lower bounds cannot be larger than upper bounds.')
    # If user has specified some par0, check that it is within bounds
    if np.any(lb[par0_passed]>par0[par0_passed]) or np.any(par0[par0_passed]>ub[par0_passed]):
        raise ValueError('The specified parameters inital values must lie between lower and upper bounds.')
    # For the rest ensure that start values are withing bounds (even for user-defined bounds)
    if np.any(lb[~par0_passed]>par0[~par0_passed]) or np.any(par0[~par0_passed]>ub[~par0_passed]):
        par0[~par0_passed] = (lb[~par0_passed] + ub[~par0_passed])/2

    return lb,ub,par0
# ==============================================================================

def _iserror(func, *args, **kw):
# ==============================================================================
    try:
        func(*args, **kw)
        return False
    except Exception:
        return True
# ==============================================================================

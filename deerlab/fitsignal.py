# fitsignal.py - Dipolar signal fit function
# ---------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
import copy
import inspect
import matplotlib.pyplot as plt
import deerlab as dl
from deerlab.classes import UncertQuant, FitResult
from deerlab.bg_models import bg_hom3d
from deerlab.ex_models import ex_4pdeer
from deerlab.utils import isempty, goodness_of_fit, Jacobian

def fitsignal(Vexp, t, r, dd_model='P', bg_model=bg_hom3d, ex_model=ex_4pdeer,
              par0=[None,None,None], lb=[None,None,None], ub=[None,None,None], verbose= False,
              weights=1, uqanalysis=True, regparam='aic', regtype = 'tikhonov'):
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

    par0 : list of array_like, optional
        Starting parameter values. Must be a 3-element list ``[par0_dd,par0_bg,par0_ex]``
        containing the start values of the distribution, background and experiment models, in that order.
        If a model does not require parameters or are to be determined automatically it must be specified 
        as an empty list ``[]``. The default is ``par0=[[],[],[]]``.
    
    lb : list of array_like, optional
        Lower bounds for parameters.  Must be a 3-element list ``[lb_dd,lb_bg,lb_ex]``
        containing the start values of the distribution, background and experiment models, in that order.
        If a model does not require parameters or are to be determined automatically it must be specified 
        as an empty list ``[]``. The default is ``lb=[[],[],[]]``.
    
    ub : list of array_like, optional
        Upper bounds for parameters, Must be a 3-element list ``[ub_dd,ub_bg,ub_ex]``
        containing the start values of the distribution, background and experiment models, in that order.
        If a model does not require parameters or are to be determined automatically it must be specified 
        as an empty list ``[]``. The default is ``ub=[[],[],[]]``.
    
    weights : array_like, optional
        Array of weighting coefficients for the individual signals in global fitting,
        the default is all weighted equally.
        If not specified all datasets are weighted equally.
    
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

    regtype : string, optional
        Regularization functional type: 
    
        * ``'tikhonov'`` - Tikhonov regularizaton
        * ``'tv'``  - Total variation regularization
        * ``'huber'`` - Huber regularization
        The default is ``'tikhonov'``.  
    
    verbose : boolean, optional
        Enable/disable printing a table of fit results, by default is disabled
    
    uqanalysis : boolean, optional
        Enable/disable the uncertainty quantification analysis, by default it is enabled.



    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    V :  ndarray or list thereof
        Fitted time-domain signal(s).
    P : ndarray
        Fitted distance distribution.
    B : ndarray or list thereof
        Fitted background decay(s).
    ddparam : ndarray
        Fitted parameters for distance distribution model.
    bgparam : ndarray or list thereof
        Fitted parameters for background model(s).
    exparam : ndarray or list thereof
        Fitted parameters for experiment model(s).
    Vuncert : :ref:`UncertQuant` or list thereof
        Uncertainty quanfitication for fitted dipolar signal(s).
    Puncert : :ref:`UncertQuant`
        Uncertainty quanfitication for fitted distance distribution.
    Buncert : :ref:`UncertQuant` or list thereof
        Uncertainty quanfitication for fitted background(s).
    ddparamUncert : :ref:`UncertQuant` 
        Uncertainty quanfitication for distribution parameters
    bgparamUncert : :ref:`UncertQuant` or list thereof
        Uncertainty quanfitication for background(s) parameters
    exparamUncert : :ref:`UncertQuant` or list thereof
        Uncertainty quanfitication for experiment(s) parameters
    scale : float int or list of float int
        Amplitude scale(s) of the dipolar signal(s).
    regparam : scalar
        Optimal regularization parameter value
    plot : callable
        Function to display the results. It will 
        display the fitted signals and distance distributions with
        confidence intervals. If requested, the function returns 
        the `matplotlib.axes` object as output. 
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

        fit = dl.fitsignal(V,t,r,dl.dd_gauss,dl.bg_hom3d,dl.ex_4pdeer)  


    Fit a 5pDEER signal with exponential background and Tikhonov regularization::

        fit = dl.fitsignal(V,t,r,'P',dl.bg_hom3d,dl.ex_5pdeer)
    
    
    Fit a 4pDEER stretched exponential background (no foreground)::
    
        fit = dl.fitsignal(V,t,r,None,dl.bg_strexp,dl.ex_4pdeer)  
    
    
    Fit a dipolar evolution function with Rician distribution::
    
        fit = dl.fitsignal(V,t,r,dl.dd_rice,None,None)
    
    
    Fit a 4pDEER form factor (no background) with bimodal Gaussian distribution:: 

        fit = dl.fitsignal(V,t,r,dl.dd_gauss2,None,dl.ex_4pdeer)   
   
     
    Fit a 4pDEER signal and a 5pDEER signal with same type of background::
    
        fit = fitsignal([V1,V2],t,r,'P',dl.bg_hom3d,[dl.ex_4pdeer,dl.ex_5pdeer])    
    """

    # Default optional settings
    normP = True

    # Make inputs into a list if just a singal is passed
    Vexp,t,bg_model,ex_model = ([var] if type(var) is not list else var for var in (Vexp,t,bg_model,ex_model))

    nSignals = len(Vexp)
    if len(t)!=nSignals:
        raise KeyError('The same number of signals V and time axes t must be provided.')
    for i in range(nSignals):
        if len(Vexp[i])!=len(t[i]):
            raise ValueError('V[{}] and t[{}] must have the same number of elements.'.format(i,i))
        if not np.all(np.isreal(Vexp[i])):
            raise ValueError('Input signals cannot be complex-valued')
    
    if len(bg_model)!=nSignals:
        bg_model = bg_model*nSignals
    
    if len(ex_model)!=nSignals:
        ex_model = ex_model*nSignals

    par0 = [[] if par0_i is None else par0_i for par0_i in par0]
    lb = [[] if lb_i is None else lb_i for lb_i in lb]
    ub = [[] if ub_i is None else ub_i for ub_i in ub]
    if type(par0) is not list or len(par0)!=3:
        raise TypeError('Initial parameters (7th input) must be a 3-element cell array.')
     
    # Get information about distance distribution parameters
    # ------------------------------------------------------
    parametricDistribution  = type(dd_model) is types.FunctionType
    parfreeDistribution = dd_model == 'P'
    includeForeground = np.array(dd_model is not None)

    par0_dd, lower_dd, upper_dd, N_dd = _getmodelparams(dd_model)

    if includeForeground and not parametricDistribution and not parfreeDistribution:
        raise TypeError('Distribution model (4th input) must either be a function handle, ''P'', or None.')

    # Get information about background parameters
    # ------------------------------------------------------
    isparamodel = np.array([type(model) is types.FunctionType for model in bg_model])
    includeBackground = np.array([model is not None for model in bg_model])

    par0_bg, lower_bg, upper_bg, N_bg = _getmodelparams(bg_model)

    isphenomenological = [_iserror(model,t,par0,1) for model,par0 in zip(bg_model,par0_bg)]


    if np.any(~isparamodel & includeBackground):
        raise TypeError('Background model (5th input) must either be a function handle, or None.')
    
    # Get information about experiment parameters
    # ------------------------------------------------------
    isparamodel = np.array([type(model) is types.FunctionType for model in ex_model])
    includeExperiment = np.array([model is not None for model in ex_model])

    par0_ex, lower_ex, upper_ex, N_ex = _getmodelparams(ex_model)

    if np.any(~isparamodel & includeExperiment):
        raise TypeError('Experiment models must either be a function handle, or None.')

    # Catch nonsensical situation
    if np.any(~includeForeground & ~includeExperiment):
        raise KeyError('Cannot fit anything without distribution model and without background model.')

    # Prepare full-parameter set
    # -------------------------------------------------------
    # Combine all initial, lower and upper bounds of parameters into vectors
    par0 = _parcombine(par0,[par0_dd,par0_bg,par0_ex])
    lb = _parcombine(lb,[lower_dd,lower_bg,lower_ex])
    ub = _parcombine(ub,[upper_dd,upper_bg,upper_ex])
    _checkbounds(lb,ub,par0)

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
        This function represent the core model of fitsignal, it takes the
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
                pathways = 1
            
            # Compute the multipathway-background
            B_ = dl.dipolarbackground(t[iSignal],pathways,Bfcn)
            # Compute the multipathway-kernel
            K_ = dl.dipolarkernel(t[iSignal],r,pathways)
            Ks.append(K_*B_[:,np.newaxis])
            Bs.append(B_)
            
        return Ks, Bs     
    # =========================================================================

    def splituq(full_uq,scales=1):
    # =========================================================================
        """ 
        Uncertainty quantification
        ---------------------------
        This function takes the full combined covariance matrix of the
        linear+nonlinear parameter sets and splits it into the different
        components of the model. 
        """
        # Pre-allocation
        paruq_bg,paruq_ex,Bfit_uq,Vfit_uq = ([],[],[],[])
        
        # Retrieve full covariance matrix
        covmat = full_uq.covmat

        Nparam = len(parfit_)
        paramidx = np.arange(Nparam)
        Pfreeidx = np.arange(Nparam,np.shape(covmat)[0])
        
        # Full parameter set uncertainty
        # -------------------------------
        subcovmat = covmat[np.ix_(paramidx,paramidx)]
        paruq = UncertQuant('covariance',parfit_,subcovmat,lb,ub)
        
        # Background parameters uncertainty
        # ---------------------------------
        for jj in range(nSignals):
            if includeBackground[jj]:
                bgsubcovmat  = paruq.covmat[np.ix_(bgidx[jj],bgidx[jj])]
                paruq_bg.append( UncertQuant('covariance',parfit_[bgidx[jj]],bgsubcovmat,lb[bgidx[jj]],ub[bgidx[jj]]))
            else:
                paruq_bg.append([None])
        
        # Experiment parameters uncertainty
        # ----------------------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                exsubcovmat  = paruq.covmat[np.ix_(exidx[jj],exidx[jj])]
                paruq_ex.append( UncertQuant('covariance',parfit_[exidx[jj]],exsubcovmat,lb[exidx[jj]],ub[exidx[jj]]))
            else:
                paruq_ex.append([None])
            
        # Distribution parameters uncertainty
        # ------------------------------------
        if parametricDistribution:
            ddsubcovmat  = paruq.covmat[np.ix_(ddidx,ddidx)]
            paruq_dd = UncertQuant('covariance',parfit_[ddidx],ddsubcovmat,lb[ddidx],ub[ddidx])
        else:
            paruq_dd = [None]
        
        # Distance distribution uncertainty
        # ----------------------------------
        nonneg = np.zeros_like(r)
        if parametricDistribution:
            Pfit_uq = paruq.propagate(Pfcn,nonneg,[])
        else:
            subcovmat = covmat[np.ix_(Pfreeidx,Pfreeidx)]
            Pfit_uq = UncertQuant('covariance',Pfit,subcovmat,nonneg,[])
        
        # Background uncertainty
        # -----------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                Bfit_uq.append( paruq.propagate(lambda par:multiPathwayModel(par)[1][jj]))
            else:
                Bfit_uq.append([None])
        
        # Dipolar signal uncertainty
        # --------------------------
        for jj in range(nSignals):
            if includeForeground and parametricDistribution:
                # Full parametric signal
                Vmodel = lambda par: scales[jj]*multiPathwayModel(par)[0][jj]@Pfcn(par[ddidx])
                Vfit_uq.append( paruq.propagate(Vmodel))
            elif includeForeground and np.all(~includeExperiment & ~includeBackground):
                # Dipola evolution function
                J = Ks[jj]
                Vcovmat = J@covmat@J.T
                Vfit_uq.append( UncertQuant('covariance',Vfit[jj],Vcovmat,[],[]))
            elif includeForeground:
                # Parametric signal with parameter-free distribution
                Jnonlin = Jacobian(lambda par: multiPathwayModel(par[paramidx])[0][jj]@Pfit,parfit_,lb[paramidx],ub[paramidx])
                J = np.concatenate((Jnonlin, Kfit[jj]),1)
                Vcovmat = J@covmat@J.T
                Vfit_uq.append( UncertQuant('covariance',Vfit[jj],Vcovmat,[],[]))
            else:
                Vfit_uq.append([None])

        return Vfit_uq,Pfit_uq,Bfit_uq,paruq_bg,paruq_ex,paruq_dd   
    # =========================================================================

    OnlyRegularization = np.all(~parametricDistribution & ~includeExperiment & ~includeBackground)
    OnlyParametric = not OnlyRegularization and (parametricDistribution or not includeForeground)

    if OnlyRegularization:
        
        # Use basic dipolar kernel
        Ks = [dl.dipolarkernel(ts,r) for ts in t]
        
        # Linear regularization fit
        fit = dl.fitregmodel(Vexp,Ks,r,regtype,regparam, weights=weights,uqanalysis=uqanalysis)
        Pfit = fit.P
        Pfit_uq = fit.uncertainty
        scales = np.atleast_1d(fit.scale)

        alphaopt = fit.regparam

        # Get fitted models
        Vfit = [scale*K@Pfit for K,scale in zip(Ks,scales)]
        Bfit = [np.ones_like(V) for V in Vexp]
        
        # No parameters
        parfit_ = np.asarray([None])
        if uqanalysis:
            Vfit_uq, Pfit_uq, Bfit_uq, paruq_bg, paruq_ex, paruq_dd = splituq(Pfit_uq)

    elif OnlyParametric:
        
        # Prepare the full-parametric model
        if includeForeground:
            Pfcn = lambda par: dd_model(r,par[ddidx])
        else:
            Pfcn = lambda _: np.ones_like(r)/np.trapz(np.ones_like(r),r)
        Vmodel =lambda par: [K@Pfcn(par) for K in multiPathwayModel(par)[0]]

        # Non-linear parametric fit
        fit = dl.fitparamodel(Vexp,Vmodel,par0,lb,ub,weights=weights,uqanalysis=uqanalysis)
        parfit_ = fit.param
        param_uq = fit.uncertainty
        scales = fit.scale
        alphaopt = None

        # Get fitted models
        Vfit = Vmodel(parfit_)
        _,Bfit = multiPathwayModel(parfit_)
        if includeForeground:
            Pfit = Pfcn(parfit_)
        else:
            Pfit = []
        if type(Vfit) is not list:
            Vfit = [Vfit]
        if type(scales) is not list:
            scales = [scales]
        Vfit = [V*scale for scale,V in zip(scales,Vfit) ]

        if uqanalysis:
            Vfit_uq, Pfit_uq, Bfit_uq, paruq_bg, paruq_ex, paruq_dd = splituq(param_uq,scales)
        
    else:
        
        # Non-negativity constraint on distributions
        lbl = np.zeros_like(r)
        
        prescales = [max(V) for V in Vexp]
        Vexp = [Vexp[i]/prescales[i] for i in range(nSignals)]

        # Separable non-linear least squares (SNNLS) 
        fit = dl.snlls(Vexp,lambda par: multiPathwayModel(par)[0],par0,lb,ub,lbl, reg=True,
                            regparam=regparam, lin_maxiter=[], uqanalysis=uqanalysis, lin_tol=[], weights=weights)
        parfit_ = fit.nonlin
        Pfit = fit.lin
        snlls_uq = fit.uncertainty
        alphaopt = fit.regparam
        scales = [prescales[i]*np.trapz(Pfit,r) for i in range(nSignals)]

        # Get the fitted models
        Kfit,Bfit = multiPathwayModel(parfit_)
        Vfit = [scale*K@Pfit for K,scale in zip(Kfit,scales)]

        if uqanalysis:
            Vfit_uq, Pfit_uq, Bfit_uq, paruq_bg, paruq_ex, paruq_dd = splituq(snlls_uq)
    
    # Normalize distribution
    # -----------------------
    scale = 1
    if normP and includeForeground:
        scale = np.trapz(Pfit,r)
        Pfit /= scale
        if uqanalysis:
            # scale CIs accordingly
            Pfit_uq_ = copy.deepcopy(Pfit_uq) # need a copy to avoid infite recursion on next step
            Pfit_uq.ci = lambda p: Pfit_uq_.ci(p)/scale

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
        modfituq['Bfit'] = [Bfit_uq[j] for j in range(nSignals)]
    else:
        paruq = dict()
        paruq['dd'] = UncertQuant('void')
        paruq['bg'] = UncertQuant('void')
        paruq['ex'] = UncertQuant('void')
        modfituq = dict()
        modfituq['Pfit'] = UncertQuant('void')
        modfituq['Vfit'] = UncertQuant('void')
        modfituq['Bfit'] = UncertQuant('void')

    Vfit_ = Vfit
    def _display_results():
    # =========================================================================
        _,axs = plt.subplots(nSignals+1,figsize=[7,3+3*nSignals])
        for i in range(nSignals):
            # Plot the signal
            axs[i].plot(t[i],Vexp[i],'.',color='grey',alpha=0.5)
            axs[i].plot(t[i],Vfit_[i],'tab:blue')
            if uqanalysis:
                # Get confidence intervals for the signal
                Vci95 = Vfit_uq[i].ci(95)
                Vci50 = Vfit_uq[i].ci(50)
                axs[i].fill_between(t[i],Vci95[:,0], Vci95[:,1],facecolor='tab:blue',linestyle='None',alpha=0.2)
                axs[i].fill_between(t[i],Vci50[:,0], Vci50[:,1],facecolor='tab:blue',linestyle='None',alpha=0.4)
            axs[i].grid(alpha=0.3)
            axs[i].set_xlabel('Time [μs]')
            axs[i].set_ylabel('V[{}]'.format(i))
            axs[i].legend(('Data','Fit','95%-CI','50%-CI'))

        # Plot the distribution
        axs[nSignals].plot(r,Pfit,'tab:blue')
        if uqanalysis:
            # Get confidence intervals for the distance distribution
            Pci95 = Pfit_uq.ci(95)
            Pci50 = Pfit_uq.ci(50)
            axs[nSignals].fill_between(r,Pci95[:,0], Pci95[:,1],facecolor='tab:blue',linestyle='None',alpha=0.2)
            axs[nSignals].fill_between(r,Pci50[:,0], Pci50[:,1],facecolor='tab:blue',linestyle='None',alpha=0.4)
        axs[nSignals].set_xlabel('Distance [nm]')
        axs[nSignals].set_ylabel('P [nm⁻¹]')
        axs[nSignals].legend(('Fit','95%-CI','50%-CI'))
        axs[nSignals].grid(alpha=0.3)
        plt.tight_layout()
        plt.show()
        return axs
    # =========================================================================


    if verbose:
        print('----------------------------------------------------------------------------')
        print('Goodness of fit')
        for i in range(nSignals):
            print('  Vexp[{}]: chi2 = {:4f}  RMSD  = {:4f}'.format(i,stats[i]['chi2red'],stats[i]['rmsd']))
        print('----------------------------------------------------------------------------')
        print('Fitted parameters and 95%-confidence intervals')
        pstr = "  {}[{:d}]:   {:5.7f}  ({:.7f}, {:.7f})  {} ({})"
        if len(parfit['dd'])>0:
            print()
            info = dd_model()
            if uqanalysis:
                ci = paruq['dd'].ci(95)
            else:
                ci = np.full((len(parfit['dd']),2),np.nan)
            for j in range(len(parfit['dd'])):
                c = parfit['dd'][j]
                print(pstr.format('ddparam',j,c,ci[j,0],ci[j,1],info['Parameters'][j],info['Units'][j]))

        for i in range(nSignals):
            print('Vfit[{}]:'.format(i))
            if includeBackground[i]:
                if len(parfit['bg'])>0:
                    info = bg_model[i]()
                if uqanalysis:
                    ci = paruq['bg'][i].ci(95)
                else:
                    ci = np.full((len(parfit['bg'][i]),2),np.nan)
                for j in range(len(parfit['bg'][i])):
                        c = parfit['bg'][i][j]
                        print(pstr.format('bgparam',j,c,ci[j,0],ci[j,1],info['Parameters'][j],info['Units'][j]))
            if includeExperiment[i]:
                if len(parfit['ex'])>0:
                    info = ex_model[i]()
                if uqanalysis:
                    ci = paruq['ex'][i].ci(95)
                else:
                    ci = np.full((len(parfit['ex'][i]),2),np.nan)
                for j in range(len(parfit['ex'][i])):
                    c = parfit['ex'][i][j]
                    print(pstr.format('exparam',j,c,ci[j,0],ci[j,1],info['Parameters'][j],info['Units'][j]))
        print('----------------------------------------------------------------------------')


    # Return numeric arrays and not lists if there is only one signal
    if nSignals==1:
        Vfit,Bfit = (var[0] if type(var) is list else var for var in [Vfit,Bfit])
        for subset in ('ex','bg'):
            parfit[subset] = parfit[subset][0]
            if uqanalysis:
                paruq[subset] = paruq[subset][0]
        for subset in ('Vfit','Bfit'):
            if uqanalysis:
                modfituq[subset] = modfituq[subset][0]
        if not isempty(stats):
            stats = stats[0]

    return FitResult(V=Vfit, P=Pfit, B=Bfit, exparam=parfit['ex'], bgparam=parfit['bg'],
                      ddparam=parfit['dd'], Vuncert = modfituq['Vfit'], Puncert = modfituq['Pfit'],
                      Buncert = modfituq['Bfit'], exparamUncert = paruq['ex'], bgparamUncert = paruq['bg'],
                      ddparamUncert = paruq['dd'], regparam = alphaopt, plot=_display_results, scale=scales,  stats=stats, cost=fit.cost,
                      residuals=fit.residuals, success=fit.success)

            
def _getmodelparams(models):
# ==============================================================================
    "Extract parameter info from parametric models"
    if type(models) is not list:
        models = [models]

    par0,lo,up,N = ([],[],[],[])
    for model in models:
        if type(model) is types.FunctionType:
            info = model()
            par0.append(info['Start'])
            lo.append(info['Lower'])
            up.append(info['Upper'])
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
    "Combine parameter subsets"
    
    for k in range(3):
        p[k] = np.atleast_1d(p[k])
        d[k] = np.atleast_1d(d[k])
        if isempty(p[k]):
            p[k] = d[k].flatten()
        if type(p[k]) is list:
            p[k] = p[k].flatten()
    pvec = np.concatenate(p,axis=0)
    return pvec
# ==============================================================================

def _checkbounds(lb,ub,par0):
# ==============================================================================
    "Boundary checking"
    nParams = len(par0)
    if len(lb)!=nParams:
        raise ValueError('Lower bounds and initial values of parameters must have the same number of elements.')
    if len(ub)!=nParams:
        raise ValueError('Lower bounds and initial values of parameters must have the same number of elements.')
    if np.any(lb>ub):
        raise ValueError('Lower bounds cannot be larger than upper bounds.')
    if np.any(lb>par0) or np.any(par0>ub):
        raise ValueError('Inital values for parameters must lie between lower and upper bounds.')
# ==============================================================================

def _iserror(func, *args, **kw):
# ==============================================================================
    try:
        func(*args, **kw)
        return False
    except Exception:
        return True
# ==============================================================================
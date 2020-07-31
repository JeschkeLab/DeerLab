# fitsignal.py - Dipolar signal fit function
# ---------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
import copy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import deerlab as dl
from deerlab import uqst
from deerlab.bg_models import bg_hom3d
from deerlab.ex_models import ex_4pdeer
from deerlab.utils import isempty, goodness_of_fit, jacobianest

def fitsignal(Vexp,t,r,dd_model='P',bg_model=bg_hom3d,ex_model=ex_4pdeer,par0=[],lb=[],ub=[], weights=1, uqanalysis=True, display = False):
    """Fits a dipolar model to the experimental signal V with time axis t, using
    distance axis r. The model is specified by the distance distribution (dd),
    the background (bg), and the experiment (ex).

    If multiple signals (V1,V2,...) and their corresponding time axes (t1,t2,...)
    are given, they will be fitted globally with a single distance distribution (dd).
    For each signal, a specific background (bg1,bg2,...) and experiment (ex1,ex2)
    models can be assigned.

    This function can handle both parametric and non-parametric distance distribution models.

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

        (1) A string 'P' to indicate a non-parametric distribution
        (2) A callable function of a DeerLab parametric distribution model (e.g. dd_gauss).
        (3) A string 'none' to indicate no distribution, i.e. only background.

        The default is 'P'.

    bg : callable or string
        Background model, the following modes are allowed:

        (1) A callable function of a DeerLab parametric background model (e.g. bg_hom3d).
        (2) A string 'none' to indicate no background decay.

        The default is bg_hom3d.

    ex : callable or string
        Experiment model, the following modes are allowed:

        (1) Function handle to experiment model (e.g. @ex_4pdeer)
        (2) 'none' to indicate simple dipolar oscillation (mod.depth = 1)

        The default is ex_4pdeer.

    par0 : list of array_like
        Starting parameter values. Must be a 3-element list ``[par0_dd,par0_bg,par0_ex]``
        containing the start values of the distribution, background and experiment models, in that order.
        If a model does not require parameters or are to be determined automatically it must be specified 
        as an empty list ``[]``. The default is ``par0=[[],[],[]]``.
    lb : list of array_like
        Lower bounds for parameters.  Must be a 3-element list ``[lb_dd,lb_bg,lb_ex]``
        containing the start values of the distribution, background and experiment models, in that order.
        If a model does not require parameters or are to be determined automatically it must be specified 
        as an empty list ``[]``. The default is ``lb=[[],[],[]]``.
    ub : list of array_like
        Upper bounds for parameters, Must be a 3-element list ``[ub_dd,ub_bg,ub_ex]``
        containing the start values of the distribution, background and experiment models, in that order.
        If a model does not require parameters or are to be determined automatically it must be specified 
        as an empty list ``[]``. The default is ``ub=[[],[],[]]``.

    Returns
    -------
    Vfit :  ndarray or list of ndarrays
        Fitted time-domain signal(s).
    Pfit : ndarray
        Fitted distance distribution.
    Bfit : ndarray or list of ndarrays
        Fitted background decay(s).
    parfit : dict
        Fitted parameters:
        
        * parfit['dd'] - Fitted parameters for distance distribution model
        * parfit['bg'] - Fitted parameters for background model
        * parfit['ex'] - Fitted parameters for experiment model
    modfituq : dict
        Uncertainty quanfitications of the fits:
        
        * modfituq['Vfit'] - Uncertainty quanfitication for fitted dipolar signal
        * modfituq['Pfit'] - Uncertainty quanfitication for fitted distance distribution
        * modfituq['Bfit'] - Uncertainty quanfitication for fitted background
    paruq : dict
        Uncertainty quanfitications of the parameters:
        
        * paruq['dd'] - Uncertainty quanfitication for distribution parameters
        * paruq['bg'] - Uncertainty quanfitication for background parameters
        * paruq['ex'] - Uncertainty quanfitication for experiment parameters
    stats : dict or list of dict
        Goodness of fit statistical estimators:

        * stats['chi2red'] - Reduced \chi^2 test
        * stats['r2'] - R^2 test
        * stats['rmsd'] - Root-mean squared deviation (RMSD)
        * stats['aic'] - Akaike information criterion
        * stats['aicc'] - Corrected Akaike information criterion
        * stats['bic'] - Bayesian information criterion

    Other Parameters
    ----------------
    weights : array_like 
        Array of weighting coefficients for the individual signals in global fitting,
        the default is all weighted equally.
        If not specified all datasets are weighted equally.
    uqanalysis : boolean
        Enable/disable the uncertainty quantification analysis, by default it is enabled.  
    display : boolean
        Enable/Disable plotting and printing of results with matplotlib library, by default it is disabled. 
        If enabled the command ``matlplotlib.pyplot.show()`` must be added to the script calling fitsignal 
        at the end to show the figure containing the plot.

    Examples
    --------
    Fit a 4pDEER signal with homogenous 3D background with Gaussian distribution::

        Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats = fitsignal(V,t,r,dd_gauss,bg_hom3d,ex_4pdeer)  


    Fit a 5pDEER signal with exponential background and Tikhonov regularization::

        Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats = fitsignal(V,t,r,'P',bg_hom3d,ex_5pdeer)        
    
    
    Fit a 4pDEER stretched exponential background (no foreground)::
    
        Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats = fitsignal(V,t,r,'none',bg_strexp,ex_4pdeer)  
    
    
    Fit a dipolar evolution function with Rician distribution::
    
        Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats = fitsignal(V,t,r,dd_rice,'none','none')          
    
    
    Fit a 4pDEER form factor (no background) with bimodal Gaussian distribution:: 

        Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats = fitsignal(V,t,r,dd_gauss2,'none',ex_4pdeer)   
   
     
    Fit a 4pDEER signal and a 5pDEER signal with same type of background::
    
        Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats = fitsignal([V1,V2],t,r,'P',bg_hom3d,[ex_4pdeer,ex_5pdeer])    
    """

    # Default optional settings
    regtype = 'tikhonov'
    regparam = 'aic'
    DisplayResults = display
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

    if isempty(lb):
        lb = [[],[],[]]
    if isempty(ub):
        ub = [[],[],[]]
    if isempty(par0):
        par0 = [[],[],[]]
    elif type(par0) is not list or len(par0)!=3:
        raise TypeError('Initial parameters (7th input) must be a 3-element cell array.')
     
    # Get information about distance distribution parameters
    # ------------------------------------------------------
    parametricDistribution  = type(dd_model) is types.FunctionType
    parfreeDistribution = dd_model is 'P'
    includeForeground = np.array(dd_model is not 'none')

    par0_dd, lower_dd, upper_dd, N_dd = _getmodelparams(dd_model)

    if includeForeground and not parametricDistribution and not parfreeDistribution:
        raise TypeError('Distribution model (4th input) must either be a function handle, ''P'', or ''none''.')

    # Get information about background parameters
    # ------------------------------------------------------
    isparamodel = np.array([type(model) is types.FunctionType for model in bg_model])
    includeBackground = np.array([model is not 'none' for model in bg_model])

    par0_bg, lower_bg, upper_bg, N_bg = _getmodelparams(bg_model)

    if np.any(~isparamodel & includeBackground):
        raise TypeError('Background model (5th input) must either be a function handle, or ''none''.')
    
    # Get information about experiment parameters
    # ------------------------------------------------------
    isparamodel = np.array([type(model) is types.FunctionType for model in ex_model])
    includeExperiment = np.array([model is not 'none' for model in ex_model])

    par0_ex, lower_ex, upper_ex, N_ex = _getmodelparams(ex_model)

    if np.any(~isparamodel & includeExperiment):
        raise TypeError('Experiment models must either be a function handle, or ''none''.')

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

            # Prepared background basis function
            if includeBackground[iSignal]:
                Bfcn = lambda t,lam: bg_model[iSignal](t,bg_par,lam)
            else:
                Bfcn = lambda _,__: np.ones_like(Vexp[iSignal])
            
            # Get pathway information
            if includeExperiment[iSignal]:
                pathinfo = ex_model[iSignal](ex_par)
            else:
                pathinfo = 1
            
            # Compute the multipathway-background
            B_ = dl.dipolarbackground(t[iSignal],pathinfo,Bfcn)
            # Compute the multipathway-kernel
            K_ = dl.dipolarkernel(t[iSignal],r,pathinfo)
            Ks.append(K_*B_[:,np.newaxis])
            Bs.append(B_)
            
        return Ks, Bs     
    # =========================================================================

    def splituq(full_uq):
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
        paruq = uqst('covariance',parfit_,subcovmat,lb,ub)
        
        # Background parameters uncertainty
        # ---------------------------------
        for jj in range(nSignals):
            if includeBackground[jj]:
                bgsubcovmat  = paruq.covmat[np.ix_(bgidx[jj],bgidx[jj])]
                paruq_bg.append( uqst('covariance',parfit_[bgidx[jj]],bgsubcovmat,lb[bgidx[jj]],ub[bgidx[jj]]))
            else:
                paruq_bg.append([])
        
        # Experiment parameters uncertainty
        # ----------------------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                exsubcovmat  = paruq.covmat[np.ix_(exidx[jj],exidx[jj])]
                paruq_ex.append( uqst('covariance',parfit_[exidx[jj]],exsubcovmat,lb[exidx[jj]],ub[exidx[jj]]))
            else:
                paruq_ex.append([])
            
        # Distribution parameters uncertainty
        # ------------------------------------
        if parametricDistribution:
            ddsubcovmat  = paruq.covmat[np.ix_(ddidx,ddidx)]
            paruq_dd = uqst('covariance',parfit_[ddidx],ddsubcovmat,lb[ddidx],ub[ddidx])
        else:
            paruq_dd = []
        
        # Distance distribution uncertainty
        # ----------------------------------
        nonneg = np.zeros_like(r)
        if parametricDistribution:
            Pfit_uq = paruq.propagate(Pfcn,nonneg,[])
        else:
            subcovmat = covmat[np.ix_(Pfreeidx,Pfreeidx)]
            Pfit_uq = uqst('covariance',Pfit,subcovmat,nonneg,[])
        
        # Background uncertainty
        # -----------------------
        for jj in range(nSignals):
            if includeExperiment[jj]:
                Bfit_uq.append( paruq.propagate(lambda par:multiPathwayModel(par)[1][jj]))
            else:
                Bfit_uq.append([])
        
        # Dipolar signal uncertainty
        # --------------------------
        for jj in range(nSignals):
            if includeForeground and parametricDistribution:
                # Simple parametric model error propagation
                Vmodel = lambda par: multiPathwayModel(par)[0][jj]@Pfcn(par[ddidx])
                Vfit_uq.append( paruq.propagate(Vmodel))
            elif includeForeground:
                # Use special structure to speed up propagation for
                # parameter-free case instead of .propagate()
                J = np.concatenate((jacobianest(lambda par: multiPathwayModel(par[paramidx])[0][jj]@Pfit,parfit_)[0], Kfit[jj]),1)
                Vcovmat = J@covmat@J.T
                Vfit_uq.append( uqst('covariance',Vfit[jj],Vcovmat,[],[]))
            else:
                Vfit_uq.append([])

        return Vfit_uq,Pfit_uq,Bfit_uq,paruq_bg,paruq_ex,paruq_dd   
    # =========================================================================

    OnlyRegularization = np.all(~includeExperiment & ~includeBackground)
    OnlyParametric = not OnlyRegularization and (parametricDistribution or not includeForeground)

    if OnlyRegularization:
        
        # Use basic dipolar kernel
        Ks = [dl.dipolarkernel(ts,r) for ts in t]
        
        # Linear regularization fit
        Pfit,Pfit_uq = dl.fitregmodel(Vexp,Ks,r,regtype,regparam, weights=weights,uqanalysis=uqanalysis)
        
        # Get fitted models
        Vfit = [K@Pfit for K in Ks]
        Bfit = [np.ones_like(V) for V in Vexp]
        
        # No parameters
        parfit_,Vfit_uq,Bfit_uq,paruq_bg,paruq_dd,paruq_ex = np.atleast_1d([[[]] for _ in range(6)])
        
    elif OnlyParametric:
        
        # Prepare the full-parametric model
        if includeForeground:
            Pfcn = lambda par: dd_model(r,par[ddidx])
        else:
            Pfcn = lambda _: np.ones_like(r)
        Vmodel =lambda par: [K@Pfcn(par) for K in multiPathwayModel(par)[0]]

        # Non-linear parametric fit
        parfit_,param_uq,_ = dl.fitparamodel(Vexp,Vmodel,par0,lb,ub,weights=weights,uqanalysis=uqanalysis)
        
        # Get fitted models
        Vfit = Vmodel(parfit_)
        _,Bfit = multiPathwayModel(parfit_)
        if includeForeground:
            Pfit = Pfcn(parfit_)
        else:
            Pfit = []
        if type(Vfit) is not list:
            Vfit = [Vfit]
        if uqanalysis:
            Vfit_uq, Pfit_uq, Bfit_uq, paruq_bg, paruq_ex, paruq_dd = splituq(param_uq)
        
    else:
        
        # Non-negativity constraint on distributions
        lbl = np.zeros_like(r)
        
        # Separable non-linear least squares (SNNLS) 
        parfit_, Pfit, snlls_uq,_ = dl.snlls(Vexp,lambda par: multiPathwayModel(par)[0],par0,lb,ub,lbl, regparam=regparam,linMaxIter=[],uqanalysis=uqanalysis,linTolFun=[],weights=weights)

        # Get the fitted models
        Kfit,Bfit = multiPathwayModel(parfit_)
        Vfit = [K@Pfit for K in Kfit]
        
        if uqanalysis:
            Vfit_uq, Pfit_uq, Bfit_uq, paruq_bg, paruq_ex, paruq_dd = splituq(snlls_uq)
    
    # Normalize distribution
    # -----------------------
    if normP and includeForeground:
        Pnorm = np.trapz(Pfit,r)
        Pfit /= Pnorm
        if uqanalysis:
            # scale CIs accordingly
            Pfit_uq_ = copy.deepcopy(Pfit_uq) # need a copy to avoid infite recursion on next step
            Pfit_uq.ci = lambda p: Pfit_uq_.ci(p)/Pnorm

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
        paruq = []
        modfituq = []


    def _display_results():
    # =========================================================================
        _,axs = plt.subplots(nSignals+1,1)
        for i in range(nSignals):
            # Get confidence intervals for the signal
            Vci95 = Vfit_uq[i].ci(95)
            Vci50 = Vfit_uq[i].ci(50)
            # Plot the signal
            axs[i].plot(t[i],Vexp[i],'k.',t[i],Vfit[i],'r')
            axs[i].fill_between(t[i],Vci95[:,0], Vci95[:,1],facecolor='r',linestyle='None',alpha=0.2)
            axs[i].fill_between(t[i],Vci50[:,0], Vci50[:,1],facecolor='r',linestyle='None',alpha=0.5)
            plt.tight_layout()
            axs[i].grid()
            axs[i].set_xlabel('Time [$\\mu s$]')
            axs[i].set_ylabel('V[{}]'.format(i))
            axs[i].legend(('Data','Fit','95%-CI','50%-CI'))

        # Get confidence intervals for the distance distribution
        Pci95 = Pfit_uq.ci(95)
        Pci50 = Pfit_uq.ci(50)
        # Plot the distribution
        axs[nSignals].plot(r,Pfit,'k')
        axs[nSignals].fill_between(r,Pci95[:,0], Pci95[:,1],facecolor='r',linestyle='None',alpha=0.2)
        axs[nSignals].fill_between(r,Pci50[:,0], Pci50[:,1],facecolor='r',linestyle='None',alpha=0.5)
        axs[nSignals].set_xlabel('Distance [nm]')
        axs[nSignals].set_ylabel('P [nm$^{-1}$]')
        axs[nSignals].legend(('Fit','95%-CI','50%-CI'))
        plt.tight_layout()
        axs[nSignals].grid()
        
        plt.show(block=False)
        print('----------------------------------------------------------------------------')
        print('Goodness of fit')
        for i in range(nSignals):
            print('  Vexp[{}]: chi2 = {:4f}  RMSD  = {:4f}'.format(i,stats[i]['chi2red'],stats[i]['rmsd']))
        print('----------------------------------------------------------------------------')
        print('Fitted parameters and 95%-confidence intervals')
        str = "  parfit['{}'][{:d}][{:d}]:   {:.7f}  ({:.7f}, {:.7f})  {} ({})"
        if len(parfit['dd'])>0:
            info = dd_model()
            ci = paruq['dd'].ci(95)
            for j in range(len(parfit['dd'])):
                c = parfit['dd'][j]
                print(str.format('dd',1,j,c,ci[j,0],ci[j,1],info['Parameters'][j],info['Units'][j]))

        if len(parfit['bg'])>0:
            for i in range(nSignals):
                info = bg_model[i]()
                ci = paruq['bg'][i].ci(95)
                for j in range(len(parfit['bg'][i])):
                    c = parfit['bg'][i][j]
                    print(str.format('bg',i,j,c,ci[j,0],ci[j,1],info['Parameters'][j],info['Units'][j]))

        if len(parfit['ex'])>0:
            for i in range(nSignals):
                info = ex_model[i]()
                ci = paruq['ex'][i].ci(95)
                for j in range(len(parfit['ex'][i])):
                    c = parfit['ex'][i][j]
                    print(str.format('ex',i,j,c,ci[j,0],ci[j,1],info['Parameters'][j],info['Units'][j]))
        print('----------------------------------------------------------------------------')
    # =========================================================================


    # Plotting
    # --------
    if DisplayResults:
        _display_results()

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

    return Vfit, Pfit, Bfit, parfit, modfituq, paruq, stats


            
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

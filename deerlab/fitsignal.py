#
# FITSIGNAL  Fit model to dipolar time-domain trace
#
#   [Vfit,Pfit,Bfit,parfit,modfitci,parci,stats] = FITSIGNAL(V,t,r,dd,bg,ex,par0,lb,ub)
#   __ = FITSIGNAL(V,t,r,dd,bg,ex,par0)
#   __ = FITSIGNAL(V,t,r,dd,bg,ex)
#   __ = FITSIGNAL(V,t,r,dd,bg)
#   __ = FITSIGNAL(V,t,r,dd)
#   __ = FITSIGNAL(V,t,r)
#   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__},par0,lb,ub)
#   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__},par0)
#   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__})
#   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},ex)
#   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd)
#   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r)
#   __ = FITSIGNAL(___,'Property',Values,___)
#
#   Fits a dipolar model to the experimental signal V with time axis t, using
#   distance axis r. The model is specified by the distance distribution (dd),
#   the background (bg), and the experiment (ex).
#
#   If multiple signals (V1,V2,...) and their corresponding time axes (t1,t2,...)
#   are given, they will be fitted globally with a single distance distribution (dd).
#   For each signal, a specific background (bg1,bg2,...) and experiment (ex1,ex2)
#   models can be assigned.
#
#   FITSIGNAL can handle both parametric and non-parametric distance
#   distribution models.
#
#  Input:
#    V      time-domain signal to fit (N-element vector)
#    t      time axis, in microseconds (N-element vector)
#    r      distance axis, in nanometers (M-element vector)
#    dd     distance distribution model (default 'P')
#           - 'P' to indicate non-parametric distribution
#           - function handle to parametric distribution model (e.g. @dd_gauss)
#           - 'none' to indicate no distribution, i.e. only background
#    bg     background model (default @bg_hom3d)
#           - function handle to parametric background model (e.g. @bg_hom3d)
#           - 'none' to indicate no background decay
#    ex     experiment model (default @ex_4pdeer)
#           - function handle to experiment model (e.g. @ex_4pdeer)
#           - 'none' to indicate simple dipolar oscillation (mod.depth = 1)
#    par0   starting parameters, 3-element cell array {par0_dd,par0_bg,par0_ex}
#           default: {[],[],[]} (automatic choice)
#    lb     lower bounds for parameters, 3-element cell array (lb_dd,lb_bg,lb_ex)
#           default: {[],[],[]} (automatic choice)
#    ub     upper bounds for parameters, 3-element cell array (ub_dd,ub_bg,ub_ex)
#           default: {[],[],[]} (automatic choice)
#
#  Output:
#    Vfit     fitted time-domain signal
#    Pfit     fitted distance distribution
#    Bfit     fitted background decay
#    parfit   structure with fitted parameters
#                .dd  fitted parameters for distance distribution model
#                .bg  fitted parameters for background model
#                .ex  fitted parameters for experiment model
#    modfitci structure with confidence intervals for Vfit, Bfit and Pfit
#    parci    structure with confidence intervals for parameter, similar to parfit
#    stats    goodness of fit statistical estimators, N-element structure array
#
#  Name-value pairs:
#
#   'Rescale'       - Enable/Disable optimization of the signal scale
#   'normP'         - Enable/Disable re-normalization of the fitted distribution
#   'Display'       - Enable/Disable plotting and printing of results
#   See "help snnls" for more options.
#
# Example:
#    Vfit = fitsignal(Vexp,t,r,@dd_gauss,@bg_hom3d,@ex_4pdeer)
#

# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
import copy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from plotly.subplots import make_subplots
import plotly.graph_objects as go


import deerlab as dl
from deerlab import uqst
from deerlab.bg_models import bg_hom3d
from deerlab.utils import isempty, goodness_of_fit, jacobianest

def fitsignal(Vexp,t,r,dd_model='P',bg_model=bg_hom3d,ex_model=[],par0=[],lb=[],ub=[], uqanalysis=True, display = False):

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
    ddidx = np.arange(N_dd)
    N_bg,N_ex = np.atleast_1d(N_bg,N_ex)
    bgidx,exidx = ([],[])
    for ii in range(nSignals):
        bgidx.append(N_dd + sum(N_bg[np.arange(0,ii)]) + np.arange(0,N_bg[ii]) )
        exidx.append(N_dd + sum(N_bg) + sum(N_ex[np.arange(0,ii)]) + np.arange(0,N_ex[ii]) )

    # Fit the dipolar multi-pathway model
    # -----------------------------------
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
            # Get parameter subsets for this signal
            ex_par = par[exidx[iSignal]]
            bg_par = par[bgidx[iSignal]]
            
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

    OnlyRegularization = np.any(~includeExperiment & ~includeBackground)
    OnlyParametric = not OnlyRegularization and (parametricDistribution or not includeForeground)

    if OnlyRegularization:
        
        # Use basic dipolar kernel
        Ks = [dl.dipolarkernel(ts,r) for ts in t]
        
        # Linear regularization fit
        Pfit,Pfit_uq = dl.fitregmodel(Vexp,Ks,r,regtype,regparam)
        
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
        parfit_,param_uq,_ = dl.fitparamodel(Vexp,Vmodel,par0,lb,ub)
        
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
        parfit_, Pfit, snlls_uq,_ = dl.snlls(Vexp,lambda par: multiPathwayModel(par)[0],par0,lb,ub,lbl, regparam=regparam,linMaxIter=[],linTolFun=[])

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

    def _display_results():
    # =========================================================================
        _,axs = plt.subplots(1, nSignals+1)
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

    if len(models)==1:
        par0 = par0[0]
        lo = lo[0]
        up = up[0]
        N = N[0]
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

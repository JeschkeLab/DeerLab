import numpy as np
import copy
import deerlab as dl
from deerlab.utils import isempty, jacobianest, hccm

def fitmultimodel(V,Kmodel,r,model,maxModels,method='aic',lb=[],ub=[], weights=1, normP = True, uqanalysis=True,**kwargs):

    # Ensure that all arrays are numpy.nparray
    lb,ub,r = [np.atleast_1d(var) for var in (lb,ub,r)]
    
    # Parse multiple datsets and non-linear operators into a single concatenated vector/matrix
    V, Kmodel, weights, Vsubsets = dl.utils.parse_multidatasets(V, Kmodel, weights)

    # Check kernel model
    if Kmodel is callable:
        # If callable, determine how many parameters the model requires
        nKparam = 0
        notEnoughParam = True
        while notEnoughParam:
            nKparam = nKparam + 1
            try:
                Kmodel(np.random.uniform(nKparam,1))
                notEnoughParam = False
            except:
                notEnoughParam = True
    else:
        # If the kernel is just a matrix make it a callable without parameters
        nKparam = 0
        K = copy.deepcopy(Kmodel) # need a copy to avoid infite recursion on next step
        Kmodel = lambda _: K
    
    # Parse boundaries
    if isempty(lb):
        lb = [[],[]]
    if isempty(ub):
        ub = [[],[]]
    # Unpack the two boundaries subsets
    lbK,lb0 = lb
    ubK,ub0 = ub
    if len(lbK) is not nKparam or len(ubK) is not nKparam:
        raise ValueError('The upper/lower bounds of the kernel parameters must be ',nKparam,'-element arrays')

    # Extract information about the model
    info = model()
    par00 = info['Start']
    nparam = len(par00)
    if isempty(lb0):
        lb0 = info['Lower']
    if isempty(ub0):
        ub0 = info['Upper']
    paramnames = info['Parameters']

    areCenterDistances = [str in ['Center','Location'] for str in paramnames]
    if any(areCenterDistances):
        # If the center of the basis function is a parameter limit it 
        # to the distance axis range (stabilizes parameter search)
        ub0[areCenterDistances] = max(r)
        lb0[areCenterDistances] = min(r)

    # Ensure that all arrays are numpy.nparray
    lb0,ub0,lbK,ubK = [np.atleast_1d(var) for var in (lb0,ub0,lbK,ubK)]

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
        K = Kmodel(par[1:nKparam])
        subset = np.atleast_1d(nKparam)
        Knonlin = np.zeros((K.shape[0],Nmodels))
        for iModel in range(Nmodels):
            subset = np.arange(subset[-1], subset[-1]+nparam)
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
        subset = np.atleast_1d(nKparam)
        Pfit = 0
        for iModel in range(len(linpar)):
            subset = np.arange(subset[-1], subset[-1]+nparam)
            # Get Gauss basis functions
            Pfit = Pfit + linpar[iModel]*model(r,nlinpar[subset])
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
        SquaredSumRes = np.sum((V - Vfit)**2)
        logprob = weights*N*np.log(SquaredSumRes/N)
        # Compute the estimators
        rmsd = weights*np.sqrt(1/N*SquaredSumRes)
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

    # Pre-allocate containers
    Vfit,Pfit,plin_,pnonlin_,nlin_ub_,nlin_lb_,lin_ub_,lin_lb_ = ([] for _ in range(8))
    logest = []

    # Loop over number of components in model
    # =======================================
    for Nmodels in np.arange(1,maxModels+1):
        
        # Prepare non-linear model with N-components
        # ===========================================
        Knonlin = lambda par: nonlinmodel(par,Nmodels)
        
        # Box constraints for the model parameters (non-linear parameters)
        nlin_lb = np.matlib.repmat(lb0,1,Nmodels)
        nlin_ub = np.matlib.repmat(ub0,1,Nmodels)
        nlin_ub = np.matlib.repmat(par00,1,Nmodels)

        # Add the box constraints on the non-linear kernel parameters
        nlin_lb = np.concatenate((lbK, nlin_lb),axis=None)
        nlin_ub = np.concatenate((ubK, nlin_ub),axis=None)
        
        # Start values of non-linear parameters
        np.random.seed(0)
        par0 = np.random.uniform(size=(len(nlin_lb)),low=nlin_lb, high=nlin_ub)
        
        # Box constraints for the components amplitudes (linear parameters)
        lin_lb = np.zeros(Nmodels)        # Non-negativity constraint
        lin_ub = np.full(Nmodels, np.inf) # Unbounded
        
        # Separable non-linear least-squares (SNLLS) fit
        # ===============================================
        pnonlin,plin,_,stats = dl.snlls(V,Knonlin,par0,nlin_lb,nlin_ub,lin_lb,[], penalty=False, uqanalysis=False,**kwargs)

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
        lin_lb_.append(lin_lb)

    Peval = Pfit

    # Select the optimal model
    # ========================
    fcnals = logest[method]
    idx = np.argmin(fcnals)
    Nopt = idx+1
    Pfit = Pfit[idx]
    Vfit = Vfit[idx]
    pnonlin = pnonlin_[idx]
    plin = plin_[idx]
    nlin_lb = nlin_lb_[idx]
    nlin_ub = nlin_ub_[idx]
    lin_lb = lin_lb_[idx]
    lin_ub = lin_ub_[idx]

    # Package the fitted parameters
    # =============================
    fitparam = [[],[],[]]
    fitparam[0] = pnonlin[0:nKparam] # Kernel parameters
    fitparam[1] = pnonlin[nKparam:-1] # Components parameters
    fitparam[2] = plin # Components amplitudes

    # Uncertainty quantification analysis (if requested)
    # ==================================================
    if uqanalysis:
        # Compute the residual vector
        Knonlin = lambda par: nonlinmodel(par,Nopt)
        res = weights*(Vfit - V)

        # Compute the Jacobian
        Jnonlin,_ = jacobianest(lambda p: weights*Knonlin(p)@plin, pnonlin)
        Jlin = weights*Knonlin(pnonlin)
        J = np.concatenate((Jnonlin, Jlin),1)
        
        # Estimate the heteroscedasticity-consistent covariance matrix
        covmatrix = hccm(J,res,'HC1')
        
        # Construct uncertainty quantification structure for fitted parameters
        paramuq = dl.uqst('covariance',np.concatenate((pnonlin, plin)),covmatrix,np.concatenate((nlin_lb, lin_lb)),np.concatenate((nlin_ub, lin_ub)))
        
        P_subset = np.arange(0, nKparam+nparam*Nopt)
        amps_subset = np.arange(P_subset[-1]+1, P_subset[-1]+1+Nopt)
        Puq = paramuq.propagate(lambda p: Pmodel(p[P_subset],p[amps_subset]), np.zeros(len(r)))
    else: 
        Puq = []
    # Goodness of fit
    # ===============
    Ndof = nKparam + nparam + Nopt
    #stats = cellfun(@(V,Vfit)gof(V,Vfit,Ndof),V,Vfit)

    # If requested re-normalize the distribution
    if normP:
        Pnorm = np.trapz(Pfit,r)
        Pfit = Pfit/Pnorm
        if uqanalysis:
            Puq_ = copy.deepcopy(Puq) # need a copy to avoid infite recursion on next step
            Puq.ci = lambda p: Puq_.ci(p)/Pnorm

    return Pfit, fitparam
    # =========================================================================
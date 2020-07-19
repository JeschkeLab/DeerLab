import numpy as np
from deerlab.nnls import fnnls, cvxnnls, nnlsbpp
from scipy.optimize import nnls
import deerlab as dl
import copy
from deerlab.utils import hccm
from deerlab.uqst import uqst

def fitregmodel(V,K,r, regtype='tikhonov', alpha='aic', regorder=2, solver='cvx', weights=1, huberparam=1.35, nonnegativity=True, obir = False, renormalize=True):

    V, K, weights = dl.utils.parse_multidatasets(V, K, weights)

    # Compute regularization matrix
    L = dl.regoperator(r,regorder)

    # Determine the type of problem to solve
    if obir:
        problem = 'obir'
    elif nonnegativity:
        problem = 'nnls'
    else:
        problem = 'unconstrained'

    # If the regularization parameter is not specified, get the optimal choice
    if type(alpha) is str:
        alpha = dl.selregparam(V,K,r,regtype,alpha,regorder=regorder,weights=weights, nonnegativity=nonnegativity,huberparam=huberparam)

    # Prepare components of the LSQ-problem
    [KtKreg, KtV] = dl.lsqcomponents(V,K,L,alpha,weights, regtype=regtype, huberparam=huberparam)

    # Unconstrained LSQ problem
    if problem == 'unconstrained':
        Pfit = np.linalg.solve(KtKreg,KtV)

    # Osher-Bregman iterated regularization
    elif problem == 'obir':
        Pfit = 0
    # Non-negative LSQ problem
    elif problem == 'nnls':

        if solver == 'fnnls':
            Pfit = fnnls(KtKreg,KtV)
        elif solver == 'nnlsbpp':
            Pfit = nnlsbpp(KtKreg,KtV,np.linalg.solve(KtKreg,KtV))
        elif solver == 'cvx':
            Pfit = cvxnnls(KtKreg, KtV)
        else:
            raise KeyError(f'{solver} is not a known non-negative least squares solver')

    # Uncertainty quantification
    # ----------------------------------------------------------------
    # Construct residual parts for for the residual and regularization terms
    res = weights*(V - K@Pfit)
    reg = alpha*L@Pfit

    # Construct Jacobians for the residual and penalty terms
    Jres = weights*K
    Jreg = alpha*L

    # Augment residual and Jacobian
    res = np.concatenate((res,reg))
    J = np.concatenate((Jres,Jreg))

    # Calculate the heteroscedasticity consistent covariance matrix 
    covmat = hccm(J,res,'HC1')
    
    # Construct confidence interval structure for P
    NonNegConst = np.zeros(len(r))
    Puq = uqst('covariance',Pfit,covmat,NonNegConst,[])

    # Re-normalization of the distributions
    # --------------------------------------
    if renormalize:
        Pnorm = np.trapz(Pfit,r)
        Pfit = Pfit/Pnorm
        Puq_ = copy.deepcopy(Puq) # need a copy to avoid infite recursion on next step
        Puq.ci = lambda p: Puq_.ci(p)/Pnorm


    return Pfit,Puq

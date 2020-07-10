import numpy as np
import math as m
import scipy.optimize as opt
from numpy.linalg import solve
from regoperator import regoperator
from selregparam import selregparam
from lsqcomponents import lsqcomponents
from fnnls import fnnls,cvxnnls
from jacobianest import jacobianest
from uqst import uqst
from cvxopt import matrix, solvers
import copy

def snlls(y,Amodel,par0,lb=[],ub=[],lbl=[],ubl=[],linsolver='fnnls'):

    getUncertainty = True
    getGoodnessOfFit = False 

    # Default optional settings
    alphaOptThreshold = 1e-3
    RegOrder = 2
    multiStarts = 1
    includePenalty = []
    weights = 1
    RegParam = 'aic'
    nonLinTolFun = 1e-5
    nonLinMaxIter = 1e4
    LinTolFun = 1e-5
    LinMaxIter = 1e4
    nonLinSolver = 'lsqnonlin'
    LinSolver = 'lsqlin'
    check = False
    regparam_prev = 0
    par_prev = [0]*len(par0)
    illConditioned = True
    linearConstrained = True
    nonNegativeOnly = True

    Nnonlin = len(par0)
    Nlin = np.shape(Amodel(par0))[1]
    linfit = []


    if not any(lb):
        lb = np.full(Nnonlin, -np.inf)

    if not any(ub):
        ub = np.full(Nnonlin, np.inf)

    if not any(lbl):
        ubl = np.full(Nlin, np.inf)

    if not any(ubl):
        ubl = np.full(Nlin, np.inf)

    includePenalty = True


    if includePenalty:
        # Use an arbitrary axis
        ax = list(range(1,Nlin+1))
        # Get regularization operator
        RegOrder = min(Nlin-1,RegOrder)
        L = regoperator(ax,RegOrder)
    else:
        L = np.eye(Nlin,Nlin)


    multiStartPar0 = par0

    # Pre-allocate containers for multi-start run
    fvals = [0]*multiStarts
    nonlinfits = []*multiStarts
    linfits = []*multiStarts


    # Residual vector function
    # ------------------------------------------------------------------
    # Function that provides vector of residuals, which is the objective
    # function for the least-squares solvers
    def ResidualsFcn(p):
        
        nonlocal par_prev, check, regparam_prev, linfit
        # Non-linear model evaluation
        # ===============================
        A = Amodel(p)
        
        # Regularization components
        # ===============================
        if includePenalty:
            if RegParam == 'aic':
                # If the parameter vector has not changed by much...
                if check and all(abs(par_prev-p)/p < alphaOptThreshold):
                    # ...use the alpha optimized in the previous iteration
                    alpha = regparam_prev
                else:
                    # ...otherwise optimize with current settings
                    alpha = selregparam(y,A,ax)
                    check = True
            else:
                # Fixed regularization parameter
                alpha =  RegParam
            
            # Store current iteration data for next one
            par_prev = p
            regparam_prev = alpha
            
            # Non-linear operator with penalty
            AtA,Aty = lsqcomponents(y,A,L,alpha)
            
        else:
            # Non-linear operator without penalty
            AtA,Aty = lsqcomponents(y,A,L,0)
        
        if not linearConstrained:
            # Unconstrained linear LSQ
            # ====================================
            linfit = solve(AtA,Aty)
            
        elif linearConstrained and not nonNegativeOnly:
            # Constrained linear LSQ
            # ====================================
            linfit = linSolverFcn(AtA,Aty,lbl,ubl,linSolverOpts)
            
        elif linearConstrained and nonNegativeOnly:
            # Non-negative linear LSQ
            # ====================================
            if linsolver=='fnnls':
                linfit = fnnls(AtA,Aty)
            elif linsolver=='nnls':
                try:
                    linfit,rnorm = opt.nnls(AtA,Aty,1e5)
                except:
                    linfit = np.zeros(np.shape(AtA)[1])
                    pass
            elif linsolver=='cvx':
                linfit =  cvxnnls(AtA, Aty, A, y)

        # Evaluate full model residual
        # ===============================
        yfit = A@linfit
        # Compute residual vector
        res = yfit - y
        if includePenalty:
            penalty = alpha*L@linfit
            # Augmented residual
            res = np.concatenate((res, penalty))

        return res
    # ------------------------------------------------------------------
        
    # Run the non-linear solver
    sol = opt.least_squares(ResidualsFcn,par0,bounds=[lb,ub])
    nonlinfit = sol.x


    # Uncertainty quantification
    # ------------------------------------------------------------------
    # Function that computes the covariance-based uncertainty quantification
    # and returns the corresponding uncertainty structure
    def uncertainty(nonlinfit):
            
        # Get full augmented residual vector
        #====================================
        res = ResidualsFcn(nonlinfit)
                
        # Compute the full augmented Jacobian
        # =====================================
        fcn = lambda p: Amodel(p)@linfit
        Jnonlin,err = jacobianest(fcn,nonlinfit)
        Jlin = Amodel(nonlinfit)
        Jreg = np.concatenate((np.zeros((np.shape(L)[0],np.shape(Jnonlin)[1])), regparam_prev*L),1)
        J = np.concatenate((Jnonlin, Jlin),1)
        J = np.concatenate((J, Jreg),0)
    
        # Estimate the heteroscedasticity-consistent covariance matrix
        sig = np.std(res)
        covmatrix = sig**2*np.linalg.pinv(J.T@J)
        
        parfit = np.concatenate((nonlinfit, linfit))
        lbs = np.concatenate((lb, lbl))
        ubs = np.concatenate((ub, ubl))
        paramuq_ = uqst('covariance',parfit,covmatrix,lbs,ubs)
        paramuq = copy.deepcopy(paramuq_)

        # Wrapper around the CI function handle of the uncertainty structure
        # ------------------------------------------------------------------
        def ci(coverage,type='full'):
            # Get requested confidence interval of joined parameter set
            paramci = paramuq_.ci(coverage)
            if type == 'nonlin':
                    # Return only confidence intervals on non-linear parameters
                    paramci = paramci[range(Nnonlin),:]
            elif type == 'lin':
                    # Return only confidence intervals on linear parameters
                    paramci = paramci[Nnonlin:,:]
            return paramci

        paramuq.ci = ci
        return paramuq
    # -------------------------------------------------------------------
        

    # Uncertainty analysis (if requested)
    if getUncertainty:
        uq = uncertainty(nonlinfit)

    # Goodness of fit (if requested)
    if getGoodnessOfFit:
            Ndof = len(y) - Nlin - Nnonlin
            stats = gof(y,yfit,Ndof)

    return nonlinfit, linfit, uq
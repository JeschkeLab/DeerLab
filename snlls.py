import numpy as np
import math as m
import scipy.optimize as opt
from numpy.linalg import solve
from regoperator import regoperator
from selregparam import selregparam
from lsqcomponents import lsqcomponents
from fnnls import fnnls

def snlls(y,Amodel,par0,lb=[],ub=[],ubl=[],lbl=[]):

    getUncertainty = False
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
    Nlin = np.shape(Amodel(par0))[1]
    linfit = []

    # Decide whether to include the regularization penalty
    if not includePenalty and illConditioned:
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
            linfit = fnnls(AtA,Aty)
        
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
    # Uncertainty analysis (if requested)
    if getUncertainty:
        paramuq = uncertainty(nonlinfit)

    # Goodness of fit (if requested)
    if getGoodnessOfFit:
            Ndof = len(y) - Nlin - Nnonlin
            stats = gof(y,yfit,Ndof)

    return nonlinfit, linfit
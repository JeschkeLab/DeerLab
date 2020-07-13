import numpy as np
import math as m
from scipy.optimize import least_squares, lsq_linear
from numpy.linalg import solve
from cvxopt import matrix, solvers
import copy
# Import DeerLab depencies
from deerlab.regoperator import regoperator
from deerlab.selregparam import selregparam
from deerlab.lsqcomponents import lsqcomponents
from deerlab.fnnls import fnnls,cvxnnls
from deerlab.jacobianest import jacobianest
from deerlab.uqst import uqst

def snlls(y,Amodel,par0,lb=[],ub=[],lbl=[],ubl=[],linsolver='cvx'):

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
    
    if len(lb) != Nnonlin or len(ub) != Nnonlin:
        raise TypeError('The lower/upper bounds of the non-linear problem must have #i elements',Nnonlin)
    if len(lbl) != Nlin or len(ubl) != Nlin:
        raise TypeError('The lower/upper bounds of the linear problem must have #i elements',Nlin)
    
    # Check if the nonlinear and linear problems are constrained
    nonLinearConstrained = (not np.all(np.isinf(lb))) or (not np.all(np.isinf(ub)))
    linearConstrained = (not np.all(np.isinf(lbl))) or (not np.all(np.isinf(ubl)))
    # Check for non-negativity constraints on the linear solution
    nonNegativeOnly = (np.all(lbl==0)) and (np.all(np.isinf(ubl)))
    
    # Check that the boundaries are valid
    if np.any(ub<lb) or np.any(ubl<lbl):
        raise ValueError('The upper bounds cannot be larger than the lower bounds.')
    # Check that the non-linear start values are inside the box constraint
    if np.any(par0>ub) or np.any(par0<lb):
        raise ValueError('The start values are outside of the specified bounds.')
        

    includePenalty = True


    if includePenalty:
        # Use an arbitrary axis
        ax = np.arange(1,Nlin+1)
        # Get regularization operator
        RegOrder = np.minimum(Nlin-1,RegOrder)
        L = regoperator(ax,RegOrder)
    else:
        L = np.eye(Nlin,Nlin)



    if not linearConstrained:
        # Unconstrained linear LSQ
        linSolver = lambda AtA,Aty: solve(AtA,Aty)
    elif linearConstrained and not nonNegativeOnly:
        # Constrained linear LSQ
        linSolver = lambda AtA,Aty: lsq_linear(AtA, Aty, bounds=[lbl,ubl], tol=LinTolFun, max_iter=LinMaxIter)
    elif linearConstrained and nonNegativeOnly and linsolver=='fnnls':
        # Non-negative linear LSQ
        linSolver = lambda AtA,Aty: fnnls(AtA,Aty, tol=LinTolFun, maxiter=LinMaxIter)
    elif linearConstrained and nonNegativeOnly and linsolver=='cvx':
        # Non-negative linear LSQ
        linSolver = lambda AtA,Aty: cvxnnls(AtA, Aty, tol=LinTolFun, maxiter=LinMaxIter)

    multiStartPar0 = par0

    # Pre-allocate containers for multi-start run
    fvals = [0]*multiStarts
    nonlinfits = []*multiStarts
    linfits = []*multiStarts



    def ResidualsFcn(p):
    #===========================================================================
        """ 
        Residuals function
        ------------------
        Provides vector of residuals, which is the objective function for the non-linear least-squares solver. 
        """

        nonlocal par_prev, check, regparam_prev, linfit
        # Non-linear model evaluation
        A = Amodel(p)
        
        # Regularization components
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
            
        else:
            # Non-linear operator without penalty
            alpha = 0

        # Components for linear least-squares        
        AtA,Aty = lsqcomponents(y,A,L,alpha)

        # Solve the linear least-squares problem
        linfit = linSolver(AtA,Aty)

        # Evaluate full model residual
        yfit = A@linfit
        # Compute residual vector
        res = yfit - y
        if includePenalty:
            penalty = alpha*L@linfit
            # Augmented residual
            res = np.concatenate((res, penalty))

        return res
    #===========================================================================
        
    # Run the non-linear solver
    sol = least_squares(ResidualsFcn,par0,bounds=[lb,ub])
    nonlinfit = sol.x

    def uncertainty(nonlinfit):
    #===========================================================================
        """ 
        Uncertainty quantification
        --------------------------
        Function that computes the covariance-based uncertainty quantification and returns the corresponding uncertainty structure
        """

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

        def ci(coverage,type='full'):
        #===========================================================================
            "Wrapper around the CI function handle of the uncertainty structure"
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
        #===========================================================================
    #===========================================================================
        
    # Uncertainty analysis (if requested)
    if getUncertainty:
        uq = uncertainty(nonlinfit)

    # Goodness of fit (if requested)
    if getGoodnessOfFit:
            Ndof = len(y) - Nlin - Nnonlin
            stats = gof(y,yfit,Ndof)

    return nonlinfit, linfit, uq


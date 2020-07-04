import numpy as np 
import scipy.optimize as opt
import math as m
from numpy.linalg import norm
from lsqcomponents import lsqcomponents
from regoperator import regoperator
from fnnls import fnnls

def selregparam(V,K,r,RegType='tikhonov',SelectionMethod='aic'):

    LcurveMethods = SelectionMethod in ['lr','lc']
    RegOrder = 2
    TolFun = 1e-9
    NoiseLevel = 0
    NonNegConstrained = True

    if LcurveMethods:
        SearchMethod = 'grid'
    else:
        SearchMethod = 'fminbnd'


    #  Preparations
    #-------------------------------------------------------------------------------
    # Get regularization operator
    L = regoperator(r,RegOrder)
    # Get range of potential alpha values candidates
    alphaRange = np.logspace(-4,4,60)

#--------------------------------------------------------------------------
    def evalalpha(alpha,SelectionMethod):
                
        #-----------------------------------------------------------------------
        #  Pseudo-Inverses and Ps
        #-----------------------------------------------------------------------
        
        KtKreg, KtV = lsqcomponents(V,K,L,alpha)
        if NonNegConstrained:
            P = fnnls(KtKreg,KtV,TolFun)
        else:
            P = np.linalg.solve(KtKreg,KtV)

        PseudoInverse = np.linalg.inv(KtKreg)@K.T
        
        Residual = norm(K@P - V)
        InfluenceMatrix = K@PseudoInverse
        
        #-----------------------------------------------------------------------
        #  Selection methods for optimal regularization parameter
        #-----------------------------------------------------------------------
        
        # If multiple selection methods are requested then process them sequentially
        Functional = 0
        nr = len(V)
        if  SelectionMethod =='cv': # Cross validation (CV)
            InfluenceDiagonal = np.diag(InfluenceMatrix)
            f_ = sum(abs(V - K*(P)/(np.ones(nr,1) - InfluenceDiagonal))**2)
                
        elif  SelectionMethod =='gcv': # Generalized Cross Validation (GCV)
            f_ = Residual**2/((1 - np.trace(InfluenceMatrix)/nr)**2)
                
        elif  SelectionMethod =='rgcv': # Robust Generalized Cross Validation (rGCV)
            TuningParameter = 0.9
            f_ = Residual**2/((1 - np.trace(InfluenceMatrix)/nr)**2)*(TuningParameter + (1 - TuningParameter)*np.trace(InfluenceMatrix**2)/nr)
                
        elif  SelectionMethod =='srgcv': # Strong Robust Generalized Cross Validation (srGCV)
            TuningParameter = 0.8
            f_ = Residual**2/((1 - np.trace(InfluenceMatrix)/nr)**2)*(TuningParameter + (1 - TuningParameter)*np.trace(PseudoInverse.T*PseudoInverse)/nr)
                
        elif  SelectionMethod =='aic': # Akaike information criterion (AIC)
            Criterion = 2
            f_ = nr*np.log(Residual**2/nr) + Criterion*np.trace(InfluenceMatrix)
                
        elif  SelectionMethod =='bic':  # Bayesian information criterion (BIC)
            Criterion = np.log(nr)
            f_ = nr*np.log(Residual**2/nr) + Criterion*np.trace(InfluenceMatrix)
                
        elif  SelectionMethod =='aicc': # Corrected Akaike information criterion (AICC)
            Criterion = 2*nr/(nr-np.trace(InfluenceMatrix)-1)
            f_ = nr*np.log(Residual**2/nr) + Criterion*np.trace(InfluenceMatrix)
                
        elif  SelectionMethod =='rm': # Residual method (RM)
            Scaling = K.T@(np.eye(np.shape(InfluenceMatrix)) - InfluenceMatrix)
            f_ = Residual**2/np.sqrt(np.trace(Scaling.T@Scaling))
                
        elif  SelectionMethod =='ee': # Extrapolated Error (EE)
            f_ = Residual**2/norm(K.T*(K@P - V))
                
        elif  SelectionMethod =='gml': # Generalized Maximum Likelihood (GML)
            Treshold = 1e-9
            EigenValues = np.linalg.eig(np.eye(np.shape(InfluenceMatrix)) - InfluenceMatrix)
            EigenValues[EigenValues < Treshold] = 0
            NonZeroEigenvalues = np.real(EigenValues[EigenValues!=0])
            f_ = V.T*(V - K@P)/m.exp(sum(map(np.log,NonZeroEigenvalues)))*(1/len(NonZeroEigenvalues))
                
        elif  SelectionMethod =='mcl':  # Mallows' C_L (MCL)
            f_ = Residual**2 + 2*NoiseLevel**2*np.trace(InfluenceMatrix) - 2*nr*NoiseLevel**2
            
        else:
            f_ = 0
        Functional = Functional + f_

        return Functional       


    # Evaluate functional over search range, using specified search method
    #-------------------------------------------------------------------------------
    if SearchMethod == 'fminbnd':
                    
        lga_min = m.log10(min(alphaRange))
        lga_max = m.log10(max(alphaRange))

        fun = lambda lga: evalalpha(10**lga,SelectionMethod)
        # Optimize alpha via the fminbnd function
        lga_opt = opt.fminbound(fun,lga_min,lga_max)
        alphaOpt = 10**lga_opt

    return alphaOpt

         



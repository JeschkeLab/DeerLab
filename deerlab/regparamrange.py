import numpy as np
from scipy.linalg import svd
from deerlab.utils import gsvd

def regparamrange(K,L,noiselvl=0,logres=0.1):
#===================================================================================================
    r"""
    Estimates an array of regularization parameter candidates from
    the generalized singular value decomposition (GSVD) of the experiment
    dipolar kernel  and regularization operator.
    
    Parameters
    ----------
    K : array_like with shape(m,n)  
        Dipolar kernel matrix
    L : array_like with shape(n-2,n)  
        Regularization operator
    noiselvl : float scalar, optional
        Estiamte of the noise standard deviation in the data. Used for scaling 
        of the singular values, the default is 0.
    logres : float scalar, optional
        Resolution of the array of alpha candidates, on a base-10 logarithmic scale, the default is 0.1 
        (corresponding to a tenth of a decade).

    Returns
    -------
    alphas : ndarray
        Regularization parameter candidates

    """
    # Set alpha range
    minmax_ratio = 16*np.finfo(float).eps*1e6  # ratio of smallest to largest alpha

    # Scaling by noise. This improves L curve corner detection for DEER.
    minmax_ratio = minmax_ratio*2**(noiselvl/0.0025)

    # Get generalized singular values of K and L
    singularValues = gsvd(K,L)

    DerivativeOrder = L.shape[1] - L.shape[0] # get order of derivative (=number of inf in singval)
    singularValues = singularValues[0:len(singularValues)-DerivativeOrder] # remove inf 
    singularValues = singularValues[::-1] # sort in decreasing order
    singularValues = singularValues[singularValues>0] # remove zeros
    lgsingularValues = np.log10(singularValues)

    # Calculate range based on singular values
    lgrangeMax = lgsingularValues[0]
    lgrangeMin = np.maximum(lgsingularValues[-1],lgsingularValues[0]+np.log10(minmax_ratio))
    lgrangeMax = np.floor(lgrangeMax/logres)*logres
    lgrangeMin = np.ceil(lgrangeMin/logres)*logres
    if lgrangeMax < lgrangeMin:
        temp = lgrangeMax
        lgrangeMax = lgrangeMin
        lgrangeMin = temp
    lgalpha = np.arange(lgrangeMin,lgrangeMax,logres)
    lgalpha = np.append(lgalpha,lgrangeMax)
    alphas = 10**lgalpha

    return alphas
#===================================================================================================

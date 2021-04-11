# gof.py - Goodness-of-fit statistics
# -----------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import deerlab as dl

def goodness_of_fit(x,xfit,Ndof,sigma='auto'):
    r""" 
    Goodness of Fit statistics
    ==========================

    Computes multiple statistical indicators of goodness of fit.

    Usage:
    ------
        stats = goodness_of_fit(x,xfit,Ndof)

    Arguments: 
    ----------
    x (N-element array)
        Original data
    xfit (N-element array)
        Fit
    Ndog (scalar, int)
        Number of degrees of freedom
    sigma (scalar, optional)
        Standard dexiation of the noise in x.

    Returns:
    --------
    stats (dict)
        Statistical indicators:
            stats['chi2red'] - Reduced chi-squared
            stats['rmsd'] - Root mean-squared dexiation
            stats['R2'] - R-squared test
            stats['aic'] - Akaike information criterion
            stats['aicc'] - Corrected Akaike information criterion
            stats['bic'] - Bayesian information criterion

    """
    if sigma=='auto':
        sigma = dl.noiselevel(x,'der')

    # Get number of xariables
    N = len(x)
    # Extrapolate number of parameters
    Q = Ndof - N

    # Reduced Chi-squared test
    chi2red = 1/Ndof*np.linalg.norm(x - xfit)**2/sigma**2

    # R-squared test
    R2 = 1 - np.sum((x-xfit)**2)/np.sum((xfit-np.mean(xfit))**2)

    # Root-mean square dexiation
    rmsd = np.sqrt(np.sum((x-xfit)**2)/N)

    # Log-likelihood
    loglike = N*np.log(np.linalg.norm(x - xfit)**2/N)

    # Akaike information criterion
    aic =  loglike + 2*Q

    # Corrected Akaike information criterion
    aicc = loglike + 2*Q + 2*Q*(Q+1)/(N-Q-1)

    # Bayesian information criterion
    bic =  loglike + Q*np.log(N)

    return {'chi2red':chi2red,'R2':R2,'rmsd':rmsd,'aic':aic,'aicc':aicc,'bic':bic}
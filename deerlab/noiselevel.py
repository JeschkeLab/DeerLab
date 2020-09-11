# noiselevel.py - Noise level estimator
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

from numpy import isreal, std, mean, shape, atleast_1d
from deerlab.utils import movmean
from deerlab import correctphase
from scipy.signal import savgol_filter
import warnings

def noiselevel(V,*args):
    r"""
    Returns the standard deviation estimation of the noise in a given signal using different methods:


    * ``sigma = noiselevel(V2D)``: If ``V`` is a 2D-dataset of different scans, the noise standard deviation
      is estimated from the deviations between scans. The second dimension of
      ``V2D`` must contain the different scans. The function returns the standard
      deviation of the averaged signal not of the individual scans

    * ``sigma = noiselevel(V,filter)``: If a 1D signal ``V`` is given, the noise level is estimated via filtering
      of the signal with a moving mean filter. The nature of the filter can
      be specified by means of a string ``filter``.

    * ``sigma = noiselevel(V,Vref)``: If a reference model signal ``Vref`` is given, the noise level is
      estimated from the difference between both signals.

    * ``sigma = noiselevel(Vco)``: If the input signal ``Vco`` contains an imaginary component, the noise
      level is estimated form the imaginary component after phase optimization.

    Parameters
    ----------
    V2D : 2D-array_like
        Datasets of single scans of a dipolar signal.
    V : array_like
        Real-valued dipolar signal
    filter : string
        Filtering method:

        * ``'movmean'`` - Moving mean filter
        * ``'savgol'``  - Savitzky-Golay filter
        The default is ``'movmean'``
        
    Vref : array_like
        Reference dipolar signal.
    Vco : array_like
        Complex-valued dipolar signal.
        
    Returns
    -------
    sigma : scalar
        Estimated noise standard deviation
    
    """

    # Parse the multiple input schemes
    # --------------------------------

    V = atleast_1d(V)

    # Input: noiselevel(V2D)
    if V.ndim == 2:
        estimationMethod = '2D'
        if args:
            raise KeyError('For 2D-datasets, only one input is required.')

    # Input: noiselevel(V)
    elif V.ndim==1 and not args and all(isreal(V)):
        estimationMethod = 'filtering'
        filterType = 'movmean'

    # Input: noiselevel(Vco)
    elif V.ndim==1 and not args and not all(isreal(V)):
        estimationMethod = 'complex'
        
    # Input: noiselevel(V,filter)
    elif V.ndim==1 and type(args[0]) is str and all(isreal(V)):
        estimationMethod = 'filtering'
        filterType = args[0]

    # Input: noiselevel(V,Vref)
    elif V.ndim==1 and type(args[0]) is not str and all(isreal(V)):
        estimationMethod = 'reference'
        Vref = args[0]
        if len(V) != len(Vref):
            raise TypeError('The input and reference signal must have the same number of elements.') 
    else:
        raise KeyError('The input is not valid.')


    # Estimation of the noise level
    # -----------------------------

    if estimationMethod == '2D':
            # Estimate standard deviations for all time point, and average over scans
            if shape(V)[1] < 10:
                raise Warning('Only a few scans are given. Noise standard deviation estimate will be inaccurate.')
            sigma = std(V,1)
            sigma = mean(sigma)
            
    elif estimationMethod == 'filtering':
        # Filter the noise in the signal    
        if filterType == 'movmean':
                Vfilt = movmean(V,3)
        elif filterType == 'savgol':
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                Vfilt = savgol_filter(V,11,3)
        else:
            raise TypeError("Filter type not found. Must be either 'savgol' or 'movmean'.")
        # And estimate the noiselevel from the resulting residual
        sigma = std(V - Vfilt)
            
    elif estimationMethod == 'complex':
            # Optimize the phase of the signal
            _,Vim,_,_ = correctphase(V,full_output=True)
            # And estimate the noiselevel from the imaginary part
            sigma = std(Vim)
            
    elif estimationMethod == 'reference':
            sigma = std(V - Vref)

    return sigma

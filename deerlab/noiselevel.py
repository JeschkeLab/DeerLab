# noiselevel.py - Noise level estimator
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

from numpy import isreal, std, mean, shape, atleast_1d
from deerlab.utils import movmean, der_snr
from deerlab import correctphase
from scipy.signal import savgol_filter
import warnings
import numpy as np
def noiselevel(V,mode='der',*args):
    r"""
    Returns the standard deviation estimation of the noise in a given signal using different methods:

    * ``sigma = noiselevel(V,'der')``: Employs the DER_SNR [1] method for estimating the noise standard deviation.

    * ``sigma = noiselevel(V,'scans')``: If ``V`` is a 2D-dataset of different scans, the noise standard deviation
      is estimated from the deviations between scans. The second dimension of
      ``V2D`` must contain the different scans. The function returns the standard
      deviation of the averaged signal not of the individual scans.

    * ``sigma = noiselevel(V,'movmean',winsize)``: The noise level is estimated via filtering
      of the signal with a moving mean filter. The size ``winsize`` of the moving mean window must be specified
      as well.

    * ``sigma = noiselevel(V,'movmean',winsize,order)``: The noise level is estimated via filtering
      of the signal with a Savitzky-Golay filter. The  window size ``winsize`` and polynomial order ``order`` 
      must be specified as well.

    * ``sigma = noiselevel(V,'reference',Vref)``: If a reference model signal ``Vref`` is given, the noise level is
      estimated from the difference between both signals.

    * ``sigma = noiselevel(V,'complex')``: If the input signal ``V`` contains an imaginary component, the noise
      level is estimated form the imaginary component after phase optimization.

    Parameters
    ----------
    V : array_like
        Input signal. Must be a real-valued 2D-matrix for the ``'scans'`` method. Otherwise it must be a real-valued 1D-array.
    method : string 

        * ``'der'`` - Estimation via the DER_SNR method [1].
        * ``'scans'`` - Estimation from deviations between scans. 
        * ``'movmean'`` - Estimation from a moving mean filtered signal.
        * ``'savgol'``  - Estimation from a Savitzky-Golay filtered signal.
        * ``'reference'`` - Estimation from comparison with a reference signal. 
        * ``'complex'`` - Estimation from phase-corrected imaginary part of a complex signal.

        The default is ``'der'``. 

    winsize : scalar integer
        Filter window size
    order : scalar integer
        Savitzky-Golay polynomial order.
    Vref : array_like
        Reference dipolar signal.
    Vco : array_like
        Complex-valued dipolar signal.
        
    Returns
    -------
    sigma : scalar
        Estimated noise standard deviation
    
    References
    ----------
    .. [1] S. Czesla, T. Molle and J. H. M. M. Schmitt
        A posteriori noise estimation in variable data sets - With applications to spectra and light curves, A&A, 609 (2018) A39 

    """

    V = atleast_1d(V)

    if mode == 'scans':
        if args:
            raise KeyError("For the 'der' method, no additional inputs are required.")

        if V.ndim!=2:
            raise KeyError("For the 'scans' method, the input signal must be a 2D-matrix.")
        # Estimate standard deviations for all time point, and average over scans
        if shape(V)[1] < 10:
            raise Warning('Only a few scans are given. Noise standard deviation estimate will be inaccurate.')
        sigma = std(V,1)
        sigma = mean(sigma)
            
    elif mode == 'movmean':
        # Filter the noise in the signal    
        if args:
            if len(args)!=1:
                raise KeyError("For the 'movmean' method, only the windows size can be specified.")
            winsize = args
        else:
            raise KeyError("For the 'movmean' method, the windows size must be specified.")
        Vfilt = movmean(V,winsize)
        sigma = std(V - Vfilt)

    elif mode == 'savgol':
        # Filter the noise in the signal    
        if args:
            if len(args)!=2:
                raise KeyError("For the 'savgol' method, only the size and order can be specified.")
            winsize,order = args
        else:
                raise KeyError("For the 'savgol' method, the size and order must be specified.")
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            Vfilt = savgol_filter(V,winsize,order)
        sigma = std(V - Vfilt)

    elif mode == 'complex':
        if all(isreal(V)):
            raise TypeError("For the 'complex' method, the input signal must be complex-valued.")
        # Optimize the phase of the signal
        _,Vim,_ = correctphase(V,full_output=True)
        # And estimate the noiselevel from the imaginary part
        sigma = std(Vim)
            
    elif mode == 'reference':
        if args: 
            Vref = args[0]
        else:
            raise KeyError("For the 'reference' method, a reference signal must be provided.")
        if len(V) != len(Vref):
            raise TypeError('The input and reference signal must have the same number of elements.') 
        sigma = std(V - Vref)

    elif mode == 'der':
        if args:
            raise KeyError("For the 'der' method, no additional inputs are required.")

        sigma = der_snr(V)

    return sigma

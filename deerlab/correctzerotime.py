
# correctzerotime.py - Zero-time correction of dipolar signals
#--------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import scipy as scp

# ==============================================================================
def correctzerotime(V,t):
    r"""
    Determines and corrects the zero-time of a dipolar signal via a first-moment integral analysis [1]_.
    
    Parameters
    ----------
    V : array_like
        Dipolar signal.
    t : array_like
        Time axis, in microseconds.

    Returns
    --------
    tcorr : ndarray
        Zero-time corrected time axis, in microseconds.

    References
    ----------
    .. [1] G. Jeschke, Appl.Magn.Reson. 30, 473-498 (2006) p.485
    """

    if len(t)!=len(V):
        raise TypeError('The dipolar signal (V) and time axis (t) must have the size.')

    # Rescale to avoid overflow during integration later
    V /= max(V)
    # Generate finely-grained interpolated signal and time axis
    resolution = 10
    t_ = np.linspace(min(t), max(t), (len(t)-1)*resolution+1)
    V_ = scp.interpolate.interp1d(t, np.real(V), kind='cubic')(t_)

    # Determine location of maximum
    idxmax = np.argmax(V_)

    if idxmax!=1 and idxmax!=len(V):
        # If maximum is not the first or last point, then do moment analysis
        # (i.e. minimize the magntitude of the first-moment integral)
        
        # Determine the interval to use for the integral
        maxDelta = np.floor(np.minimum(idxmax, len(V_)-idxmax)/2)
        offsetrange = np.arange(-maxDelta,maxDelta+1, dtype=np.int32)

        # Look around the maximum for lowest-magnitude moment integral
        minIntegral = np.inf
        absInt = lambda i: abs(np.sum(V_[offsetrange+i]*(t_[offsetrange+i]-t_[i])))
        # Run over all zero-position candidates and evalulate the integral
        for idx in np.arange(idxmax-maxDelta, idxmax+maxDelta, dtype=np.int32):
            absIntegral = absInt(idx)
            # If |integral| is smaller than prior best, then update candidate
            if absIntegral < minIntegral:
                minIntegral = absIntegral
                idxt0 = idx
    else:
        # Return maximum if it's the first or the last point
        idxt0 = idxmax

    t0 = t_[idxt0]

    # Correct time axis
    tcorr = t - t0

    # Determine index of time axis point closest to t0
    idxt0 = np.argmin(abs(tcorr))

    return tcorr
# ==============================================================================
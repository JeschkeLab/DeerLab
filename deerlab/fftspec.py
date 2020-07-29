
# fftspec.py - Fast-Fourier transform spectrum
# ---------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from numpy.fft import fft, fftshift

def fftspec(V,t,mode='abs',zerofilling='auto',apodization=True):
    """
    Fast-Fourier transform spectrum
    ===============================

    Computes the FFT spectrum of a signal.

    Usage:
    -------
        nu,spec = fftspec(t,S)
    
    Arguments:
    ----------
    V (array)
        Signal to be processed.
    t (array)
        Time axis, in microseconds.

    Returns:
    --------
    nu (array)
        Frequency axis, in megahertz
    spec (array)
        FFT spectrum.

    Additional keyword arguments:
    -----------------------------
    mode (string, default='abs') 
        Type of spectrum to be returned ('real','imag','abs')
    zerofilling (scalar, default = 2*len(V) )
        Number of elements in the output FFT spectrum
    aposization (boolean, default=True)
        Use of a Hamming apodization window
    """
    
    if zerofilling is 'auto':
        zerofilling = 2*len(V)

    #If requested apply Hamming apodization window
    if apodization:
        arg = np.linspace(0,np.pi,len(V))
        ApoWindow = 0.54 + 0.46*np.cos(arg)
        V = V*ApoWindow

    #Compute fft spectrum
    spec = fftshift(fft(V,zerofilling))

    #Get the requested component/type of spectrum
    if mode is 'abs':
            spec = np.abs(spec)
    elif mode is 'real':
            spec  = np.real(spec)
    elif mode is 'imag':
            spec = np.imag(spec)
    else:
        raise KeyError("Invalid spectrum mode. Must be 'abs', 'real', or 'imag'. ")

    #Get frequency axis and switch output order
    dt = np.mean(np.diff(t))
    nu = np.linspace(-1/(2*dt),1/(2*dt),zerofilling)

    return nu, spec 

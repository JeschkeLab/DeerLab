# correctscale.py - Scale correction pre-processing  of dipolar signals
# ---------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab import fitparamodel, dipolarkernel
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp
from deerlab.utils import isempty

def correctscale(V,t,tmax=[],model='deer'):
    r""" 
Amplitude scale correction

Takes the experimental dipolar signal (V) on a given time axis (t) and
rescales it so that it is 1 at time zero, taking into account the noise.
It fits a model specified in (model) over the interval ``abs(t)<tmax`` around time zero.
It returns the rescaled signal (Vc) and the scaling factor (V0).

Parameters
----------
V : array_like
    Experimenal signal.
    
t : array_like 
    Time axis, in microseconds. 

tmax : scalar, optional
    Time cutoff for fit range, in microseconds.

model : string, optional
    Model to fit over ``abs(t)<tmax``
    * ``'deer'`` - DEER model with Gaussian distribution.
    * ``'gauss'`` - Raised Gaussian.
    
    The default is ``'deer'``

Returns
-------
Vc : ndarray
    Rescaled signal.
    """
    if not all(np.isreal(V)):
        raise ValueError('Input signal cannot be complex.') 

    if isempty(tmax):
        tmax = max(t)
    # Approximate scaling
    Amp0 = max(V)
    V_ = V/Amp0

    # Limit fit to region around time zero
    idx = abs(t)<=tmax
    V_ = V_[idx]
    t_ = t[ idx]

    if model == 'deer':
        # DEER model with Gaussian distribution
        if len(V_)<5:
            raise KeyError('Number of points in fit range cannot be smaller than number of model parameters. Increase tmax.')
        r = np.linspace(0.5,20,300)
        fitmodel = lambda p: p[0]*(dipolarkernel(t_,r,mod=p[1],bg=bg_exp(t_,p[4]))@dd_gauss(r,p[[2,3]]))
        # parameters: amplitude, mod depth, Gauss position, Gauss std, conc
        par0 = [1, 0.5, 3, 0.3, 0.2]
        lb = [1e-3, 0.01, 0, 0.01, 0]
        ub = [10, 1, 20, 5, 100]
    elif model == 'gauss':
        # Gaussian function
        if len(V_)<3:
            raise KeyError('Number of points in fit range cannot be smaller than number of model parameters. Increase tmax.')
        fitmodel = lambda p: (p[0]-p[1])+p[1]*np.exp(-t_**2/p[2]**2)
        par0 = [1, 1, 1]
        lb = [1e-3, 1e-3, 1e-3]
        ub = [10, 1000, 10000]
    else: 
        raise KeyError(f"Unknown model '{model}'")

    # Run the parametric model fitting
    fit = fitparamodel(V_,fitmodel,par0,lb,ub,fitscale=False,uq=False)

    # Get the fitted signal amplitude and scale the signal
    V0 = Amp0*fit.param[0]
    Vc = V/V0

    return Vc
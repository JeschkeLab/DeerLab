
# dipolarbackground.py - Multipathway background generator 
# --------------------------------------------------------
 # This file is a part of DeerLab. License is MIT (see LICENSE.md).
 # Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types

def dipolarbackground(t,pathinfo,Bmodel,renormalize = True,overtonecoeff = 1):
    r""" Constructs background decay functions according to the multi-pathway model.
    
    Parameters
    ----------
    t : array_like
        Time axis, in microseconds.
    pathinfo : array_like with shape (p,2) or (p,3)
        Array of pathway amplitudes (lambda), refocusing points (T0), and harmonics (n) 
        for multiple (p) dipolar pathways. Each row contains ``[lambda T0 n]`` or ``[lambda T0]`` 
        for one pathway. If n is not given it is assumed to be 1. For a pathway with unmodulated contribution, set the refocusing time to ``numpy.nan``.
    Bmodel : callable
        Background basis function. A callable function accepting a time-axis array as first input and a pathway amplitude as a second, i.e. ``B = lambda t,lam: bg_model(t,par,lam)``
     
    Returns
    -------
    B : ndarray
        Multipathway background decay

    Other parameters
    ----------------
    renormalize : boolean
        The multi-pathway background does not necessarily satisfy ``V(0) == 1``. This option enables/disables a re-normalization of the background function to ensure that equality is satisfied, by default it is enabled.
    overtonecoeff : array_like
        Array containing the overtone coefficients for RIDME experiments. 

    Notes
    -----
    Computes the multipathway background ``B`` for the time axis ``t`` corresponding to the dipolar pathways specified in ``pathinfo``. The total background is computed from the basis background function model specified in ``Bmodel``. For a multi-pathway DEER signal (e.g, 4-pulse DEER with 2+1 contribution; 5-pulse DEER with 4-pulse DEER residual signal, and more complicated experiments), ``pathinfo`` is a 2D-array that contains a list of modulation depths (amplitudes), refocusing times (in microseconds), and optional harmonics for all modulated pathway signals

    Examples
    --------
    To calculate the background for the standard model for 4-pulse DEER with an unmodulated offset and a single dipolar pathway that refocuses at time 0, use::

        from numpy import linspace
        from deerlab import dipolarbackground

        t = linspace(-5,20,501) # time axis (us)
        lam = 0.4 # modulation depth main signal
        conc = 200   # spin concentration (uM)

        pathinfo = [[],[]]
        pathinfo[0] = [1-lam, NaN] # unmodulated part, gives offset
        pathinfo[1] = [lam, 0] # main modulation, refocusing at time zero

        Bmodel = @(t,lam) bg_hom3d(t,conc,lam)

        B = dipolarbackground(t,pathinfo,Bmodel)


    """

    # Ensure that all inputs are numpy arrays
    t = np.atleast_1d(t)
    pathinfo = np.atleast_1d(pathinfo)
    overtonecoeff = np.atleast_1d(overtonecoeff)


    if type(Bmodel) is not types.LambdaType:
        raise TypeError('For a model with multiple modulated pathways, B must be a function handle of the type: @(t,lambda) bg_model(t,par,lambda)')

    if not np.isreal(pathinfo).all:
        raise TypeError('lambda/pathinfo must be a numeric array.')

    if len(pathinfo) == 1:
        lam = pathinfo[0]
        pathinfo = np.array([[1-lam, np.NaN], [lam, 0]])

    if not np.any(np.shape(pathinfo)[1] != np.array([2, 3])):
        raise TypeError('pathinfo must be a numeric array with two or three columns.')

    if np.any(np.isnan(pathinfo[:, 0])):
        raise ValueError('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path[1,:] = [Lam0 NaN]')

    # Normalize the pathway amplitudes to unity
    pathinfo[:, 0] = pathinfo[:, 0]/sum(pathinfo[:, 0])
    lam = pathinfo[:, 0]
    T0 = pathinfo[:, 1]
    if np.shape(pathinfo)[1] == 2:
        n = np.ones(np.shape(T0))
    else:
        n = pathinfo[:, 2]
    

    # Combine all unmodulated components, and eliminate from list
    unmodulated = np.where(np.isnan(T0))
    lam = np.delete(lam, unmodulated)
    T0 = np.delete(T0, unmodulated)
    n = np.delete(n, unmodulated)

    # Fold overtones into pathway list
    nCoeffs = len(overtonecoeff)
    lam_,T0_,n_ = (np.empty(0) for _ in range(3))
    for i in range(nCoeffs):
        lam_ = np.concatenate((lam_, lam*overtonecoeff[i]))
        T0_ = np.concatenate((T0_,T0))
        n_ = np.concatenate((n_,n*(i+1)))
    lam,T0,n = (lam_,T0_,n_)
    nModPathways = len(lam)

    # Construction of multi-pathway background function 
    #-------------------------------------------------------------------------------
    Bnorm = 1
    B = 1
    for pathway in range(nModPathways):        
        B = B*Bmodel(n[pathway]*(t-T0[pathway]),lam[pathway])
        Bnorm = Bnorm*Bmodel(-T0[pathway]*n[pathway],lam[pathway])
    
    if renormalize:
        B = B/Bnorm

    return B
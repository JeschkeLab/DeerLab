# dipolarbackground.py - Multipathway background generator 
# --------------------------------------------------------
 # This file is a part of DeerLab. License is MIT (see LICENSE.md).
 # Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types

def dipolarbackground(t, pathinfo, Bmodel, renormalize=True):
    r""" Constructs background decay functions according to the multi-pathway model.
    
    Parameters
    ----------
    t : array_like
        Time axis, in microseconds.
    pathinfo : list of lists or scalar
        List of pathways. Each pathway is defined as a list of the pathway's amplitude (lambda), refocusing time (T0), 
        and harmonic (n), i.e. ``[lambda, T0, n]`` or ``[lambda, T0]`` for one pathway. If n is not given it is assumed to be 1. 
        For a pathway with unmodulated contribution, only the amplitude must be specified, i.e. ``[Lambda0]``.
        If a single value is specified, it is interpreted as the 4-pulse DEER pathway amplitude (modulation depth).  
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

    Notes
    -----
    Computes the multipathway background ``B`` for the time axis ``t`` corresponding to the dipolar pathways specified in ``pathinfo``. The total background is computed from the basis background function model specified in ``Bmodel``. For a multi-pathway DEER signal (e.g, 4-pulse DEER with 2+1 contribution; 5-pulse DEER with 4-pulse DEER residual signal, and more complicated experiments), ``pathinfo`` is a 2D-array that contains a list of modulation depths (amplitudes), refocusing times (in microseconds), and optional harmonics for all modulated pathway signals

    Examples
    --------
    To calculate the background for the standard model for 4-pulse DEER with an unmodulated offset and a single dipolar pathway that refocuses at time 0, use::

        import numpy as np
        import deerlab as dl

        t = np.linspace(-5,20,501) # time axis (us)
        lam = 0.4 # modulation depth main signal
        conc = 200   # spin concentration (uM)

        pathinfo = []
        pathinfo.append([1-lam]) # unmodulated part, gives offset
        pathinfo.append([lam, 0]) # main modulation, refocusing at time zero

        def Bmodel(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t,pathinfo,Bmodel)
    """

    # Ensure that all inputs are numpy arrays
    t = np.atleast_1d(t)

    if type(Bmodel) is not types.LambdaType:
        raise TypeError('For a model with multiple modulated pathways, B must be a function handle of the type: @(t,lambda) bg_model(t,par,lambda)')

    if not isinstance(pathinfo,list): pathinfo = [pathinfo] 
    if len(pathinfo) == 1:
        lam = pathinfo[0]
        pathinfo = [[1-lam], [lam, 0]]

    pathinfo = [np.atleast_1d(path) for path in pathinfo]
    
    # Get unmodulated pathways    
    unmodulated = [pathinfo.pop(i) for i,path in enumerate(pathinfo) if len(path)==1]
    # Combine all unmodulated contributions
    Lambda0 = sum(np.concatenate([path for path in unmodulated]))

    # Check structure of pathways
    for i,path in enumerate(pathinfo):
        if len(path) == 2:
            # If harmonic is not defined, append default n=1
            pathinfo[i] = np.append(path,1) 
        elif len(path) != 3:
            # Otherwise paths are not correctly defined
            raise KeyError('The pathway #{} must be a list of two or three elements [lam, T0] or [lam, T0, n]'.format(i))

    # Normalize the pathway amplitudes to unity
    lamsum = Lambda0 + sum([path[0] for path in pathinfo])
    Lambda0 /= lamsum
    for i in range(len(pathinfo)):
        pathinfo[i][0] /= lamsum 

    # Construction of multi-pathway background function 
    #-------------------------------------------------------------------------------
    Bnorm = 1
    B = 1
    for pathway in pathinfo:        
        n,T0,lam = pathway
        B = B*Bmodel(n*(t-T0),lam)
        Bnorm = Bnorm*Bmodel(-T0*n,lam)
    
    if renormalize:
        B = B/Bnorm

    return B
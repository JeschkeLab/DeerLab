# dipolarbackground.py - Multipathway background generator 
# --------------------------------------------------------
 # This file is a part of DeerLab. License is MIT (see LICENSE.md).
 # Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import inspect
import types

def dipolarbackground(t, pathways, Bmodel):
    r""" Constructs background decay functions according to the multi-pathway model.
    
    Parameters
    ----------
    t : array_like
        Time axis, in microseconds.

    pathways : list of lists or scalar
        List of pathways. Each pathway is defined as a list of the pathway's amplitude (lambda), refocusing time (T0), 
        and harmonic (n), i.e. ``[lambda, T0, n]`` or ``[lambda, T0]`` for one pathway. If n is not given it is assumed to be 1. 
        For a pathway with unmodulated contribution, only the amplitude must be specified, i.e. ``[Lambda0]``. 

    Bmodel : callable
        Background basis function. Its use depends on the type of background model. 

        * Physical background models: A callable function accepting a time-axis array as first input and a pathway amplitude as a second, e.g. ``B = lambda t,lam: bg_hom3d(t,par,lam)``
        * Phenomenological background models: A callable function accepting a time-axis array as input, e.g. ``B = lambda t: bg_exp(t,par)``

    Returns
    -------
    B : ndarray
        Multipathway background decay

    Notes
    -----
    Computes the multipathway background ``B`` for the time axis ``t`` corresponding to the dipolar pathways specified in ``pathways`` [1]_. 
    The total background is computed from the basis background function model specified in ``Bmodel``. For a multi-pathway DEER signal
    (e.g, 4-pulse DEER with 2+1 contribution; 5-pulse DEER with 4-pulse DEER residual signal, and more complicated experiments), ``pathways``
    is a 2D-array that contains a list of modulation depths (amplitudes), refocusing times (in microseconds), and optional harmonics for all modulated pathway signals

    References
    ----------
    .. [1] L. Fábregas Ibáñez, G. Jeschke, and S. Stoll. 
        DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson., 1, 209–224, 2020 

    Examples
    --------
    To calculate the background for the standard model for 4-pulse DEER with an unmodulated offset and a single dipolar pathway that refocuses at time 0, use::

        import numpy as np
        import deerlab as dl

        t = np.linspace(-5,20,501) # time axis (µs)
        lam = 0.4 # modulation depth main signal
        conc = 200   # spin concentration (µM)

        pathways = []
        path0 = [1-lam]             # unmodulated part, gives offset
        path1 = [lam, 0]            # main modulation, refocusing at time zero
        pathways = [path0, path1]

        def Bmodel(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t,pathways,Bmodel)
    """

    # Ensure that all inputs are numpy arrays
    t = np.atleast_1d(t)

    if type(Bmodel) is not types.LambdaType:
        raise TypeError('For a model with multiple modulated pathways, B must be a function handle of the type: @(t,lambda) bg_model(t,par,lambda)')

    # Identify phenomenological background models
    takes_lambda = len(inspect.signature(Bmodel).parameters)>1
    if not takes_lambda:
        _Bmodel = Bmodel
        Bmodel = lambda t,_: _Bmodel(t)

    pathways = [np.atleast_1d(path) for path in pathways]
    
    # Remove unmodulated pathways    
    unmodulated = [pathways.pop(i) for i,path in enumerate(pathways) if len(path)==1]
 
    # Check structure of pathways
    for i,path in enumerate(pathways):
        if len(path) == 2:
            # If harmonic is not defined, append default n=1
            pathways[i] = np.append(path,1) 
        elif len(path) != 3:
            # Otherwise paths are not correctly defined
            raise KeyError(f'The pathway #{i} must be a list of two or three elements [lam, T0] or [lam, T0, n]')

    # Construction of multi-pathway background function 
    #-------------------------------------------------------------------------------
    B = 1
    for pathway in pathways:        
        n,T0,lam = pathway
        B = B*Bmodel(n*(t-T0),lam)

    return B
    
# dipolarbackground.py - Multipathway background generator
# --------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2023: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import inspect

def dipolarbackground(t, pathways, basis):
    r"""Computes the (multi-pathway) background function.

    The function computes the two-spin interaction dipolar background which accounts for the structure and shape of the echo modulation 
    due to intermolecular dipolar interactions over a sample (background contribution). The dipolar background is constructed from the 
    contributions arising from different dipolar pathways [1]_ during a pulse sequence. The two-spin dipolar background reads:

    .. math::

        B(\mathbf{t}) = \prod_k B_\mathrm{basis} \Big( \lambda_k, \mathbf{\delta}_k (\mathbf{t} - \mathbf{t}_{\mathrm{ref},k}) \Big)

    where `\mathbf{t}=(t_1,\dots,t_D)` are the experimental time incrementation vectors of a `D`-dimensional experiment. 
    The dipolar pathways are characterized by their amplitudes `\lambda_k`, their refocusing times `\mathbf{t}_{\mathrm{ref},k}`, 
    and their (sub)harmonic factors `\mathbf{\delta}_k`. The function `B_\mathrm{basis}(\lambda,t)` represents the background basis function that
    describes the shape of the intermolecular contributions for a single dipolar pathway.

    The function can also compute three-spin interaction dipolar background functions that account for the additional the structure and shape
    of the echo modulation in multi-spin systems [2]_. The three-spin dipolar background is constructed from the contributions arising from 
    the combinations of different pair dipolar pathways [2]_ during a pulse sequence. 
    The three-spin dipolar background function reads:

    .. math::

        B^{(3)}(\mathbf{t}) = \prod_k\prod_q^3 B_\mathrm{basis}\Big(\lambda_k, \mathbf{\delta}_{k,q} (\mathbf{t} - \mathbf{t}_{\mathrm{ref},k,q}) \Big)

    The three-spin dipolar pathways are now characterized by an amplitude `\lambda_k`, and by a set of three refocusing times
     `\mathbf{t}_{\mathrm{ref},k,q}` and three (sub)harmonic factors `\mathbf{\delta}_{k,q}`.   

    Parameters
    ----------
    t : ndarray or list thereof
        Experimental time incrementation vectors (time coordinates) `t`, in microseconds. For multi-dimensional experiments, 
        specify a list ``[t1,...,tD]`` of the time coordinates `\mathbf{t}` along each of the ``D`` different experimental dimensions. 

    pathways : list of dictionaries  or ``None``, optional
        List of dipolar pathways. 
        
        For *two-spin dipolar interactions*, each pathway is defined as a dictionary containing the pathway's amplitude :math:`\lambda_k` (``'amp'``), 
        refocusing time :math:`t_{\mathrm{ref},k}` (``'reftime'``) in microseconds, and (sub)harmonic :math:`\delta_k` (``'harmonic'``).
        For example, for a one-dimensional experiment ``{'amp:'lambda, 'reftime':tref, 'harmonic':delta}`` or ``{'amp:'lambda, 'reftime':tref}``. If ``'harmonic'`` is not
        specified, it is assumed to be 1 (modulated). For an unmodulated pathway, specify only the amplitude, i.e. ``{'amp':Lambda0}``.
        
        For *multi-dimensional experiments*, the pathway refocusing times \mathbf{\t}_{\mathrm{ref},k}` and harmonics `\mathbf{\delta}_k` 
        along the different dimensions must be specified as a list, e.g. for a pathway in a two-dimensional experiment 
        ``{'amp:'lambda, 'reftime':[tref1,tref2], 'harmonic':[delta1,delta2]}``
        `
        For *three-spin interactions*, each pathway is defined as a dictionary containing the pathway's amplitude :math:`\lambda_k` (``'amp'``),
        a tuple of refocusing times :math:`t_{\mathrm{ref},k,q}` (``'reftime'``), and a tuple of (sub)harmonics :math:`\delta_{k,q}` (``'harmonic'``). 
        For example, for a one-dimensional experiment ``{'amp:'lambda, 'reftime':(tref1,tref2,tref3), 'harmonic':(delta1,delta2,delta3)}``. 
        
        If neither ``pathways`` or ``mod`` are specified (or set to ``None``), ``pathways=[{'amp':1,'reftime':0}]`` is used as the default.
    
    basis : callable or array_like or ``None``, optional
        Dipolar background basis function. Must be a callable function, accepting a time-vector array as first input and a pathway amplitude 
        as a second, i.e. ``bg = lambda t,lam: bg_model(t,par,lam)``. If set to ``None``, no background function is included.
        For a single-pathway model, the numerical background decay can be passed as an array. 

    Returns
    -------
    B : ndarray
        Background decay. 

    References
    ----------
    .. [1] L. Fábregas Ibáñez, G. Jeschke, and S. Stoll.
        DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson., 1, 209–224, 2020

    .. [2] L. Fábregas Ibáñez, M. H. Tessmer, G. Jeschke, and S. Stoll. 
        Dipolar pathways in multi-spin and multi-dimensional dipolar EPR spectroscopy, Phys. Chem. Chem. Phys., 24 2022, 22645-22660

    Examples
    --------
    To specify single-pathway background based on a homogenous 3D distribution of spins that refocuses at time `t=0`, use::

        lam = 0.4  # modulation depth
        pathways = [
            {'amp': lam, 'reftime': 0}
            ]

        def basis(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t, pathways, basis)

    To specify a three-pathway kernel of the 4-pulse DEER experiment that, e.g that includes the 2+1 contribution, use::

        tau1 = 0.5  # First interpulse delay (µs) 
        tau2 = 4.0  # Second interpulse delay (µs)
        
        lam1 = 0.40  # Pathway #1 amplitude
        lam2 = 0.10  # Pathway #2 amplitude
        lam3 = 0.10  # Pathway #3 amplitude

        tref1 = tau1       # Pathway #1 refocusing time
        tref2 = tau1+tau2  # Pathway #2 refocusing time
        tref3 = 0          # Pathway #3 refocusing time
        
        # Define the list of dipolar pathways
        pathways = [
            {'amp': lam1, 'reftime': tref1},
            {'amp': lam2, 'reftime': tref2},
            {'amp': lam3, 'reftime': tref3},
        ]

        def basis(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t, pathways, basis)


    To specify a dipolar background for an experiment on a three-spin system based on a single-pathway model use::

        tau1 = 0.5  # First interpulse delay (µs)         
        lam = 0.40  # Pathway #1 amplitude
        
        # Three dipolar pathways needed, one for each dipolar interaction in the three-spin system 
        pathways = [
            {'amp': lam, 'reftime': (tau1,None,None), 'harmonic': (1,0,0)},
            {'amp': lam, 'reftime': (None,tau1,None), 'harmonic': (0,1,0)},
            {'amp': lam, 'reftime': (None,None,tau1), 'harmonic': (0,0,1)},
        ]
        
        def basis(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t, pathways, basis)


    To construct a background function for a three-spin system studied with a two-dimensional experiment, with a single dipolar pathway refocusing at, e.g.,  `t=(\tau_1,\tau_2)` use:: 

        # Three dipolar pathways needed, one for each dipolar interaction in the three-spin system 
        pathways = [
            {'amp': 1-3*lam},
            {'amp': lam, 'reftime': ([tau1,tau2],None,None), 'harmonic': ([1,1],0,0)},
            {'amp': lam, 'reftime': (None,[tau1,tau2],None), 'harmonic': (0,[1,1],0)},
            {'amp': lam, 'reftime': (None,None,[tau1,tau2]), 'harmonic': (0,0,[1,1])},
        ]
        
        def basis(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t, pathways, basis)

    """

    # Mare sure that t is a list
    if not isinstance(t,(list)): t = [t]
    D  = len(t) # Number of experimental time-coordinates 
    t = [np.expand_dims(np.atleast_1d(t[d]), axis=tuple(np.delete(np.arange(D), d))) for d in range(D)]

    if not callable(basis):
        raise TypeError(
            "basis must be a function that can be called as basis(t,lam) or basis(t)."
        )
    nArgs = len(inspect.signature(basis).parameters)
    if nArgs > 2 or nArgs == 0:
        raise TypeError(
            "basis must be a function that can be called as basis(t,lam) or basis(t)."
        )

    # Identify physical background model
    takes_lambda = nArgs == 2

    # Ensure correct data types of the pathways and its components
    #pathways = [{key: np.atleast_1d(value).astype(float) for key,value in pathway.items()} for pathway in pathways]

    # Check structure of pathways
    for k,pathway in enumerate(pathways):
        if not 'amp' in pathway.keys():
            raise SyntaxError('All pathways must contain at least an \'amp\' field.')
        if not 'reftime' in pathway.keys():
            pathways[k]['reftime'] = ([None]*len(t)) 
            pathways[k]['harmonic'] = (np.zeros(D))
        for key in ['reftime','harmonic']:
            if key in pathway.keys():
                if not isinstance(pathway[key], tuple):
                    pathway[key] = tuple([pathway[key]])
                pathway[key] = tuple([np.atleast_1d(tref) for tref in pathway[key]])
        if 'amp' in pathway.keys() and 'reftime' in pathway.keys() and not 'harmonic' in pathway.keys():
            # If harmonic is not defined, append default n=1
            pathways[k]['harmonic'] = tuple([np.ones(D)]*len(pathways[k]['reftime']))
    # Remove unmodulated pathways
    [pathways.pop(k) for k,pathway in enumerate(pathways) if np.all(pathway['harmonic']==np.zeros(len(t)))]

    # Construction of multi-pathway background function
    B = 1
    for pathway in pathways:
        λ,tref,δ = [pathway['amp'],pathway['reftime'],pathway['harmonic']]
        for δq,trefq in zip(δ,tref):
            if len(trefq)!=len(t): 
                raise SyntaxError(f"Some pathway's number of refocusing times and harmonics do not match the number of time coordinates.")
            # Set trefs defined as None to an arbitrary numerical value
            trefq = [trefq[d] if δqd!=0 else 0 for d,δqd in enumerate(δq)]
            tdip = np.sum(np.array([δ_qd*(t_d-tref_qd) for t_d,δ_qd,tref_qd in zip(t,δq,trefq)],dtype=object), axis=0).astype(float)
            if takes_lambda:
                B *= basis(tdip, λ)
            else:
                B *= basis(tdip)

    return B

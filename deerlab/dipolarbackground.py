# dipolarbackground.py - Multipathway background generator
# --------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import inspect


def dipolarbackground(t, pathways, Bmodel):
    r"""Calculate multi-pathway dipolar background decays.

    Parameters
    ----------
    t : array_like
        Time axis, in microseconds.

    pathways : list of lists, or scalar
        List of dipolar pathways. Each pathway is represented by ``[lam, t0, delta]``, where ``lam`` is the pathway amplitude, ``t0``
        the refocusing time and ``delta`` the harmonic. If ``[lam, t0]`` is given, ``delta`` is set to zero. To specify the
        unmodulated contribution, use ``[lam]`` or just ``lam``.

    Bmodel : callable
        Background function, called either as ``Bmodel(t,lam)`` (for physical background models) or ``Bmodel(t)`` (for
        phenomenological background models), where ``t`` is the time axis and ``lam`` is the pathway amplitude.

    Returns
    -------
    B : 1D ndarray
        Background decay

    Notes
    -----
    Computes the total multi-pathway background decay function ``B`` for the time axis ``t`` for the set of dipolar pathways
    specified in ``pathways`` [1]_. The background for each pathway is calculated by calling ``Bmodel`` with the appropriate
    ampltidue, zero time and harmonic. The total background is the product over these single-pathway backgrounds.

    References
    ----------
    .. [1] L. Fábregas Ibáñez, G. Jeschke, and S. Stoll.
        DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson., 1, 209–224, 2020

    Examples
    --------
    To calculate the background for 4-pulse DEER with an unmodulated offset and a single  pathway that refocuses at time 0, use::

        import numpy as np
        import deerlab as dl

        t = np.linspace(-5,20,501)  # time axis (µs)
        lam = 0.4                   # modulation depth of main modulated signal
        conc = 200                  # spin concentration (µM)

        path0 = [1-lam]             # unmodulated part, gives offset
        path1 = [lam, 0]            # main modulation, refocusing at time zero
        pathways = [path0, path1]

        def Bmodel(t,lam):
            return dl.bg_hom3d(t,conc,lam)

        B = dl.dipolarbackground(t, pathways, Bmodel)
    """

    # Ensure that all inputs are numpy arrays
    t = np.atleast_1d(t)

    if not callable(Bmodel):
        raise TypeError(
            "Bmodel must be a function that can be called as Bmodel(t,lam) or Bmodel(t)."
        )
    nArgs = len(inspect.signature(Bmodel).parameters)
    if nArgs > 2 or nArgs == 0:
        raise TypeError(
            "Bmodel must be a function that can be called as Bmodel(t,lam) or Bmodel(t)."
        )

    # Identify physical background model
    takes_lambda = nArgs == 2

    # Extract modulated pathways, skipping unmodulated ones
    pathways = [path for path in pathways if len(np.atleast_1d(path)) > 1]

    # Check structure of pathways
    for i, path in enumerate(pathways):
        if len(path) == 2:
            # If harmonic is not defined, append default delta=1
            pathways[i] = np.append(path, 1)
        elif len(path) != 3:
            # Otherwise paths are not correctly defined
            raise KeyError(
                f"The pathway #{i} must be a list of two or three elements [lam, T0] or [lam, T0, delta]"
            )

    # Construction of multi-pathway background function
    B = 1
    if takes_lambda:
        # Physical background models
        for lam, t0, delta in pathways:
            B *= Bmodel(delta * (t - t0), lam)
    else:
        # Phenomenological background models
        for _, t0, delta in pathways:
            B *= Bmodel(delta * (t - t0))

    return B

# distancerange.py - DeerAnalysis distance axis constructor
# ------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.
import numpy as np
from deerlab.utils import isempty

def distancerange(t, nr=None):
    r""" 
    Empirical distance range given a DEER time axis

    Parameters
    ----------
    t : array_like
        Time axis, in microseconds.
    nr : scalar, integer
        Length of output distance axis. If not given, only min and max distance are returned.

    Returns
    -------
    r : ndarray or tuple
        Distance axis, in nanometers, running between empirical lower and upper limits ``rmin`` and ``rmax``.
        Either an ndarray (if ``nr`` is given) or a tuple (if ``nr`` is not given).

    Notes
    -----
    The minimal and maximal distances, ``rmin`` and ``rmax``, are empirical values that determine the
    minimal and maximal distance for which the given time trace can provide reliable information.

    The minimum distance is determined by the time step :math:`\Delta t` and the Nyquist criterion:

    .. math:: 
        
        r_\text{min} = \left( \frac{4\Delta t \nu_0}{0.85} \right)^{1/3}

    The maximum distance is determined by the requirement that at least half an oscillation
    should be observable over the measured time window from :math:`t_\text{min}`
    to :math:`t_\text{max}`.

    .. math:: 
        
        r_\text{max} = 6\left( \frac{t_\text{max}-t_\text{min}}{2} \right)^{1/3}

    where :math:`\nu_0` = 52.04 MHz nm^3.
    
    See Jeschke et al, Appl. Magn. Reson. 30, 473-498 (2006), https://doi.org/10.1007/BF03166213

    """

    t = np.atleast_1d(t)

    D = 52.04  # MHz nm^3
    
    # Minimum distance is determined by maximum frequency detectable with
    # the given time increment, based on Nyquist
    dt = np.mean(np.diff(t)) # time increment
    nupara_max = 1/2/dt  # maximum parallel dipolar frequency (Nyquist)
    nu_max = nupara_max/2  # maximum perpendicular dipolar frequency
    nu_max = nu_max*0.85  # add a bit of buffer
    rmin = (D/nu_max)**(1/3)
    
    # At least half a period of the oscillation
    # should be observable in the time window.
    trange = np.max(t) - np.min(t)
    Tmax = trange*2  # maximum period length
    rmax = (D*Tmax)**(1/3)

    if nr is not None:
        r = np.linspace(rmin, rmax, nr)
    else:
        r = (rmin, rmax)

    return r

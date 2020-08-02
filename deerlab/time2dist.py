# time2dist.py - DeerAnalysis distance axis constructor
# ------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.
import numpy as np
from deerlab.utils import isempty

def time2dist(t, M=[]):
    r""" 
    DeerAnalysis conversion from time-axis to distance-axis

    .. warning:: This is a legacy function. Its use is not recommended for routine or accurate data analysis.

    Parameters
    ----------
    t : array_like
        Time axis, in microseconds.
    M : int scalar
        Length of output distance axis, by default equal length as time axis.

    Returns
    -------
    r : ndarray
        Distance axis, in nanometers.

    Notes
    -----
    The minimal and maximal distances (rmin,rmax) are determined by the empirical
    approximations derived by Gunnar Jeschke as implemented in DeerAnalysis.

    These empirical equation approximate the minimal and maximal detectable distances
    given a certain time step :math:`\\Delta t` and trace length :math:`t_\\text{max}`.

    .. math:: r_\\text{min} = 4\left( \\frac{4\Delta t \\nu_0}{0.85} \\right)^{1/3}

    .. math:: r_\\text{max} = 6\left( \\frac{t_\\text{max}}{2} \\right)^{1/3}

    where :math:`\\nu_0` = 52.04 MHz is the dipolar frequency of between two nitroxide electron spins separated by 1 nm.

    """

    t = np.atleast_1d(t)

    if isempty(M):
        M = len(t)

    t = np.abs(t)
    dt = np.mean(np.abs(np.diff(t)))
    tmax = np.max(t)
    ny0 = 52.04
    rmin = (4*dt*ny0/0.85)**(1/3)
    rmax = 6*(tmax/2)**(1/3)

    r = np.linspace(rmin,rmax,M)

    return r

import numpy as np
from deerlab.utils import isempty

def correctphase(V, phase='posrealint', full_output=False):
# ==========================================================================
    r"""
    Phase correction of complex-valued data

    Performs a phase optimization on the complex-valued data ``V`` by determining a phase
    rotation of ``V`` that minimizes its imaginary component.
    
    Two-dimensional datasets ``V2D``, e.g. from measurements with multiple scans, can be provided, 
    and the phase correction will be done on each trace individually. The first dimension ``V2D[:,i]``
    must contain the single traces. An array of phases ``phases`` can be specified to manually correct the traces.

    Parameters
    ----------
    V : array_like, or list of array_like
        Complex-valued signals or list thereof.

    phase : string, optional
        Criterion for selection of correction phase. 

        * ``'posrealint'`` - Select the phase that gives the largest positive integral of the real part.
        * ``'negrealint'`` - Select the phase that gives the largest negative integral of the real part.
        * ``'close'`` - Select the phase closest to the average phase of the original data.

        The default behaviour is ``'posrealint'``.

    full_output : boolean, optional
        If True, the function will return additional output arguments, by default False.

    Returns
    -------
    Vr : ndarray
        Real part of the phase corrected dataset.

    Vi : ndarray (if full_output==True)
        Imaginary part of the phase corrected dataset.

    phase : float scalar or ndarray (if full_output==True)
        Fitted phase used for correction, in radians.

    """

    V_2d = V.copy()
    if V_2d.ndim == 1:
        V_2d = V_2d[:, np.newaxis]

    V_2d = V_2d.T

    # Calculate 3 points of cost function which should be a smooth, continuous sine wave with a frequency of 2 * phi
    phis = np.array([0, np.pi / 2, np.pi]) / 2
    costs = np.imag(V_2d[:, None] * np.exp(1j * phis)[None, :, None])
    costs = (costs * costs).sum(axis=-1)

    # Calculate sine function fitting 3 points
    offset = (costs[:, 0] + costs[:, 2]) / 2
    phase_shift = np.arctan2(costs[:, 0] - offset, costs[:, 1] - offset)
    amp = np.sqrt((costs[:, 0] - offset) ** 2 + (costs[:, 1] - offset) ** 2)

    # Calculate minima by setting the first derivative 0
    phis = (3 * np.pi / 2 - phase_shift) / 2
    phis[amp < 0] -= np.pi / 2
    phis[phis < 0] += np.pi

    tempspec = V_2d * np.exp(1j * phis)[:, None]
    if phase == 'posrealint':
        phis[tempspec.sum(axis=1) < 0] += np.pi
    elif phase == 'negrealint':
        phis[tempspec.sum(axis=1) > 0] -= np.pi

    V_2d *= np.exp(1j * phis)[:, None]

    V_2d = np.squeeze(V_2d.T)

    # Output
    Vreal = np.real(V_2d)
    Vimag = np.imag(V_2d)

    # Map phase angle to [-pi,pi) interval
    phase = phis

    if full_output:
        return Vreal, Vimag, phase
    else:
        return Vreal

# ==========================================================================
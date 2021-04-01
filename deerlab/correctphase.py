import numpy as np
import scipy.optimize as opt
from deerlab.utils import isempty

def correctphase(V, phase='posrealint', full_output=False):
# ==========================================================================
    r"""
    Phase correction of complex-valued data

    Performs a phase optimization on the complex-valued data ``V`` by determining a phase
    rotation of ``V`` that minimizes the imaginary component.
    
    Two-dimensional datasets ``V2D``, e.g. from multiple scans measurements, can be provided, 
    and the phase correction will be done on each trace individually. The first dimension ``V2D[:,i]``
    must contain the single traces. An array of phases ``phases`` can be specified to manually correct the traces.

    Parameters
    ----------
    V : array_like or list of array_like
        Complex-valued signals or list thereof.

    phase : string, optional
        Criterion for selection of correction phase. 

        * ``'posrealint'`` - Select the phase that gives the largest positive integral of the real part.
        * ``'negrealint'`` - Select the phase that gives the largest negative integral of the real part.
        * ``'close'`` - Select the phase closest to the average phase of the original data.

        The default behaviour is ``'posrealint'``.

    full_output : boolean, optional
        If enabled, the function will return additional output arguments, by default disabled.

    Returns
    -------
    Vr : ndarray
        Real part of the phase corrected dataset.

    Vi : ndarray (if full_output==True)
        Imaginary part of the phase corrected dataset.

    phase : float scalar (if full_output==True)
        Fitted phase used for correction, in radians.    

    """

    V_2d = V.copy()
    if V_2d.ndim == 1:
        Ntraces = 1
        V_2d = V_2d[:,np.newaxis]
    else:
        Ntraces = V_2d.shape[1]
    n = V_2d.shape[0]
 
    phaseopt = np.zeros(Ntraces)
    for i in range(Ntraces):
        V_ = V_2d[:,i]

        # Get the average phase in the signal
        phaseav = np.mean(np.angle(V_))

        # Objective function for optimization of the phase
        phase_objfcn = lambda phase: np.linalg.norm(np.imag(V_*np.exp(-1j*phase)))

        pars = opt.fmin(phase_objfcn,phaseav,maxfun=1e5,maxiter=1e5,disp=False)
        
        # Get the two phases that minimize the imaginary part
        phases = [pars[0], pars[0]+np.pi]
        realint = [np.sum(V_*np.exp(-1j*phase)) for phase in phases]
            
        if phase == 'posrealint':
            phaseopt[i] =  phases[np.argmax(realint)]
        elif phase == 'negrealint':
            phaseopt[i] =  phases[np.argmin(realint)]
        elif phase == 'close':
            phaseopt[i] =  phases[np.argmin(abs(phases - phaseav))]

        # Apply phase
        V_2d[:,i] = V_*np.exp(-1j*phaseopt[i])

    if Ntraces==1:
        V_2d = np.squeeze(V_2d)

    # Output
    Vreal = np.real(V_2d)
    Vimag = np.imag(V_2d)
    # Map phase angle to [-pi,pi) interval
    phase = np.angle(np.exp(-1j*phaseopt)) 

    if full_output:
        return Vreal,Vimag,phase
    else:
        return Vreal

# ==========================================================================
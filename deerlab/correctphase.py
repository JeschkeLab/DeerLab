import numpy as np

def correctphase(V, full_output=False):
    r"""
    Phase correction of complex-valued data.

    Rotates the phase of complex-valued data ``V`` to minimize the imaginary component.
    Among the two phases that minimize the imaginary part, the one that gives a real
    part with a positive average is used.
    
    For two-dimensional datasets ``V2D``, e.g. from measurements with multiple scans,
    each slice ``V2D[:,i]`` is phase-rotated independently.
    

    Parameters
    ----------
    V : array_like, or list of array_like
        Complex-valued 1D or 2D signal.

    full_output : boolean, optional
        If ``True``, return additional output arguments. (default: ``False``)

    Returns
    -------
    Vr : ndarray
        Real part of the phase-corrected data.

    Vi : ndarray (only if ``full_output==True``)
        Imaginary part of the phase-corrected data.

    phase : float scalar or ndarray (only if ``full_output==True``)
        Fitted phase, or list of phases for 2D data, used for correction, in radians.

    """
    
    if not np.iscomplexobj(V):
        raise ValueError("Data set must be complex-valued.")

    data1d = V.ndim==1

    V_2d = V.copy()
    if data1d:
        V_2d = V_2d[:, np.newaxis]

    # The follwing determines the phase that minimizes the cost
    # function = sum of squares of imaginary part of V*exp(1j*phi)
    # This cost function has the analytical form
    #
    #   (A+B) + (B-A)*cos(2*phi) + C*sin(2*phi)
    #    = offset + amp*cos(2*phi-phi0)
    #
    # where
    #    A = sum_k real(V_k)^2 / 2
    #    B = sum_k imag(V_k)^2 / 2
    #    C = sum_k real(V_k)*imag(V_k)
    #
    #    offset = A+B
    #    amp = sqrt((B-A)^2+C^2)
    #    phi0 = atan2(C, B-A)
    #
    # The cost function has two minima:
    #    phi = phi0/2 + pi/2   and   phi = phi0/2 + 3*pi/2
    
    # Calculate phase that minimizes cost function
    Vr = np.real(V_2d)
    Vi = np.imag(V_2d)
    A = np.sum(Vr**2, axis=0)/2
    B = np.sum(Vi**2, axis=0)/2
    C = np.sum(Vr*Vi, axis=0)
    phi0 = np.arctan2(C, B-A)
    phimin = phi0/2 + np.pi/2  # one of the two minimizers
    
    # Apply phase rotation
    V_2d *= np.exp(1j*phimin)[None,:]
    
    # Pick minimizer that yields positive average of real part
    reAvg = np.average(V_2d, axis=0)
    idx = reAvg < 0
    phimin[idx] += np.pi
    V_2d[:,idx] = -V_2d[:,idx]

    # Assemble output
    if data1d:
        V_2d = np.squeeze(V_2d, axis=1)
    Vreal = np.real(V_2d)
    Vimag = np.imag(V_2d)
    if full_output:
        return Vreal, Vimag, phimin
    else:
        return Vreal

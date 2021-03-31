import numpy as np
import scipy.optimize as opt
from deerlab.utils import isempty

def correctphase(V, phase=None, imagoffset=False, selphase='maxrealint', full_output=False):
# ==========================================================================
    r"""
    Phase correction of complex-valued data

    Performs a phase optimization on the complex-valued data ``V`` by determining a phase
    rotation of ``V`` that minimizes the imaginary component.
    
    Also, the phase can be corrected manually by specifying  a phase ``phase``, in radians.

    Two-dimensional datasets ``V2D``, e.g. from multiple scans measurements, can be provided, 
    and the phase correction will be done on each trace individually. The first dimension ``V2D[:,i]``
    must contain the single traces. An array of phases ``phases`` can be specified to manually correct the traces.

    Parameters
    ----------
    V : array_like or list of array_like
        Complex-valued signals or list thereof.

    phase : float scalar, optional
        Phase shift for manual correction, in radians. 

    imagoffset : boolean, optional
        Enables/Disables the fitting and correction of an imaginary offset, by default disabled.

    selphase : string, optional 
        Selection criterion for phase optimization. 

        * ``'maxrealint'`` - Maximization of the integral of the real component of ``V``.
        * ``'minphase'`` - Minimization of the imaginary component of ``V``.

        The default behaviour is ``'maxrealint'``.

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

    imoffset : float scalar (if full_output==True)
        Fitted imaginary offset used for correction.

    """

    V_ = V.copy()
    if V_.ndim == 1:
        Ntraces = 1
        V_ = V_[:,np.newaxis]
    else:
        Ntraces = V_.shape[1]
    n = V_.shape[0]

    # Determine if phase must be fitted or has been passed
    if phase is None:
        fitPhase = True
    else: 
        fitPhase = False
        phase = np.atleast_1d(phase)
        if len(phase) != Ntraces:
            raise ValueError('The number of input phases must agree with the number of traces.') 
    if imagoffset:
        selphase=='minphase'

    def phase_objfcn(params,V):
    # ==========================================================================
        " Objective function for optimization of the phase"

        phase = params[0]
        if len(params)>1:
            imoffsets = params[1]
        else:
            imoffsets = 0

        Vcorr = (V-1j*imoffsets)*np.exp(-1j*phase)
        if selphase=='minphase':
            # Minimize the imaginary part / Maximize the real part
            objfcn = np.linalg.norm(np.imag(Vcorr))
        elif selphase=='maxrealint':
            # Maximize the integral of the real opart
            objfcn = np.sum(np.abs(Vcorr))-np.sum(np.real(Vcorr))
        else: 
            raise KeyError("The requested method is not valid. It must be either 'maxrealint' or 'minphase'.")

        return objfcn
    # ==========================================================================

    # Phase/offset fitting
    #-------------------------------------------------------------------------------
    ImagOffset = np.zeros(Ntraces)
    if fitPhase:
        phase = np.zeros(Ntraces)
        for i in range(Ntraces):
            par0 = []
            FitRange = np.arange(round(n/8),n) # use only last 7/8 of data for phase/offset correction
            V_cut = V_[FitRange,i]
            par0.append(np.mean(np.angle(V_cut))) # use average phase as initial value
            if imagoffset:
                par0.append(np.mean(np.imag(V_cut))) # use average offset as initial value

            pars = opt.fmin(lambda par: phase_objfcn(par,V_cut),par0,maxfun=1e5,maxiter=1e5,disp=False)
            phase[i] = pars[0]
            if imagoffset:
                ImagOffset[i] = pars[1]
    else:
        if imagoffset:
            # Fit only imaginary offset        
            for i in range(Ntraces):
                par0 = 0
                ImagOffset[i] = opt.fmin(lambda offset: phase_objfcn([phase, offset],V_cut[:,i]),par0,maxfun=1e5,maxiter=1e5,disp=False)

    ImagOffset = ImagOffset*1j

    # Apply phase/offset correction
    ph = np.exp(-1j*phase)
    Vc = (V_ - ImagOffset)*ph

    if Ntraces==1:
        Vc = np.squeeze(Vc)

    # Output
    Vreal = np.real(Vc)
    Vimag = np.imag(Vc)
    # Map phase angle to [-pi,pi) interval
    phase = np.angle(ph) 

    if full_output:
        return Vreal,Vimag,phase,ImagOffset
    else:
        return Vreal

# ==========================================================================
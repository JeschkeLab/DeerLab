
import numpy as np
import scipy.optimize as opt
from deerlab.utils import isempty

def correctphase(V, Phase = [], fitImagOffset=False, full_output=False):
# ==========================================================================
    r"""
    Phase correction of complex-valued data

    Performs a phase optimization on the complex-valued data ``V`` by minimization of the
    imaginary component of the data. The phase can be corrected manually by specifying 
    a phase ``phase``, in radians.

    Two-dimensional datasets ``V2D``, e.g. from multiple scans measurements, can be provided, 
    and the phase correction will be done on each trace individually. The first dimension ``V2D[:,i]``
    must contain the single traces. An array of phases ``phases`` can be specified to manually correct the traces.

    Parameters
    ----------
    V : array_like or list of array_like
        Complex-valued signals or list thereof.

    Returns
    -------
    Vr : ndarray
        Real part of the phase corrected dataset.
    Vi : ndarray (if full_output==True)
        Imaginary part of the phase corrected dataset.
    Phase : float scalar (if full_output==True)
        Fitted phase used for correction, in radians.    
    ImOffset : float scalar (if full_output==True)
        Fitted imaginary offset used for correction.

    Other Parameters
    ----------------
    Phase  : float scalar
        Phase shift for manual correction, in radians. 
    fitImagOffset : boolean
        Enables/Disables the fitting and correction of an imaginary offset, by default disabled.
    full_output : boolean
        If enabled the function will return additional output arguments in a tuple, by default disabled.

    """

    if V.ndim == 1:
        Ntraces = 1
        V = V[:,np.newaxis]
    else:
        Ntraces = V.shape[1]
    n = V.shape[0]

    # Determine if phase must be fitted or has been passed
    if isempty(Phase):
        fitPhase = True
    else: 
        fitPhase = False
        Phase = np.atleast_1d(Phase)
        if len(Phase) != Ntraces:
            raise ValueError('The number of input phases must agree with the number of traces.') 

    # Phase/offset fitting
    #-------------------------------------------------------------------------------
    ImagOffset = np.zeros(Ntraces)
    if fitPhase:
        Phase = np.zeros(Ntraces)
        for i in range(Ntraces):
            par0 = []
            FitRange = np.arange(round(n/8),n) # use only last 7/8 of data for phase/offset correction
            V_ = V[FitRange,i]
            par0.append(np.mean(np.angle(V_))) # use average phase as initial value
            if fitImagOffset:
                par0.append(np.mean(np.imag(V_))) # use average offset as initial value
            fun = lambda par: _imaginarynorm(par,V_)

            pars = opt.fmin(fun,par0,maxfun=1e5,maxiter=1e5,disp=False)
            Phase[i] = pars[0]
            if fitImagOffset:
                ImagOffset[i] = pars[1]
    else:
        if fitImagOffset:
            # Fit only imaginary offset        
            for i in range(Ntraces):
                par0 = 0
                fun = lambda offset: _imaginarynorm([Phase, offset],V[:,i])
                ImagOffset[i] = opt.fmin(fun,par0)

    ImagOffset = ImagOffset*1j

    # Apply phase/offset correction
    ph = np.exp(1j*Phase)
    Vc = (V - ImagOffset)/ph

    if Ntraces==1:
        Vc = np.squeeze(Vc)

    # Output
    Vreal = np.real(Vc)
    Vimag = np.imag(Vc)
    Phase = np.angle(ph) # map phase angle to [-pi,pi) interval

    if full_output:
        return Vreal,Vimag,Phase,ImagOffset
    else:
        return Vreal

# ==========================================================================


def _imaginarynorm(params,V):
# ==========================================================================
        """
        Computes norm of the imaginary part of phase-corrected data from zero before
        phase correction, an offset can be subtracted from the imaginary part.
        """

        phase = params[0]
        if len(params)>1:
            imoffsets = params[1]
        else:
            imoffsets = 0

        Vcorr = (V-1j*imoffsets)*np.exp(-1j*phase)
        ImagNorm = np.linalg.norm(np.imag(Vcorr))

        return ImagNorm
# ==========================================================================
    
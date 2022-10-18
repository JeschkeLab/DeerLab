# dipolarkernel.py - Dipolar kernel operator
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

# Numpy + SciPy
import numpy as np
from numpy import inf
from scipy.integrate import quad_vec, quad
from scipy.special import fresnel
from scipy.interpolate import make_interp_spline
# Other packages
import types
import warnings
from memoization import cached
# DeerLab dependencies
from deerlab.dipolarbackground import dipolarbackground
from deerlab.constants import *

def dipolarkernel(t, r, *, pathways=None, mod=None, bg=None, method='fresnel', excbandwidth=inf, orisel=None, g=None, 
                  integralop=True, nKnots=5001, complex=False, clearcache=False, memorylimit=8):
#===================================================================================================
    r"""Compute the (multi-pathway) dipolar kernel operator which enables the linear transformation from
    distance-domain to time-domain data. 

    Parameters
    ----------
    t : ndarray or list thereof
        Dipolar evolution time axis, in microseconds.
    
    r : ndarray or list thereof
        Distance axis, in nanometers.
    
    pathways : list of lists  or ``None``, optional
        List of dipolar pathways [1]_. Each pathway is defined as a list of the pathway's amplitude (lambda), refocusing time in microseconds (T0), 
        and harmonic (n), i.e. ``[lambda, T0, n]`` or ``[lambda, T0]``. If n is not given, it is assumed to be 1. 
        For a unmodulated pathway, specify only the amplitude, i.e. ``[Lambda0]``. If neither ``pathways`` or ``mod`` are specified
        (or ``None``), ``pathways=[[1,0]]`` is used as the default.
    
    mod : scalar  or ``None``, optional
        Modulation depth for the simplified 4-pulse DEER model. If neither ``pathways`` or ``mod`` are specified (or ``None``),
        ``mod=1`` is used as the default.
        
    bg : callable or array_like or ``None``, optional
        For a single-pathway model, the numerical background decay can be passed as an array. 
        For multiple pathways, a callable function must be passed, accepting a time-axis array as first input and a pathway amplitude as a second, i.e. ``B = lambda t,lam: bg_model(t,par,lam)``. If set to ``None``, no background decay is included.

    method : string, optional
        Numerical method for kernel matrix calculation: 

            * ``'fresnel'`` - uses Fresnel integrals for the kernel (fast, accurate)
            * ``'integral'`` - uses explicit integration function (slow, accurate)
            * ``'grid'`` - powder average via explicit grid integration (slow, inaccurate)
        
        The default is ``'fresnel'``.

    orisel : callable  or ``None``, optional 
        Probability distribution of possible orientations of the interspin vector to account for orientation selection. Must be 
        a function taking a value of the angle θ∈[0,π/2] between the interspin vector and the external magnetic field and returning
        the corresponding probability density. If specified as ``None`` (by default), a uniform distribution is assumed. 
        Requires the ``'grid'`` or ``'integral'`` methods.

    excbandwidth : scalar, optional
        Excitation bandwidth of the pulses in MHz to account for limited excitation bandwidth [5]_.
        If not specified, an infinite excitation bandwidth is used.
        Requires the ``'grid'`` or ``'integral'`` methods.

    g : scalar, 2-element array, optional
        Electron g-values of the spin centers ``[g1, g2]``. If a single g is specified, ``[g, g]`` is used.
        If not specified, g = 2.002319... is used for both spins.

    complex : boolean, optional 
        Return the complex-valued kernel such that the matrix operation ``V=K@P`` with a distance distribution yields the in-phase 
        and out-of-phase components of the dipolar signal ``V``. Disabled by default. 
        Requires the ``'fresnel'`` or ``'grid'`` methods.

    integralop : boolean, optional
        Return the kernel as an integral operator (i.e ``K = K*dr``) or not (``K``). Usage as an integral operator means that the 
        matrix operation ``V=K@P`` with a normalized distance distribution (i.e. ``trapz(r,P)==1``) leads to a signal ``V`` with 
        amplitude ``V(t=0)=1``. The default is ``True``.
    
    nKnots : scalar, optional
        Number of knots for the grid of powder orientations to be used in the ``'grid'`` kernel calculation method. By default set to 5001 knots.
        
    clearcache : boolean, optional
        Clear the cached dipolar kernels at the beginning of the function to save memory. Disabled by default.

    memorylimit : scalar, optional
        Memory limit to be allocated for the dipolar kernel. If the requested kernel exceeds this limit, the 
        execution will be stopped. The default is 12GB.  

    Returns
    --------
    K : ndarray
        Dipolar kernel operator matrix, such that for a distance distribution ``P``, the dipolar signal is ``V = K@P``

    Notes
    -----
    For a multi-pathway DEER [1]_, [2]_ signal (e.g, 4-pulse DEER with 2+1 contribution 5-pulse DEER with 4-pulse DEER residual signal, and more complicated experiments), ``pathways`` contains a list of pathway amplitudes and refocusing times (in microseconds).
    The background function specified as ``B`` is used as basis function, and the actual multipathway background included into the kernel is computed using :ref:`dipolarbackground`. The background in included in the dipolar kernel definition [3]_. 
    Optionally, the harmonic (1 = fundamental, 2 = first overtone, etc.) can be given as a third value in each row. This can be useful for modeling RIDME signals [4]_. If not given, the harmonic is 1 for all pathways. 


    Examples
    --------
    To specify the standard model for 4-pulse DEER with an unmodulated offset and a single dipolar pathway that refocuses at time 0, use::
    
        lam = 0.4  # modulation depth main signal
        pathways = [[1-lam], [lam, 0]]
        
        K = dl.dipolarkernel(t,r,pathways=pathways)


    A shorthand input syntax equivalent to this input::

        lam = 0.4
        K = dl.dipolarkernel(t,r,mod=lam)


	To specify a more complete 4-pulse DEER model that, e.g includes the 2+1 contribution, use::

        Lam0 = 0.5  # unmodulated part
        lam = 0.4  # modulation depth main signal
        lam21 = 0.1  # modulation depth 2+1 contribution
        tau2 = 4  # refocusing time (µs) of 2+1 contribution

        path0 = Lam0  # unmodulated pathway
        path1 = [lam1, 0]  # main dipolar pathway, refocusing at time zero
        path1 = [lam2, tau2]  # 2+1 dipolar pathway, refocusing at time tau2
        pathways = [path0, path1, path2]  
        
        K = dl.dipolarkernel(t,r,pathways=pathways)


    References
    ----------
    .. [1] L. Fábregas Ibáñez, M. H. Tessmer, G. Jeschke, and S. Stoll. 
        Dipolar pathways in dipolar EPR spectroscopy, Phys. Chem. Chem. Phys., 2022, Advance Article

    .. [2] L. Fábregas Ibáñez, G. Jeschke, and S. Stoll. 
        DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson., 1, 209–224, 2020 

    .. [3] L. Fábregas Ibáñez, and G. Jeschke
        Optimal background treatment in dipolar spectroscopy, Physical Chemistry Chemical Physics, 22, 1855–1868, 2020.

    .. [4] K. Keller, V. Mertens, M. Qi, A. I. Nalepa, A. Godt, A. Savitsky, G. Jeschke, and M. Yulikov
        Computing distance distributions from dipolar evolution data with overtones: RIDME spectroscopy with Gd(III)-based spin labels, Physical Chemistry Chemical Physics, 19

    .. [5] J. E. Banham, C. M. Baker, S. Ceola, I. J. Day, G.H. Grant, E. J. J. Groenen, C. T. Rodgers, G. Jeschke, C. R. Timmel
        Distance measurements in the borderline region of applicability of CW EPR and DEER: A model study on a homologous series of spin-labelled peptides, Journal of Magnetic Resonance, 191, 2, 2008, 202-218
    """
    # Clear cache of memoized function is requested
    if clearcache:
        elementarykernel.cache_clear()

    # If not specified, assume the free-electron g-values
    if g is None:
        g = [ge, ge]

    # Mare sure that r and t are lists
    if not isinstance(t,(list)): t = [t]
    if not isinstance(r,(list)): r = [r]

    D  = len(t) # Number of experimental time-coordinates 
    Q = len(r) # Number of dipolar interactions

    # Parse the individual time and distance dimensions 
    t = [np.expand_dims(np.atleast_1d(t[d]), axis=tuple(np.delete(np.arange(D), d))) for d in range(D)]
    for q,r_ in enumerate(r):
        r[q] = np.atleast_1d(r_)    
        if np.any(r[q]<=0):
            raise ValueError("All distances must be nonnegative and nonzero.")

    # Check that requested kernel does not exceed the memory limit
    memory_requirement = 8/1e9 # GB/float
    memory_requirement *= np.prod([len(t_) for t_ in t])*np.prod([len(r_) for r_ in r]) 
    if method=='grid': 
        memory_requirement = memory_requirement*nKnots
    if memory_requirement > memorylimit: 
        raise MemoryError(f'The requested kernel requires {memory_requirement:.2f}GB, exceeding the current memory limit {memorylimit:.2f}GB.')

    # Check validity of the g-values
    g = np.atleast_1d(g)
    if len(g) == 1:
        g = [g, g]
    elif len(g) != 2:
        raise TypeError('The g-value must be specified as a scalar or a two-element array.')

    # Check that the kernel construction method is compatible with requested options
    if method=='fresnel':
        if Q>1: 
            raise KeyError("The 'fresnel' method is not compatible with multiple dipolar interactions. Use the 'grid' or 'integral' methods.")
        if excbandwidth != inf:
            raise KeyError("Excitation bandwidths can only be specified with the 'grid' or 'integral' methods.")
        if orisel is not None:
            raise KeyError("Orientation selection weights can only be specified with  the 'grid' or 'integral' methods.")
    if method=='integral':
        if complex:
            raise KeyError("Complex-valued kernels cannot be computed by the 'integral' method.")

    # Check whether the full pathways or the modulation depth have been passed
    pathways_passed =  pathways is not None
    moddepth_passed =  mod is not None
    if pathways_passed and moddepth_passed: 
        raise KeyError("Modulation depth and dipolar pathways cannot be specified simultaneously.")    

    if moddepth_passed:
        # Construct single-pathway kernel if modulation depth is specified
        unmod = {'amp': 1 - mod}
        singpathway = {'amp': mod, 'reftime': tuple([np.zeros(D)]*Q), 'harmonic': tuple([np.ones(D)]*Q)}
        pathways = [unmod,singpathway]
    elif not pathways_passed:
        # Full modulation if neither pathways or modulation depth are specified
        unmod = {'amp': 0}
        singpathway = {'amp': 1, 'reftime': tuple([np.zeros(D)]*Q), 'harmonic': tuple([np.ones(D)]*Q)}
        pathways = [unmod,singpathway]

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

        for trefq in pathways[k]['reftime']:
            if len(trefq)!=len(t): 
                raise SyntaxError(f"Some pathway's number of refocusing times and harmonics do not match the number of time coordinates.")


    # Get unmodulated pathways    
    unmodulated = [pathways.pop(k) for k,pathway in enumerate(pathways) if np.all(pathway['harmonic']==np.zeros(len(t)))]
    # Combine all unmodulated contributions
    if unmodulated:
        Λ0 = sum([pathway['amp'] for pathway in unmodulated])
    else: 
        Λ0 = 0

    # Define kernel matrix auxiliary function
    K0 = lambda tdip,r,method: elementarykernel(tdip,r,method,excbandwidth,nKnots,g,orisel,complex)

    # Build dipolar kernel matrix, summing over all pathways
    K = Λ0
    for pathway in pathways:
        λ,tref,δ = [pathway['amp'],pathway['reftime'],pathway['harmonic']]

        # Set trefs of unmodulated patwhays to an arbitrary numerical value
        trefq = [[trefq[d] if δqd!=0 else 0 for d,δqd in enumerate(δq)] for δq,trefq in zip(δ,tref)]

        # Construct multi-dimensional effective dipolar evolution time
        tdip = [np.sum(np.array([δ_qd*(t_d-tref_qd) for t_d,δ_qd,tref_qd in zip(t,δq,trefq)],dtype=object), axis=0).astype(float) for δq,trefq in zip(δ,tref)]

        # Determine number of spins participating in the pathway
        Nspin = np.sum(np.abs(δ))+1 

        # Two-spin dipolar interactions
        if Nspin==2 and method!='grid':
            q = int(np.where([np.any(δq!=0) for δq in δ])[0])
            K += λ*K0(tdip[q],r[q],method)
        else:
            K += λ*K0(tdip,r,'grid')

    # Multiply by background if given
    if bg is not None:
        if type(bg) is types.LambdaType:
            B = dipolarbackground(t,pathways,bg)
        else:
            B = np.atleast_1d(bg)
        K *= B[:,np.newaxis]

    # Include Δr factor for integration
    for r_q in r:
        if integralop and len(r_q)>1:
            # Vectorized matrix form of trapezoidal integration method
            dr = np.zeros_like(r_q)
            dr[0] = r_q[1] - r_q[0]
            dr[-1] = r_q[-1] - r_q[-2]
            dr[1:-1] = r_q[2:] - r_q[0:-2]
            dr = dr/2
            K *= dr
    
    return K
#==============================================================================

def _echomodulation(ωr,tdip,cosθ,ωex,complex):
#==============================================================================
    "Spin echo modulation due to dipolar coupling(s)"
    echomod = 1
    # Loop over the different dipolar interactions
    for ωr_q, tdip_q in zip(ωr,tdip):
        # Differentiate between complex and real valued cases for a speedup
        if complex: 
            echomod *= np.exp(-1j*ωr_q*tdip_q*(1 - 3*cosθ**2))
        else:
            echomod *= np.cos(ωr_q*abs(tdip_q)*(1 - 3*cosθ**2))
        # If given, include limited excitation bandwidth
        if not np.isinf(ωex):
            echomod *= np.exp(-(ωr_q*(1 - 3*cosθ**2))**2/ωex**2)
    return echomod
#==============================================================================

@cached(max_size=100)
def elementarykernel(tdip,r,method,ωex,nKnots,g,Pθ,complex):
#==============================================================================
    "Calculates the elementary dipolar kernel (cached for speed)"

    # Work with vectors
    r = np.atleast_1d(r)

    # Dipolar frequencies
    ωr = [(μ0/2)*μB**2*g[0]*g[1]/h*1e21/(r_q**3) for r_q in r]  # rad μs^-1

    # Orientation selection
    orientationselection = Pθ is not None 
    if orientationselection:
        # Ensure zero-derivatives at [0,π/2]
        θ = np.linspace(0,π/2,50) # rad
        Pθ_ = make_interp_spline(θ, Pθ(θ),bc_type="clamped")
        # Ensure normalization of probability density function (∫P(cosθ)dcosθ=1)
        Pθnorm,_ = quad(lambda cosθ: Pθ_(np.arccos(cosθ)),0,1,limit=1000)
        Pθ = lambda θ: Pθ_(θ)/Pθnorm

    def elementarykernel_fresnel(tdip):
    #==========================================================================
        """Calculate kernel using Fresnel integrals (fast and accurate)"""

        # Calculation using Fresnel integrals
        if complex:
            ɸ = np.multiply.outer(tdip,ωr)
            κ = np.lib.scimath.sqrt(6*ɸ/π)
            S, C = fresnel(κ)
            K0 = np.exp(-1j*ɸ)*(C + 1j*S)

        else:
            ɸ = np.multiply.outer(np.abs(tdip),ωr)
            κ = np.lib.scimath.sqrt(6*ɸ/π)
            S, C = fresnel(κ)
            K0 = C*np.cos(ɸ) + S*np.sin(ɸ)

        # Supress divide by 0 warning       
        with warnings.catch_warnings(): 
            warnings.simplefilter("ignore")
            K0 = K0/κ

        K0 = np.where(ɸ==0,1,K0)
        return K0
    #==========================================================================

    def elementarykernel_grid(tdip,ωex,nKnots,Pθ):
    #==========================================================================
        """Calculate kernel using grid-based powder integration (converges very slowly with nKnots)"""

        # Orientational grid
        cosθ = np.linspace(0,1,int(nKnots))

        # Vectorized 3D-grid evaluation (t,r,orientations)
        D,Q = len(tdip[0].shape),len(ωr)
        dims = np.arange(D + Q + 1)
        # Expand dimensions of all variables on a multi-dimensional grid
        tdip_ = [np.expand_dims(tdip[q], axis=tuple(np.delete(dims, range(D)))) for q in range(Q)]
        ωr_   = [np.expand_dims(ωr[q],   axis=tuple(np.delete(dims, D + q))) for q in range(Q)]
        cosθ_ =  np.expand_dims(cosθ,    axis=tuple(np.delete(dims, D + Q)))

        # Integrate over the grid of echo modulations
        if orientationselection:
            # With orientation distribution Pθ
            K0 = np.dot(_echomodulation(ωr_,tdip_,cosθ_,ωex,complex),Pθ(np.arccos(cosθ)))/nKnots
        else: 
            # Without orientation selection
            K0 = np.sum(_echomodulation(ωr_,tdip_,cosθ_,ωex,complex),axis=-1)/nKnots
        return K0
    #==========================================================================

    def elementarykernel_integral(tdip,ωex,Pθ):
    #==========================================================================
        """Calculate kernel using explicit numerical integration """
        K0 = np.zeros((np.size(tdip),np.size(ωr)))
        for n,ωr_ in enumerate(ωr):
            #==================================================================
            def integrand(cosθ):
                integ = _echomodulation([ωr_,0,0],[tdip,0,0],cosθ,ωex,complex)
                # If given, include orientation selection
                if orientationselection:
                    integ = integ*Pθ(np.arccos(cosθ))  
                return integ
            #==================================================================   
            K0[:,n] = quad_vec(integrand,0,1,limit=1000)[0].squeeze()
        return K0
    #==========================================================================

    if method=='fresnel':
        K0 = elementarykernel_fresnel(tdip)
    elif method=='integral': 
        K0 = elementarykernel_integral(tdip,ωex,Pθ)
    elif method=='grid': 
        K0 = elementarykernel_grid(tdip,ωex,nKnots,Pθ)
    return K0
#==============================================================================

# dipolarkernel.py - Dipolar kernel operator
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2023: Luis Fabregas, Stefan Stoll and other contributors.

# Numpy + SciPy
import numpy as np
from numpy import inf
from scipy.integrate import quad_vec, quad
from scipy.special import fresnel
from scipy.interpolate import make_interp_spline,interp1d
# Other packages
import numexpr as ne
import types
import warnings
from memoization import cached
# DeerLab dependencies
from deerlab.dipolarbackground import dipolarbackground
from deerlab.utils import sophegrid
from deerlab.constants import *

def dipolarkernel(t, r, *, pathways=None, mod=None, bg=None, method='fresnel', excbandwidth=inf, orisel=None, g=None, 
                  integralop=True, gridsize=1000, complex=False, clearcache=False, memorylimit=8, tinterp=None):
#===================================================================================================
    r"""
    Compute the (multi-pathway) dipolar kernel operator. 

    The function computes the two-spin interaction dipolar kernel which accounts for the structure and shape of the echo modulation 
    due to intra- and inter-molecular dipolar interactions over a sample. The dipolar kernel is constructed from the 
    contributions arising from different dipolar pathways [1]_ during a pulse sequence. The two-spin dipolar kernel reads:

    .. math::

        K(\mathbf{t},r) = \Lambda_0 + \sum_k \lambda_k 
        \int_0^\pi \mathrm{d}\theta \sin(\theta) 
        \exp\Bigg(-i \frac{D\mathbf{\delta}_k (\mathbf{t} - \mathbf{t}_{\mathrm{ref},k})}{r^3}(1-3\cos(\theta)^2) \Bigg)
        \prod_k B(\lambda_k, \mathbf{\delta}_k (\mathbf{t} - \mathbf{t}_{\mathrm{ref},k}))

    where `\mathbf{t}=(t_1,\dots,t_D)` are the experimental time incrementation vectors of a `D`-dimensional experiment
    and `r` represents the interspin distance. The dipolar pathways are characterized by their amplitudes `\lambda_k`,  
    their refocusing times `\mathbf{t}_{\mathrm{ref},k}`, and their (sub)harmonic factors `\mathbf{\delta}_k`. The constant factor 
    `\Lambda_0` represents the total amplitude of all unmodulated contributions. The function `B(\lambda,t)` represents
    the background function due to the intermolecular contributions. The dipolar signal due to two-spin contributions given a
    distance distribution `P(r)` is given by the linear Fredholm integral 

    .. math:: 

        V(\mathbf{t}) = \int_0^\infty \mathrm{d}r P(r) K(\mathbf{t},r)

    The function can also compute three-spin interaction dipolar kernels that account for the additional the structure and shape
    of the echo modulation in multi-spin systems [2]_. The three-spin dipolar kernel is constructed from the 
    contributions arising from the combinations of different pair dipolar pathways [2]_ during a pulse sequence. 
    The three-spin dipolar kernel reads:

    .. math::

        K^{(3)}(\mathbf{t},\mathbf{r}) = & \Lambda_0 + \sum_k \lambda_k \frac{1}{2\pi}
        \int_0^{2\pi}\mathrm{d}\varphi 
        \int_0^{\pi}\mathrm{d}\theta_1
        \sin(\theta_1)
        \exp\Bigg(-i \sum_q^3 \frac{D \delta_{k,q}(\mathbf{t} - \mathbf{t}_{\mathrm{ref},k,q})}{r_q^3}
        (1 - 3\cos(\theta_{1q})^2 \Bigg) 
        \prod_k\prod_q^3 B(\lambda_k, \mathbf{\delta}_{k,q} (\mathbf{t} - \mathbf{t}_{\mathrm{ref},k,q}))

    with `\cos(\theta_{1q}) = \sin\chi_{1q} \cos\varphi \sin\theta_1 + \cos\chi_{1q} \cos\theta_1`,
    where `\mathbf{r} = (r_1,r_2,r_3)` represents the vector of the three interspin distances between the three interacting spins, 
    and `\chi_{1q}` are the angles between the first and `q`-th interspin vector. The three-spin dipolar pathways are now characterized 
    by an amplitude `\lambda_k`, and by a set of three refocusing times `\mathbf{t}_{\mathrm{ref},k,q}` and three (sub)harmonic factors
    `\mathbf{\delta}_{k,q}`.   

    Due to the higher dimensionality, the three-spin kernel is not constructed as a full matrix but folds the three distance dimensions 
    into a single one, such that the contribution to the dipolar echo modulation is given directly by

    .. math:: 

        V^{(3)}(\mathbf{t}) = \int_0^\infty\int_0^\infty \int_0^\infty  \mathrm{d}\mathbf{r} K^{(3)}(\mathbf{t},\mathbf{r})

    without the specification of a distance distribution.

    Parameters
    ----------
    t : ndarray or list thereof
        Experimental time incrementation vectors (time coordinates) `t`, in microseconds. For multi-dimensional experiments, 
        specify a list ``[t1,...,tD]`` of the time coordinates `\mathbf{t}` along each of the ``D`` different experimental dimensions. 
    
    r : ndarray or list thereof
        Distance vector `r`, in nanometers. For three-spin interactions, a list ``[r1,r2,r3]`` of three vectors of interspin distances 
        must be specified.
    
    pathways : list of dictionaries  or ``None``, optional
        List of dipolar pathways. 
        
        For *two-spin dipolar interactions*, each pathway is defined as a dictionary containing the pathway's amplitude :math:`\lambda_k` (``'amp'``), 
        refocusing time :math:`t_{\mathrm{ref},k}` (``'reftime'``) in microseconds, and (sub)harmonic :math:`\delta_k` (``'harmonic'``).
        For example, for a one-dimensional experiment ``{'amp:'lambda, 'reftime':tref, 'harmonic':delta}`` or ``{'amp:'lambda, 'reftime':tref}``. If ``'harmonic'`` is not
        specified, it is assumed to be 1 (modulated). For an unmodulated pathway, specify only the amplitude, i.e. ``{'amp':Lambda0}``.
        
        For *multi-dimensional experiments*, the pathway refocusing times \mathbf{\t}_{\mathrm{ref},k}` and harmonics `\mathbf{\delta}_k` 
        along the different dimensions must be specified as a list, e.g. for a pathway in a two-dimensional experiment 
        ``{'amp:'lambda, 'reftime':[tref1,tref2], 'harmonic':[delta1,delta2]}``
        
        For *three-spin interactions*, each pathway is defined as a dictionary containing the pathway's amplitude :math:`\lambda_k` (``'amp'``),
        a tuple of refocusing times :math:`t_{\mathrm{ref},k,q}` (``'reftime'``), and a tuple of (sub)harmonics :math:`\delta_{k,q}` (``'harmonic'``). 
        For example, for a one-dimensional experiment ``{'amp:'lambda, 'reftime':(tref1,tref2,tref3), 'harmonic':(delta1,delta2,delta3)}``. 
        
        If neither ``pathways`` or ``mod`` are specified (or set to ``None``), ``pathways=[{'amp':1,'reftime':0}]`` is used as the default.
    
    mod : scalar  or ``None``, optional
        Modulation depth for single-pathway dipolar models. This definition is equivalent to ``pathways=[{'amp':1-mod},{'amp':mod,'reftime':0}]``.
        If neither ``pathways`` or ``mod`` are specified (or ``None``), ``mod=1`` is used as the default. 
        
    bg : callable or array_like or ``None``, optional
        Dipolar background basis function. Must be a callable function, accepting a time-vector array as first input and a pathway amplitude 
        as a second, i.e. ``bg = lambda t,lam: bg_model(t,par,lam)``. If set to ``None``, no background function is included.
        For a single-pathway model, the numerical background decay can be passed as an array. 

    method : string, optional
        Numerical method for kernel matrix calculation: 

            * ``'fresnel'`` - uses Fresnel integrals for the kernel (fast, accurate)
            * ``'integral'`` - uses explicit integration function (slow, accurate)
            * ``'grid'`` - powder average via explicit grid integration (slow, inaccurate)
        
        The default is ``'fresnel'``. For all three-spin dipolar kernels, the ``'grid'`` method is used automatically. 

    orisel : callable  or ``None``, optional 
        Probability distribution of possible orientations of the interspin vector to account for orientation selection. Must be 
        a function taking a value of the angle θ∈[0,π/2] between the interspin vector and the external magnetic field and returning
        the corresponding probability density. If specified as ``None`` (by default), a uniform distribution is assumed. 
        Requires the ``'grid'`` or ``'integral'`` methods.

    excbandwidth : scalar, optional
        Excitation bandwidth of the pulses in MHz to account for limited excitation bandwidth.
        If not specified, an infinite excitation bandwidth is used.
        Requires the ``'grid'`` or ``'integral'`` methods.

    g : scalar, 2-element array, optional
        Electron g-values of the spin centers ``[g1, g2]``. If a single g is specified, ``[g, g]`` is used.
        If not specified, g = 2.002319... is used for both spins.

    complex : boolean, optional 
        If enabled, return the complex-valued dipolar kernel, otherwise return the real-valued dipolar kernel.
        Disabled by default. Requires the ``'fresnel'`` or ``'grid'`` methods.

    integralop : boolean, optional
        Return the kernel as an integral operator (i.e ``K = K*dr``) or not (``K``). Usage as an integral operator means that the 
        matrix operation ``V=K@P`` with a normalized distance distribution (i.e. ``trapz(r,P)==1``) leads to a signal ``V`` with 
        amplitude ``V(t=0)=1``. The default is ``True``.
    
    gridsize : scalar, optional
        Number of points on the grid of powder orientations to be used in the ``'grid'`` kernel calculation method. By default set to 1000 points.
        If the expected memory costs of the combined sizes of the ``t`` vector(s), ``r`` vector(s) and grid exceed the ``memorylimit`` value, 
        a ``MemoryError`` will be raised. 

    clearcache : boolean, optional
        Clear the cached dipolar kernels at the beginning of the function to save memory. Disabled by default.

    memorylimit : scalar, optional
        Memory limit to be allocated for the dipolar kernel. If the requested kernel exceeds this limit, the 
        execution will be stopped. The default is 8GB.  

    tinterp : array_like, optional 
        Effective dipolar evolution time vector used for interpolation. If specified, the elementary dipolar kernel is computed based on ``tinterp`` 
        and all contributions from the different pathways are interpolated based on it. This results in a significant reduction of the computation time. 
        The vector ``tinterp`` must cover all possible time shifts of the original ``t`` by the different refocusing times of the dipolar pathways.     

    Returns
    -------
    K : ndarray
        Dipolar kernel matrix.

    Notes
    -----
    For a multi-pathway DEER [1]_, [2]_, [3]_ signal (e.g, 4-pulse DEER with 2+1 contribution, 5-pulse DEER with 
    4-pulse DEER residual signal, and more complicated experiments), ``pathways`` contains a list of pathway amplitudes and refocusing times (in microseconds).
    The background function specified as ``bg`` is used as basis function, and the actual multipathway background 
    included into the kernel is computed using :ref:`dipolarbackground`. The background in included in the dipolar 
    kernel definition [3]_. Optionally, the harmonic (1 = fundamental, 2 = first harmonic, etc.) can be specified for modeling RIDME signals or
    subharmonics (1/2 = first subharmonic) can be specified for modeling SIFTER signals. 
    

    Examples
    --------
    To specify a model with an unmodulated offset and a single dipolar pathway that refocuses at time `t=0`, use::

        lam = 0.4  # modulation depth
        pathways = [
            {'amp': [1-lam]},
            {'amp': lam, 'reftime': 0}
            ]
        
        K = dl.dipolarkernel(t,r,pathways=pathways)


    A shorthand input syntax equivalent to this input::

        lam = 0.4 # modulation
        K = dl.dipolarkernel(t,r,mod=lam)


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
            {'amp': 1-lam1-lam2-lam3},
            {'amp': lam1, 'reftime': tref1},
            {'amp': lam2, 'reftime': tref2},
            {'amp': lam3, 'reftime': tref3},
        ]
        
        K = dl.dipolarkernel(t,r,pathways=pathways)


    To specify a dipolar kernel for an experiment on a three-spin system based on a single-pathway model use::

        tau1 = 0.5  # First interpulse delay (µs)         
        lam = 0.40  # Pathway #1 amplitude
        
        # Three dipolar pathways needed, one for each dipolar interaction in the three-spin system 
        pathways = [
            {'amp': 1-3*lam},
            {'amp': lam, 'reftime': (tau1,None,None), 'harmonic': (1,0,0)},
            {'amp': lam, 'reftime': (None,tau1,None), 'harmonic': (0,1,0)},
            {'amp': lam, 'reftime': (None,None,tau1), 'harmonic': (0,0,1)},
        ]
        
        K = dl.dipolarkernel(t, [r1,r2,r3], pathways=pathways)


    To construct a kernel for a three-spin system studied with a two-dimensional experiment, with a single dipolar pathway refocusing at, e.g.,  `t=(\tau_1,\tau_2)` use:: 

        # Three dipolar pathways needed, one for each dipolar interaction in the three-spin system 
        pathways = [
            {'amp': 1-3*lam},
            {'amp': lam, 'reftime': ([tau1,tau2],None,None), 'harmonic': ([1,1],0,0)},
            {'amp': lam, 'reftime': (None,[tau1,tau2],None), 'harmonic': (0,[1,1],0)},
            {'amp': lam, 'reftime': (None,None,[tau1,tau2]), 'harmonic': (0,0,[1,1])},
        ]
        
        K = dl.dipolarkernel(t, [r1,r2,r3], pathways=pathways)


    References
    ----------
    .. [1] L. Fábregas Ibáñez, M. H. Tessmer, G. Jeschke, and S. Stoll. 
        Dipolar pathways in dipolar EPR spectroscopy, Phys. Chem. Chem. Phys., 24 2022, 2504-2520

    .. [2] L. Fábregas Ibáñez, M. H. Tessmer, G. Jeschke, and S. Stoll. 
        Dipolar pathways in multi-spin and multi-dimensional dipolar EPR spectroscopy, Phys. Chem. Chem. Phys., 24 2022, 22645-22660

    .. [3] L. Fábregas Ibáñez, G. Jeschke, and S. Stoll. 
        DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson., 1, 209–224, 2020 

    .. [4] L. Fábregas Ibáñez, and G. Jeschke
        Optimal background treatment in dipolar spectroscopy, Physical Chemistry Chemical Physics, 22, 1855–1868, 2020.
    """
    # Clear cache of memoized function is requested
    if clearcache:
        elementarykernel_twospin.cache_clear()

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
    memory_requirement *= np.prod([len(t_) for t_ in t])*len(r[0])
    if method=='grid': 
        memory_requirement = memory_requirement*gridsize
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
        elif hasattr(pathway['amp'], "__len__"): 
            if len(np.atleast_1d(pathway['amp']))>1: 
                raise SyntaxError('A dipolar patwhay can only have one amplitude value.')
            else: pathway['amp'] = pathway['amp'].item()
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

    if tinterp is not None:
        # Construct interpolator, this way elementarykernel_twospin is executed only once independently of how many pathways there are
        Kinterpolator = [interp1d(tinterp,elementarykernel_twospin(tinterp,r_q,method,excbandwidth,gridsize,g,orisel,complex),axis=0,kind='cubic') for r_q in r]
        withinInterpolation = lambda tdip: np.all((np.max(tinterp) >= np.max(tdip)) & (np.min(tinterp) <= np.min(tdip)))

    # Define kernel matrix auxiliary functions
    def K0_2spin(tdip,q):
    #------------------------------------------------------------------------
        if tinterp is not None:
            if withinInterpolation(tdip):
                # Interpolated elementary dipolar kernel
                K0 = Kinterpolator[q](tdip)
                return K0
        # Exact elementary dipolar kernel construction
        K0 = elementarykernel_twospin(tdip,r[q],method,excbandwidth,gridsize,g,orisel,complex)
        return K0
    #------------------------------------------------------------------------
    
    def K0_3spin(tdip):
    #------------------------------------------------------------------------
        return elementarykernel_threespin(tdip,r,'grid',excbandwidth,gridsize,g,orisel,complex)
    #------------------------------------------------------------------------

    # Build dipolar kernel matrix, summing over all pathways
    K = Λ0
    for pathway in pathways:
        λ,tref,δ = [pathway['amp'],pathway['reftime'],pathway['harmonic']]

        # Set trefs of unmodulated patwhays to an arbitrary numerical value
        tref = [[trefq[d] if δqd!=0 else 0 for d,δqd in enumerate(δq)] for δq,trefq in zip(δ,tref)]

        # Determine number of spins participating in the pathway
        Nspin = np.sum([np.any(δq!=0) for δq in δ])+1 

        # Two-spin dipolar interactions
        if Nspin==2:

            # Determine the modulated dipolar interaction
            q = int(np.where([np.any(δq!=0) for δq in δ])[0])
            # Construc the effective dipolar evolution time for the modulated dipolar interactions
            tdip = np.sum(np.array([δ_qd*(t_d-tref_qd) for t_d,δ_qd,tref_qd in zip(t,δ[q],tref[q])],dtype=object), axis=0).astype(float)
            # Compute and accumulate the two-spin dipolar pathway contribution to the dipolar kernel
            K += λ*K0_2spin(tdip,q)

        elif Nspin==3:
            # Construct effective dipolar evolution time for all dipolar interactions
            tdip = [np.sum(np.array([δ_qd*(t_d-tref_qd) for t_d,δ_qd,tref_qd in zip(t,δq,trefq)],dtype=object), axis=0).astype(float) for δq,trefq in zip(δ,tref)]
            K += λ*K0_3spin(tdip)

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



#==============================================================================
#   TWO-SPIN ELEMENTARY DIPOLAR KERNEL
#==============================================================================

def _echomodulation(ωr,tdip,cosθ,ωex,complex):
#==============================================================================
    "Spin echo modulation due to dipolar coupling(s)"
    # Differentiate between complex and real valued cases for a speedup
    if complex: 
        echomod = np.exp(-1j*ωr*tdip*(1 - 3*cosθ**2))
    else:
        echomod = np.cos(ωr*abs(tdip)*(1 - 3*cosθ**2))
    # If given, include limited excitation bandwidth
    if not np.isinf(ωex):
        echomod *= np.exp(-(ωr*(1 - 3*cosθ**2))**2/ωex**2)
    return echomod
#==============================================================================

@cached(max_size=100)
def elementarykernel_twospin(tdip,r,method,ωex,gridsize,g,Pθ,complex):
#==============================================================================
    "Calculates the elementary two-spin dipolar interaction kernel (cached for speed)"

    # Work with vectors
    r = np.atleast_1d(r)

    # Dipolar frequencies
    ωr = (μ0/2)*μB**2*g[0]*g[1]/h*1e21/(r**3)  # rad μs^-1

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
    #------------------------------------------------------------------------
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
    #------------------------------------------------------------------------

    def elementarykernel_grid(tdip,ωex,gridsize,Pθ):
    #------------------------------------------------------------------------
        """Calculate kernel using grid-based powder integration (converges very slowly with gridsize)"""

        # Orientational grid
        cosθ = np.linspace(0,1,int(gridsize))

        D = len(tdip.shape)
        dims = np.arange(D + 2)

        # Expand dimensions of all variables on a multi-dimensional grid
        tdip_ =  np.expand_dims(tdip, axis=tuple(np.delete(dims, range(D)))) 
        ωr_   =  np.expand_dims(ωr,   axis=tuple(np.delete(dims, D)))
        cosθ_ =  np.expand_dims(cosθ, axis=tuple(np.delete(dims, D + 1)))

        # Integrate over the grid of echo modulations
        if orientationselection:
            # With orientation distribution Pθ
            K0 = np.dot(_echomodulation(ωr_,tdip_,cosθ_,ωex,complex),Pθ(np.arccos(cosθ)))/gridsize
        else: 
            # Without orientation selection
            K0 = np.sum(_echomodulation(ωr_,tdip_,cosθ_,ωex,complex),axis=-1)/gridsize
        return K0
    #------------------------------------------------------------------------

    def elementarykernel_integral(tdip,ωex,Pθ):
    #------------------------------------------------------------------------
        """Calculate kernel using explicit numerical integration """
        K0 = np.zeros((np.size(tdip),np.size(ωr)))
        for n,ωr_ in enumerate(ωr):
            def integrand(cosθ):
            #------------------------------------------------------------------------
                integ = _echomodulation(ωr_,tdip,cosθ,ωex,complex)
                # If given, include orientation selection
                if orientationselection:
                    integ = integ*Pθ(np.arccos(cosθ))  
                return integ
            #------------------------------------------------------------------------
            K0[:,n] = quad_vec(integrand,0,1,limit=1000)[0].squeeze()
        return K0
    #------------------------------------------------------------------------

    if method=='fresnel':
        K0 = elementarykernel_fresnel(tdip)
    elif method=='integral': 
        K0 = elementarykernel_integral(tdip,ωex,Pθ)
    elif method=='grid': 
        K0 = elementarykernel_grid(tdip,ωex,gridsize,Pθ)
    return K0
#==============================================================================

#==============================================================================
#   THREE-SPIN ELEMENTARY DIPOLAR KERNEL
#==============================================================================


def _echomodulation_threespin(tdips,ωrs,rs,θ1,φ1,weights,complex):
#==============================================================================
    "Spin echo modulation due to dipolar coupling(s)"

    # Orientation of the first interspin vector relative to the magnetic field
    cosθ1 = np.cos(θ1)
    sinθ1 = np.sin(θ1)
    cosφ1 = np.cos(φ1)

    # Law of cosines
    r1,r2,r3 = rs
    Δθ12 = np.arccos((r1**2 + r2**2 - r3**2)/(2*r1*r2))
    Δθ13 = np.arccos((r1**2 + r3**2 - r2**2)/(2*r1*r3))
    sinΔθ12, cosΔθ12 = np.sin(Δθ12), np.cos(Δθ12)
    sinΔθ13, cosΔθ13 = np.sin(Δθ13), np.cos(Δθ13)

    # Compute orientations of the other two interspin vectors relative to the magnetic field
    cosθ2 = sinΔθ12*cosφ1*sinθ1 + cosΔθ12*cosθ1
    cosθ3 = sinΔθ13*cosφ1*sinθ1 + cosΔθ13*cosθ1 
    
    # Precompute partial broadcasting for efficiency
    ωdip1_tdip1,ωdip2_tdip2,ωdip3_tdip3 = [ωdip_q*tdip_q for ωdip_q,tdip_q in zip(ωrs,tdips)]

    # Differentiate between complex and real valued cases for a speedup
    if complex: 
        echomod = ne.evaluate('exp(-1j*ωdip1_tdip1*(1 - 3*cosθ1**2) )* \
                               exp(-1j*ωdip2_tdip2*(1 - 3*cosθ2**2) )* \
                               exp(-1j*ωdip3_tdip3*(1 - 3*cosθ3**2) )')
    else:
        echomod = ne.evaluate('cos(ωdip1_tdip1*(1 - 3*cosθ1**2) )* \
                               cos(ωdip2_tdip2*(1 - 3*cosθ2**2) )* \
                               cos(ωdip3_tdip3*(1 - 3*cosθ3**2) )')
    # Orientation weighting
    echomod = echomod*weights

    return echomod
#==============================================================================

def elementarykernel_threespin(tdips,rs,method,ωex,gridsize,g,Pθ,complex):
#==============================================================================
    "Calculates the elementary three-spin dipolar interaction kernel"

    # Work with vectors
    rs = np.atleast_1d(*rs)

    # Dipolar frequencies
    ωrs = [(μ0/2)*μB**2*g[0]*g[1]/h*1e21/(r_q**3) for r_q in rs] # rad μs^-1

    # Orientation selection
    orientationselection = Pθ is not None 
    if orientationselection:
        # Ensure zero-derivatives at [0,π/2]
        θ = np.linspace(0,π/2,50) # rad
        Pθ_ = make_interp_spline(θ, Pθ(θ),bc_type="clamped")
        # Ensure normalization of probability density function (∫P(cosθ)dcosθ=1)
        Pθnorm,_ = quad(lambda cosθ: Pθ_(np.arccos(cosθ)),0,1,limit=1000)
        Pθ = lambda θ: Pθ_(θ)/Pθnorm


    def elementarykernel_threespin_grid(tsdip,gridsize,Pθ):
    #------------------------------------------------------------------------
        """Calculate kernel using grid-based powder integration (converges very slowly with gridsize)"""

        # Orientational grid
        φ,θ,weights = sophegrid(octants=4,maxphi=2*np.pi,size=int(np.floor(1/2*(np.sqrt(2*gridsize - 1) +1 ))))

        D = len(tsdip[0].shape)
        dims = np.arange(D + 2)
        # Expand dimensions of all variables on a multi-dimensional grid
        tsdip_ = [np.expand_dims(tdip, axis=tuple(np.delete(dims, range(D)))) for tdip in tsdip]
        ωr_   =  [np.expand_dims(ωr_q,   axis=tuple(np.delete(dims, D))) for ωr_q in ωrs]
        rs_   =  [np.expand_dims(r_q,   axis=tuple(np.delete(dims, D))) for r_q in rs]
        θ_ =  np.expand_dims(θ, axis=tuple(np.delete(dims, D + 1)))
        φ_ =  np.expand_dims(φ, axis=tuple(np.delete(dims, D + 1)))
        weights_ =  np.expand_dims(weights, axis=tuple(np.delete(dims, D + 1)))

        # Integrate over the grid of echo modulations
        if orientationselection:
            # With orientation distribution Pθ
            K0 = np.dot(_echomodulation_threespin(tsdip_,ωr_,rs_,θ_,φ_,weights_,complex),Pθ(θ))
        else:
            # Without orientation selection
            K0 = np.sum(_echomodulation_threespin(tsdip_,ωr_,rs_,θ_,φ_,weights_,complex),axis=-1)
        return K0
    #------------------------------------------------------------------------


    if  method=='grid': 
        K0 = elementarykernel_threespin_grid(tdips,gridsize,Pθ)
    return K0
#==============================================================================

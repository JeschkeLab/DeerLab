# dipolarkernel.py - Dipolar kernel operator
# -----------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

# Numpy + SciPy
import numpy as np
from numpy import pi,inf
import scipy.integrate
from scipy.special import fresnel
# Other packages
import types
import warnings
from memoization import cached
# DeerLab dependencies
from deerlab.dipolarbackground import dipolarbackground

# Fundamental constants (CODATA 2018)
ge = 2.00231930436256 # free-electron g factor
muB = 9.2740100783e-24 # Bohr magneton, J/T 
mu0 = 1.25663706212e-6 # magnetic constant, N A^-2 = T^2 m^3 J^-1 
h = 6.62607015e-34 # Planck constant, J/Hz
def w0(g):
    return (mu0/2)*muB**2*g[0]*g[1]/h*1e21 # Hz m^3 -> MHz nm^3 -> Mrad s^-1 nm^3

def dipolarkernel(t, r, pathways = 1, B = 1, method = 'fresnel', excbandwidth = inf, g = [ge, ge], 
                  integralop = True, nKnots = 5001, renormalize = True, renormpaths=True, clearcache = False):
#===================================================================================================
    r"""Compute the dipolar kernel operator which enables the linear transformation from
    distance-domain to time-domain data. 

    Parameters
    ----------
    t : array_like
        Dipolar time axis, in microseconds.
    
    r : array_like
        Distance axis, in nanometers.
    
    pathways : list of lists or scalar
        List of pathways. Each pathway is defined as a list of the pathway's amplitude (lambda), refocusing time (T0), 
        and harmonic (n), i.e. ``[lambda, T0, n]`` or ``[lambda, T0]`` for one pathway. If n is not given it is assumed to be 1. 
        For a pathway with unmodulated contribution, only the amplitude must be specified, i.e. ``[Lambda0]``.
        If a single value is specified, it is interpreted as the 4-pulse DEER pathway amplitude (modulation depth).  
    
    B : callable or array_like
        For a single-pathway model, the numerical background decay can be passed as an array. 
        For multiple pathways, a callable function must be passed, accepting a time-axis array as first input and a pathway amplitude as a second, i.e. ``B = lambda t,lam: bg_model(t,par,lam)``

    method : string, optional
        Numerical method for kernel matrix calculation: 

            * ``'fresnel'`` - uses Fresnel integrals for the kernel (default)
            * ``'integral'`` - uses explicit integration function (slow, accurate)
            * ``'grid'`` - powder average via explicit grid integration (slow, inaccurate)
        The default is ``'fresnel'``.

    excbandwidth : scalar, optional
        Excitation bandwidth of the pulses in MHz to account for limited excitation bandwidth [4]_.
    
    g : scalar, 2-element array, optional
        Electron g-values of the spin centers ``[g1, g2]``. If a single g is specified, ``[g, g]`` is assumed 
    
    integralop : boolean, optional
        Whether to return K as an integral operator (i.e ``K = K*dr``) or not (``K``). Usage as an integral operator means that the 
        matrix operation ``V=K@P`` with a normalized distance distribution (i.e. ``trapz(r,P)==1``) leads to a signal ``V`` with 
        ampliude ``V(t=0)=1``. Enabled by default.
    
    nKnots : scalar, optional
        Number of knots for the grid of powder orientations to be used in the ``'grid'`` kernel calculation method.
    
    renormalize : boolean, optional
        Re-normalization of multi-pathway kernels to ensure the equality ``K(t=0,r)==1`` is satisfied. Enabled by default.
    
    renormpaths: boolean, optional
        Normalization of the pathway amplitudes such that ``Lam0 + lam1 + ... + lamN = 1``. Enabled by default. 
    
    clearcache : boolean, optional
        Clear the cached dipolar kernels at the beginning of the function. Disabled by default.

    Returns
    --------
    K : ndarray
        Dipolar kernel operator, such that for a distance distribution (P), the dipolar signal is ``V = K@P``

    Notes
    -----
    For a multi-pathway DEER [1]_ signal (e.g, 4-pulse DEER with 2+1 contribution 5-pulse DEER with 4-pulse DEER residual signal, and more complicated experiments), ``pathways`` contains a list of modulation depths (amplitudes) and refocusing times (in microseconds).
    The background function specified as ''B'' is used as basis function, and the actual multipathway background included into the kernel is compued using :ref:`dipolarbackground`. The background in included in the dipolar kernel definition [2]_. 
    Optionally, the harmonic (1 = fundamental, 2 = first overtone, etc.) can be given as a third value in each row. This can be useful for modeling RIDME signals [3]_. If not given, the harmonic is 1 for all pathways. 


    Examples
    --------
    To specify the standard model for 4-pulse DEER with an unmodulated offset and a single dipolar pathway that refocuses at time 0, use::
    
        lam = 0.4  # modulation depth main signal
        pathways = [[1-lam], [lam, 0]]
        
        K = dl.dipolarkernel(t,r,pathways)


    A shorthand input syntax equivalent to this input::

        lam = 0.4
        K = dl.dipolarkernel(t,r,lam)


	To specify a more complete 4-pulse DEER model that, e.g includes the 2+1 contribution, use::

        Lam0 = 0.5  # unmodulated part
        lam = 0.4  # modulation depth main signal
        lam21 = 0.1  # modulation depth 2+1 contribution
        tau2 = 4  # refocusing time (us) of 2+1 contribution

        path0 = Lam0  # unmodulated pathway
        path1 = [lam1, 0]  # main dipolar pathway, refocusing at time zero
        path1 = [lam2, tau2]  # 2+1 dipolar pathway, refocusing at time tau2
        pathways = [path0, path1, path2]  
        
        K = dl.dipolarkernel(t,r,pathways)


    References
    ----------
    .. [1] L. Fábregas Ibáñez, G. Jeschke, and S. Stoll. 
        DeerLab: A comprehensive toolbox for analyzing dipolar EPR spectroscopy data, Magn. Reson., 1, 209–224, 2020 

    .. [2] L. Fábregas Ibáñez, and G. Jeschke
        Optimal background treatment in dipolar spectroscopy, Physical Chemistry Chemical Physics, 22, 1855–1868, 2020.

    .. [3] K. Keller, V. Mertens, M. Qi, A. I. Nalepa, A. Godt, A. Savitsky, G. Jeschke, and M. Yulikov
        Computing distance distributions from dipolar evolution data with overtones: RIDME spectroscopy with Gd(III)-based spin labels, Physical Chemistry Chemical Physics, 19

    .. [4] J. E. Banham, C. M. Baker, S. Ceola, I. J. Day, G.H. Grant, E. J. J. Groenen, C. T. Rodgers, G. Jeschke, C. R. Timmel
        Distance measurements in the borderline region of applicability of CW EPR and DEER: A model study on a homologous series of spin-labelled peptides, Journal of Magnetic Resonance, 191, 2, 2008, 202-218
    
    """
    # Clear cache of memoized function is requested
    if clearcache:
        calckernelmatrix.cache_clear()

    # Ensure that all inputs are numpy arrays
    r = np.atleast_1d(r)
    t = np.atleast_1d(t)

    if type(B) is types.LambdaType:
        ismodelB = True
    else:
        ismodelB = False
        B = np.atleast_1d(B)

    g = np.atleast_1d(g)
    if len(g) == 1:
        g = [g, g]
    elif len(g) != 2:
        raise TypeError('The g-value must be specified as a scalar or a two-element array.')

    if np.any(r<=0):
        raise ValueError("All elements in r must be nonnegative and nonzero.")

    if not isinstance(pathways,list): pathways = [pathways] 
    if len(pathways) == 1:
        lam = pathways[0]
        pathways = [[1-lam], [lam, 0]]

    paths = [np.atleast_1d(path) for path in pathways]

    # Get unmodulated pathways    
    unmodulated = [paths.pop(i) for i,path in enumerate(paths) if len(path)==1]
    # Combine all unmodulated contributions
    Lambda0 = sum(np.concatenate([path for path in unmodulated]))

    # Check structure of pathways
    for i,path in enumerate(paths):
        if len(path) == 2:
            # If harmonic is not defined, append default n=1
            paths[i] = np.append(path,1) 
        elif len(path) != 3:
            # Otherwise paths are not correctly defined
            raise KeyError('The pathway #{} must be a list of two or three elements [lam, T0] or [lam, T0, n]'.format(i))

    # Normalize the pathway amplitudes to unity
    if renormpaths:
        lamsum = Lambda0 + sum([path[0] for path in paths])
        Lambda0 /= lamsum
        for i in range(len(paths)):
            paths[i][0] /= lamsum 

    # Define kernel matrix auxiliary function
    kernelmatrix = lambda t: calckernelmatrix(t,r,method,excbandwidth,nKnots,g)

    # Build dipolar kernel matrix, summing over all pathways
    K = Lambda0
    for pathway in paths:
        lam,T0,n = pathway
        K = K + lam*kernelmatrix(n*(t-T0))

    # Renormalize if requested
    if renormalize:
        Knorm = Lambda0
        for pathway in paths:
            lam,T0,n = pathway
            Knorm = Knorm + lam*kernelmatrix(-T0*n)
        K = K/Knorm

    # Multiply by background
    if ismodelB:
        B = dipolarbackground(t,pathways,B)
    K = K*B[:,np.newaxis]

    # Include delta-r factor for integration
    if integralop and len(r)>1:
            
        # Vectorized matrix form of trapezoidal integration method
        dr = np.zeros(len(r))
        dr[0] = r[1] - r[0]
        dr[-1] = r[-1] - r[-2]
        dr[1:-1] = r[2:] - r[0:-2]
        dr = dr/2
        K = K*dr
    
    return K
#==============================================================================

@cached
def calckernelmatrix(t,r,method,wex,nKnots,g):
#==============================================================================

    t = np.atleast_1d(t)
    r = np.atleast_1d(r)

    # Pre-allocate K matrix
    nr = np.size(r)
    nt = np.size(t)  
    K = np.zeros((nt,nr))
    wr = w0(g)/(r**3)  # rad s^-1

    def kernelmatrix_fresnel(t,r):
    #==========================================================================
        """Calculate kernel using Fresnel integrals (fast and accurate)"""

       # Calculation using Fresnel integrals
        ph = np.outer(np.abs(t), wr)
        kappa = np.sqrt(6*ph/pi)
        
        # Supress divide by 0 warning        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S, C = fresnel(kappa)/kappa
        
        K = C*np.cos(ph) + S*np.sin(ph)
        K[t==0] = 1 

        return K
    #==========================================================================
    
    def kernelmatrix_grid(t,r,wex,nKnots):
    #==========================================================================
        """Calculate kernel using grid-based powder integration (converges very slowly with nKnots)"""

        # Orientational grid
        costheta = np.linspace(0,1,int(nKnots))
        q = 1 - 3*costheta**2
        # Vectorized 3D-grid evaluation (t,r,powder averaging)
        C = np.cos(wr[np.newaxis,:,np.newaxis]*q[np.newaxis,np.newaxis,:]*abs(t[:,np.newaxis,np.newaxis]))
        # If given, include limited excitation bandwidth
        if not np.isinf(wex):
            C = C*np.exp(-(wr[np.newaxis,:,np.newaxis]*q[np.newaxis,np.newaxis,:])**2/wex**2)
        K = np.sum(C,2)/nKnots
        return K
    #==========================================================================

    def kernelmatrix_integral(t,r,wex):
    #==========================================================================
        """Calculate kernel using explicit numerical integration """
        for ir in range(len(wr)):
            for it in range(len(t)):
                #==================================================================
                def integrand(costheta):
                    integ = np.cos(wr[ir]*abs(t[it])*(1-3*costheta**2))
                    # If given, include limited excitation bandwidth
                    if not np.isinf(wex):
                        integ = integ*np.exp(-(wr[ir]*(1-3*costheta**2))**2/wex**2)
                    return integ
                #==================================================================   
                K[it,ir],_ = scipy.integrate.quad(integrand,0,1,limit=1000)
        return K
    #==========================================================================

    if method=='fresnel':
        K = kernelmatrix_fresnel(t,r)
    elif method=='integral': 
        K = kernelmatrix_integral(t,r,wex)
    elif method=='grid': 
        K = kernelmatrix_grid(t,r,wex,nKnots)
    return K
#==============================================================================
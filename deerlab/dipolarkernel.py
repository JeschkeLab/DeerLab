# dipolarkernel.py - Dipolar kernel operator
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

def dipolarkernel(t,r,pathinfo = 1, B = 1, method = 'fresnel', excbandwidth = inf, g = [ge, ge], 
                  integralop = True, nknots = 5001, renormalize = True, clearcache = False):
#===================================================================================================
    """
    Dipolar kernel operator
    =========================

    Computes the dipolar kernel operator which enables the linear transformation from
    distance-domain to time-domain data. 

    Usage: 
    -----------
        K = dipolarkernel(t,r)
        K = dipolarkernel(t,r,lambda)
        K = dipolarkernel(t,r,lambda,B)
        K = dipolarkernel(t,r,pathinfo)
        K = dipolarkernel(t,r,pathinfo,Bmodel)

    Arguments: 
    -----------
    t (N-element array)  
        Dipolar time axis, in microseconds
    r (M-element array) 
        Distance axis, in nanometers
    lambda (scalar)
        Modulation depth (value between 0 and 1) (default=1)
        This is equivalent to pathinfo = [1-lambda NaN 0; lambda 0 1]
    pathinfo (px2-array, px3-array)
        Array of pathway amplitudes (lambda), refocusing points (T0), and harmonics (n) 
        for multiple (p) dipolar pathways. Each row contains [lambda T0 n] or [lambda T0] 
        for one pathway. If n is not given it is assumed to be 1.
    B/Bmodel (lambda-function, N-element array)
        For a single-pathway model, the numerical background decay can be passed as an array. 
        For multiple pathways, a lambda-function must be specified with the following structure:
            Bmodel = lambda par,lambda: bg_model(t,par,lambda)
        By default, no background decay is included.
 
    Keyword arguments:
    ------------------
    method (string, optional)
        Numerical method for kernel matrix calculation:
                    'fresnel' - uses Fresnel integrals for the kernel (default)
                    'integral' - uses explicit integration function (slow, accurate)
                    'grid' - powder average via explicit grid integration (slow, inaccurate)
    excbandwidth (scalar, optional)
        Excitation bandwidth of the pulses in MHz to account for limited excitation bandwidth
    g (scalar, 2-element array, optional)
        Electron g-values of the spin centers [g1, g2]. If a single g is specified, [g, g] is assumed 
    integralop (boolean, optional)
        Whether to return K as an integral operator (i.e K = K*dr) or not (K). Default is True.
    nknots (scalar, optional)
        Number of knots for the grid of powder orientations to be used in the 'grid' kernel calculation method.
    renormalize (boolean, optional)
        Re-normalization of multi-pathway kernels to ensure the equality K(t=0,r)==1 is satisfied. Default is True.
    clearcache (boolean, optional)
        Clear the cached dipolar kernels at the beginning of the function. Default is False.

    Returns:
    -----------
    K (NxM-array)
        Dipolar kernel operator, such that for a distance distribution (P), the dipolar signal is V = K@P


    """

    # Clear cache of memoized function is requested
    if clearcache:
        calckernelmatrix.cache_clear()

    # Ensure that all inputs are numpy arrays
    r = np.atleast_1d(r)
    t = np.atleast_1d(t)
    pathinfo = np.atleast_1d(pathinfo)

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

    if not np.isreal(pathinfo).all:
        raise TypeError('lambda/pathinfo must be a numeric array.')

    if len(pathinfo) == 1:
        lam = pathinfo[0]
        pathinfo = np.array([[1-lam, np.NaN], [lam, 0]])

    if not np.any(np.shape(pathinfo)[1] != np.array([2, 3])):
        raise TypeError('pathinfo must be a numeric array with two or three columns.')

    if np.any(np.isnan(pathinfo[:,0])):
        raise ValueError('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path[1,:] = [Lam0 NaN]')

    # Normalize the pathway amplitudes to unity
    pathinfo[:,0] = pathinfo[:,0]/sum(pathinfo[:,0])
    lam = pathinfo[:,0]
    T0 = pathinfo[:,1]
    if np.shape(pathinfo)[1] == 2:
        n = np.ones(np.shape(T0))
    else:
        n = pathinfo[:,2]
    
    # Combine all unmodulated components into Lambda0, and eliminate from list
    unmodulated = np.where(np.isnan(T0))
    Lambda0 = np.sum(lam[unmodulated])
    lam = np.delete(lam, unmodulated)
    T0 = np.delete(T0, unmodulated)
    n = np.delete(n, unmodulated)
    nModPathways = len(lam)

    kernelmatrix = lambda t: calckernelmatrix(t,r,method,excbandwidth,nknots,g)
    
    # Build dipolar kernel matrix, summing over all pathways
    K = Lambda0
    for pathway in range(nModPathways):
        K = K + lam[pathway]*kernelmatrix(n[pathway]*(t-T0[pathway]))

    # Renormalize if requested
    if renormalize:
        Knorm = Lambda0
        for pathway in range(nModPathways):
            Knorm = Knorm + lam[pathway]*kernelmatrix(-T0[pathway]*n[pathway])
        K = K/Knorm


    # Multiply by background
    if ismodelB:
        B = dipolarbackground(t,pathinfo,B)
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
def calckernelmatrix(t,r,method,wex,nknots,g):
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
    
    def kernelmatrix_grid(t,r,wex,nknots):
    #==========================================================================
        """Calculate kernel using grid-based powder integration (converges very slowly with nknots)"""

        # Orientational grid
        costheta = np.linspace(0,1,int(nknots))
        q = 1 - 3*costheta**2
        # Vectorized 3D-grid evaluation (t,r,powder averaging)
        C = np.cos(wr[np.newaxis,:,np.newaxis]*q[np.newaxis,np.newaxis,:]*abs(t[:,np.newaxis,np.newaxis]))
        # If given, include limited excitation bandwidth
        if not np.isinf(wex):
            C = C*np.exp(-(wr[np.newaxis,:,np.newaxis]*q[np.newaxis,np.newaxis,:])**2/wex**2)
        K = np.sum(C,2)/nknots
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
        K = kernelmatrix_grid(t,r,wex,nknots)
    return K
#==============================================================================
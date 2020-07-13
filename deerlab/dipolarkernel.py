import numpy as np
from numpy import pi,inf
import types
from deerlab.dipolarbackground import dipolarbackground
import scipy.integrate
from scipy.special import fresnel
import warnings
from memoization import cached

# Fundamental constants (CODATA 2018)
ge = 2.00231930436256 # free-electron g factor
muB = 9.2740100783e-24 # Bohr magneton, J/T 
mu0 = 1.25663706212e-6 # magnetic constant, N A^-2 = T^2 m^3 J^-1 
h = 6.62607015e-34 # Planck constant, J/Hz
def w0(g):
    return (mu0/2)*muB**2*g[0]*g[1]/h*1e21 # Hz m^3 -> MHz nm^3 -> Mrad s^-1 nm^3

def dipolarkernel(t,r,
        pathinfo=1,
        B=1,
        method = 'fresnel',
        excbandwidth = inf,
        g = [ge, ge],
        nknots = 5001,
        renormalize = True,
        cache = True
    ):
    """
    K = dipolarkernel(t,r)
    Calculate dipolar kernel matrix.
    Assumes t in microseconds and r in nanometers
    """

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
    print(len(g))
    if len(g) == 1:
        g = [g, g]
    elif len(g) != 2:
        raise TypeError('The g-value must be specified as a scalar or a two-element array.')

    if np.any(r<0):
        raise ValueError("All elements in r must be nonnegative.")

    if not np.isreal(pathinfo).all:
        raise TypeError('lambda/pathinfo must be a numeric array.')

    if len(pathinfo) == 1:
        lam = pathinfo[0]
        pathinfo = np.array([[1-lam, np.NaN], [lam, 0]])

    if not np.any(np.shape(pathinfo)[1] != np.array([2, 3])):
        raise TypeError('pathinfo must be a numeric array with two or three columns.')

    if np.any(np.isnan(pathinfo[:, 0])):
        raise ValueError('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path[1,:] = [Lam0 NaN]')

    # Normalize the pathway amplitudes to unity
    pathinfo[:, 0] = pathinfo[:, 0] / sum(pathinfo[:, 0])
    lam = pathinfo[:, 0]
    T0 = pathinfo[:, 1]
    if np.shape(pathinfo)[1] == 2:
        n = np.ones(np.shape(T0))
    else:
        n = pathinfo[:, 1]
    

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
    if len(r)>1:
        dr = np.mean(np.diff(r))
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
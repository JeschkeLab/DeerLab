import math as m
import warnings
import types
import numpy as np
import pandas as pd
from dipolarbackground import dipolarbackground
from scipy.special import fresnel

# Fundamental constants (CODATA 2018)
ge = 2.00231930436256 # free-electron g factor
muB = 9.2740100783e-24 # Bohr magneton, J/T 
mu0 = 1.25663706212e-6 # magnetic constant, N A^-2 = T^2 m^3 J^-1 
h = 6.62607015e-34 # Planck constant, J/Hz
nu0 = (mu0/4/np.pi)*muB**2*ge*ge/h*1e21 # Hz m^3 -> MHz nm^3
w0 = 2*np.pi*nu0 # Mrad s^-1 nm^3

def dipolarkernel(t,r,pathinfo=1,B=1):
    """
    K = dipolarkernel(t,r)
    Calculate dipolar kernel matrix.
    Assumes t in microseconds and r in nanometers
    """
    
    Renormalize = True

    # Ensure that all inputs are numpy arrays
    r = np.atleast_1d(r)
    t = np.atleast_1d(t)
    pathinfo = np.atleast_1d(pathinfo)

    if type(B) is types.LambdaType:
        ismodelB = True
    else:
        ismodelB = False
        B = np.atleast_1d(B)


    if np.any(r<0):
        raise ValueError("All elements in r must be nonnegative.")

    if not np.isreal(pathinfo).all:
        raise TypeError('lambda/pathinfo must be a numeric array.')

    if len(pathinfo)==1:
        lam = pathinfo
        pathinfo = np.array([[1-lam, np.NaN], [lam, 0]], dtype=object)

    if not any(np.shape(pathinfo)[1]!=np.array([2, 3])):
        raise TypeError('pathinfo must be a numeric array with two or three columns.')

    if any(pd.isnull(pathinfo[:,0])):
        raise ValueError('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path[1,:] = [Lam0 NaN]')

    # Normalize the pathway amplitudes to unity
    pathinfo[:,0] = pathinfo[:,0]/sum(pathinfo[:,0])
    lam = pathinfo[:,0]
    T0 = pathinfo[:,1]
    if np.shape(pathinfo)[1]==2:
        n = np.ones(np.shape(T0))
    else:
        n = pathinfo[:,1]
    

    # Combine all unmodulated components into Lambda0, and eliminate from list
    unmodulated = np.where(pd.isnull(T0))
    Lambda0 = np.sum(lam[unmodulated])
    lam = np.delete(lam,unmodulated)
    T0 = np.delete(T0,unmodulated)
    n = np.delete(n,unmodulated)  

    nModPathways = len(lam)

    kernelmatrix = lambda t: calckernelmatrix('fresnel',t,r)

    # Build dipolar kernel matrix, summing over all pathways
    K = Lambda0
    for pathway in range(nModPathways):
        K = K + lam[pathway]*kernelmatrix(n[pathway]*(t-T0[pathway]))

    # Renormalize if requested
    if Renormalize:
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

def calckernelmatrix(method,t,r):

    # Pre-allocate K matrix
    nr = np.size(r)
    nt = np.size(t)  
    K = np.zeros((nt,nr))

    def kernelmatrix_fresnel(t,r):
        # Calculation using Fresnel integrals
        wr = w0/(r**3)  # rad s^-1
        ph = np.outer(np.abs(t), wr)
        kappa = np.sqrt(6 * ph / np.pi)
        
        # Supress divide by 0 warning        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            S, C = fresnel(kappa) / kappa
        
        K = C * np.cos(ph) + S * np.sin(ph)
        K[t==0] = 1 

        return K

    if method=='fresnel':
        K = kernelmatrix_fresnel(t,r)

    return K

import numpy as np
import math as m
from scipy.special import fresnel

# Fundamental constants (CODATA 2018)
ge = 2.00231930436256 # free-electron g factor
muB = 9.2740100783e-24 # Bohr magneton, J/T 
mu0 = 1.25663706212e-6 # magnetic constant, N A^-2 = T^2 m^3 J^-1 
h = 6.62607015e-34 # Planck constant, J/Hz
nu0 = (mu0/4/np.pi)*muB**2*ge*ge/h*1e21 # Hz m^3 -> MHz nm^3
w0 = 2*np.pi*nu0 # Mrad s^-1 nm^3

def dipolarkernel(t,r):
    """
    K = dipolarkernel(t,r)
    Calculate dipolar kernel matrix.
    Assumes t in microseconds and r in nanometers
    """
    
    r = np.atleast_1d(r)
    t = np.atleast_1d(t)


    if np.any(r<=0):
        raise ValueError("All elements in r must be nonnegative.")
    
    # Pre-allocate K matrix
    nr = np.size(r)
    nt = np.size(t)
    K = np.zeros((nt,nr))
    
    # Calculation using Fresnel integrals
    wr = w0/(r**3)  # rad s^-1
    for ir in range(nr):
        ph = wr[ir]*(np.abs(t))
        kappa = np.sqrt(6*ph/m.pi)
        S, C = fresnel(kappa)
        K[:,ir] = (np.cos(ph)*C+np.sin(ph)*S)/kappa
    
    K[t==0,:] = 1   # fix div by zero
    
    # Include delta-r factor for integration
    if len(r)>1:
        dr = np.mean(np.diff(r))
        K = K*dr
    
    return K

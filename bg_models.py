import numpy as np
import math as m
import matrixtools as mat

def bg_exp(t,param,lam=1):
    kappa = param
    B = np.exp(-lam*kappa*np.abs(t))
    B = mat.column(B)
    return B

def bg_hom3d(t,param,lam=1):
    conc = param            # concentration, uM
    
    NA = 6.02214076e23      # Avogadro constant, mol^-1
    muB = 9.2740100783e-24  # Bohr magneton, J/T (CODATA 2018 value);
    mu0 = 1.25663706212e-6  # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
    h = 6.62607015e-34      # Planck constant, J/Hz (CODATA 2018)
    ge = 2.00231930436256   # free-electron g factor (CODATA 2018 value)
    hbar = h/2/m.pi         # reduced Planck constant, J/(rad/s)

    D = (mu0/4/m.pi)*(muB*ge)**2/hbar   # dipolar constant, m^3 s^-1
    
    conc = conc*1e-6*1e3*NA # umol/L -> mol/L -> mol/m^3 -> spins/m^3
    
    # Compute background function
    B = np.exp(-8*m.pi**2/9/m.sqrt(3)*lam*conc*D*np.abs(t*1e-6))
    B = mat.column(B)
    return B

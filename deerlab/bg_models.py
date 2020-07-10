import numpy as np
import math as m

def bg_exp(*args):

    if not args:
        info = dict(
            Parameters = ('Decay Rate'),
            Units = ('us-1'),
            Start = np.asarray([0.35]),
            Lower = np.asarray([0]),
            Upper = np.asarray([200])
        )
        return info
    elif len(args)!=2 and len(args)!=3:
        raise TypeError('Two or three input arguments required: bg_model(t,p) or bg_model(r,p,lam)')
    else:
        t = args[0]
        param = args[1]

    if len(args)>2:
        lam = args[2]
    else:
        lam = 1    
    
    t = np.atleast_1d(t)
    param = np.atleast_1d(param)

    kappa = param[0]
    B = np.exp(-lam*kappa*np.abs(t))
    return B

def bg_hom3d(*args):

    if not args:
        info = dict(
            Parameters = ('Spin concentration'),
            Units = ('uM'),
            Start = np.asarray([50]),
            Lower = np.asarray([0.01]),
            Upper = np.asarray([5000])
        )
        return info
    elif len(args)!=2 and len(args)!=3:
        raise TypeError('Two or three input arguments required: bg_model(t,p) or bg_model(r,p,lam)')
    else:
        t = args[0]
        param = args[1]

    if len(args)>2:
        lam = args[2]
    else:
        lam = 1  

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
    return B

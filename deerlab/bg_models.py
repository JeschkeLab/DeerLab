# bg_models.py - Background parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.
 
import numpy as np
import math as m
import scipy as scp
from numpy import pi
import inspect
from deerlab.utils import load_exvolume_redfactor
from deerlab.model import Model

# Natural constants
Nav = 6.02214076e23      # Avogadro constant, mol^-1
μB = 9.2740100783e-24  # Bohr magneton, J/T (CODATA 2018 value)
μ0 = 1.25663706212e-6  # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
h = 6.62607015e-34      # Planck constant, J/Hz (CODATA 2018)
ge = 2.00231930436256   # free-electron g factor (CODATA 2018 value)
hbar = h/2/pi         # reduced Planck constant, J/(rad/s)

D = (μ0/4/pi)*(μB*ge)**2/hbar   # dipolar constant, m^3 s^-1


def _docstring(model,notes):
#---------------------------------------------------------------------------------------
    args = model._parameter_list(order='vector').tolist()
    args.insert(model.axis_argument[1],model.axis_argument[0])

    parameters = ''
    for arg in args:
        if arg==model.axis_argument[0]:
            type = 'array_like'
            parameters += f'\n    {arg} : {type} \n        Time vector, in microseconds.'
        elif len(np.atleast_1d(getattr(model,arg).idx))>1:
            type = 'array_like'
            parameters += f'\n    {arg} : {type} \n        {str(getattr(model,arg).description):s}'
        else: 
            type = 'scalar'
            parameters += f'\n    {arg} : {type} \n        {str(getattr(model,arg).description):s}'

        
    string = inspect.cleandoc(f"""

    {model.description}

    Parameters
    ----------
    {parameters}
    
    Returns
    -------

    B : ndarray
        Dipolar background vector.

    Notes
    -----

    **Parameter List**

    ============ ========= ========== =========== ========== ==========================
        Name       Lower     Upper      Type        Units     Description  
    ============ ========= ========== =========== ========== ==========================""")
    for n,paramname in enumerate(model._parameter_list(order='vector')): 
        string += f'\n   {paramname:7s}'
        string += f'     {getattr(model,paramname).lb:5.3g}'
        string += f'     {getattr(model,paramname).ub:5.3g}'
        string += f'      {"linear" if getattr(model,paramname).linear else "nonlin"}'
        string += f'       {str(getattr(model,paramname).units):6s}'
        string += f'   {str(getattr(model,paramname).description):s}'
    string += f'\n============ ========= ========== =========== ========== =========================='

    string += f'\n{notes}'

    model.__doc__ = string
    return model
#---------------------------------------------------------------------------------------

#=======================================================================================
#                                     bg_hom3d
#=======================================================================================
notes = r"""
**Model:**

This model describes the inter-molecular interaction of one observer spin with a 3D homogenous distribution of spins of concentration `c_s`

.. image:: ../images/model_scheme_bg_hom3d.png
   :width: 350px

The expression for this model is

.. math::
   B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\lambda c_s D |t|\right)`

where `c_s` is the spin concentration (entered in spins/m\ :sup:`3` into this expression) and D is the dipolar constant

.. math::
   D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}
"""  
def hom3d_fcn(t,conc,lam):
#---------------------------------------------------------------------------------------
    # Units conversion    
    conc = conc*1e-6*1e3*Nav # umol/L -> mol/L -> mol/m^3 -> spins/m^3
    # Compute background function
    B = np.exp(-8*pi**2/9/m.sqrt(3)*lam*conc*D*np.abs(t*1e-6))
    return B
#---------------------------------------------------------------------------------------
# Create model
bg_hom3d = Model(hom3d_fcn,axis='t')
bg_hom3d.description = 'Background from homogeneous distribution of spins in a 3D medium'
# Parameters
bg_hom3d.conc.set(description='Spin concentration', lb=0.01, ub=5000, par0=50, units='μM')
bg_hom3d.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, units='')
# Add documentation
bg_hom3d.__doc__ = _docstring(bg_hom3d,notes)


#=======================================================================================
#                                     bg_hom3dex
#=======================================================================================
notes = r"""
**Model:**

.. image:: ../images/model_scheme_bg_hom3dex.png
   :width: 350px

This implements a hard-shell excluded-volume model, with spin concentration ``c_s`` (first parameter, in μM) and distance of closest approach ``R`` (second parameter, in nm).

The expression for this model is

.. math:: B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\alpha(R) \lambda c_s D |t|\right)`

where :math:`c_s` is the spin concentration (entered in spins/m\ :sup:`3` into this expression) and :math:`D` is the dipolar constant

.. math:: D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

The function :math:`\alpha(R)` of the exclusion distance :math:`R` captures the excluded-volume effect. It is a smooth function, but doesn't have an analytical representation. For details, see `Kattnig et al, J.Phys.Chem.B 2013, 117, 16542 <https://pubs.acs.org/doi/abs/10.1021/jp408338q>`_.
"""  
def hom3dex(t,conc,rexcl,lam):    
    # Load precalculated reduction factor look-up table (Kattnig Eq.(18))
    dR_tab,alphas_tab = load_exvolume_redfactor()

    # Conversion: umol/L -> mol/L -> mol/m^3 -> spins/m^3
    conc = conc*1e-6*1e3*Nav 

    A = (μ0/4/pi)*(ge*μB)**2/hbar # Eq.(6) m^3 s^-1

    # Calculate reduction factor (Eq.(18))
    if rexcl==0:
        alpha = 1
    else:
        dR = A*abs(t*1e-6)/(rexcl*1e-9)**3 # unitless
        
        # Use interpolation of look-up table for small dR
        small = dR < max(dR_tab)
        alpha = np.zeros(np.shape(dR))
        alpha[small] = np.interp(dR[small], dR_tab, alphas_tab)
        
        # For large dR, use limiting dR->inf expression
        alpha[~small] = 1 - (3/2/pi)*np.sqrt(3)/dR[~small]

    K = 8*pi**2/9/np.sqrt(3)*A*abs(t*1e-6)*alpha # Eq.(17)
    B = np.exp(-lam*conc*K) # Eq.(13)
    return B
# Create model
bg_hom3dex = Model(hom3dex,axis='t')
bg_hom3dex.description = 'Background from homogeneous distribution of spins with excluded-volume effects'
# Parameters
bg_hom3dex.conc.set(description='Spin concentration', lb=0.01, ub=5000, par0=50, units='μM')
bg_hom3dex.rexcl.set(description='Exclusion distance', lb=0.01, ub=20, par0=1, units='nm')
bg_hom3dex.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, units='')
# Add documentation
bg_hom3dex.__doc__ = _docstring(bg_hom3dex,notes)


#=======================================================================================
#                                     bg_homfractal
#=======================================================================================
notes = r"""
**Model:**

This implements the background due to a homogeneous distribution of spins in a d-dimensional space, with d-dimensional spin concentration ``c_d``.
"""  
def homfractal(t,fconc,fdim,lam):
    # Unpack model paramters
    d = float(fdim)     # fractal dimension    

    # Units conversion of concentration    
    conc = fconc*1e-6*(np.power(10,d))*Nav # umol/dm^d -> mol/m^d -> spins/m^d

    # Compute constants
    if d==3:
        c = -pi/2
        Λ = 4/3/np.sqrt(3)
    else:
        c = np.cos(d*pi/6)*scp.special.gamma(-d/3)
        integrand = lambda z: abs(1-3*z**2)**(d/3)
        Λ,_ = scp.integrate.quad(integrand,0,1,limit=1000)

    # Compute background function
    B = np.exp(4*pi/3*c*Λ*lam*conc*D**(d/3)*abs(t*1e-6)**(d/3))
    return B
 # ======================================================================
# Create model
bg_homfractal = Model(homfractal,axis='t')
bg_homfractal.description = 'Background from homogeneous distribution of spins in a fractal medium'
# Parameters
bg_homfractal.fconc.set(description='Fractal concentration of spins', lb=0.01, ub=5000, par0=50, units='μmol/dmᵈ')
bg_homfractal.fdim.set(description='Fractal dimensionality', lb=0.01, ub=5.99, par0=3, units='')
bg_homfractal.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, units='')
# Add documentation
bg_homfractal.__doc__ = _docstring(bg_homfractal,notes)


#=======================================================================================
#                                     bg_exp
#=======================================================================================
notes= r"""
**Model:**

.. math::

   B(t) = \exp\left(-\kappa \vert t \vert\right)

Although the ``bg_exp`` model has the same functional form as ``bg_hom3d``, it is distinct since its 
parameter is a decay rate constant and not a spin concentration like for ``bg_hom3d``. 
    """  
def exp(t,decay):
    return np.exp(-decay*np.abs(t))
# Create model
bg_exp = Model(exp,axis='t')
bg_exp.description = 'Exponential background model'
# Parameters
bg_exp.decay.set(description='Decay rate', lb=0, ub=200, par0=0.35, units='μs⁻¹')
# Add documentation
bg_exp.__doc__ = _docstring(bg_exp,notes)



#=======================================================================================
#                                     bg_strexp
#=======================================================================================
notes = r"""
**Model:**

.. math::

    B(t) = \exp\left(-\kappa \vert t\vert^{d}\right)

Although the ``bg_strexp`` model has the same functional form as ``bg_homfractal``, it is distinct since its
first parameter is a decay rate constant and not a spin concentration like for ``bg_homfractal``.
"""  
def strexp(t,decay,stretch):
    return np.exp(-decay*abs(t)**stretch)
# Create model
bg_strexp = Model(strexp,axis='t')
bg_strexp.description = 'Stretched exponential background model'
# Parameters
bg_strexp.decay.set(description='Decay rate', lb=0, ub=200, par0=0.25, units='μs⁻¹')
bg_strexp.stretch.set(description='Stretch factor', lb=0, ub=6, par0=1, units='')
# Add documentation
bg_strexp.__doc__ = _docstring(bg_strexp,notes)



#=======================================================================================
#                                     bg_prodstrexp
#=======================================================================================
notes = r"""
**Model:**

:math:`B(t) = \exp\left(-\kappa_1 \vert t \vert^{d_1}\right) \exp\left(-\kappa_2 \vert t\vert^{d_2}\right)`
"""  
def prodstrexp(t,decay1,stretch1,decay2,stretch2):
    strexp1 = np.exp(-decay1*abs(t)**stretch1)
    strexp2 = np.exp(-decay2*abs(t)**stretch2)
    return strexp1*strexp2
# Create model
bg_prodstrexp = Model(prodstrexp,axis='t')
bg_prodstrexp.description = 'Product of two stretched exponentials background model'
# Parameters
bg_prodstrexp.decay1.set(description='Decay rate of 1st component', lb=0, ub=200, par0=0.25, units='μs⁻¹')
bg_prodstrexp.decay2.set(description='Decay rate of 2nd component', lb=0, ub=200, par0=0.25, units='μs⁻¹')
bg_prodstrexp.stretch1.set(description='Stretch factor of 1st component', lb=0, ub=6, par0=1, units='')
bg_prodstrexp.stretch2.set(description='Stretch factor of 2nd component', lb=0, ub=6, par0=1, units='')
# Add documentation
bg_prodstrexp.__doc__ = _docstring(bg_prodstrexp,notes)



#=======================================================================================
#                                     bg_sumstrexp
#=======================================================================================
notes = r"""
**Model:**

:math:`B(t) = A_1\exp \left(-\kappa_1 \vert t \vert^{d_1}\right) + (1-A_1)\exp\left(-\kappa_2 \vert t \vert^{d_2}\right)`
""" 
def sumstrexp(t,decay1,stretch1,weight1,decay2,stretch2):
    strexp1 = np.exp(-decay1*abs(t)**stretch1)
    strexp2 = np.exp(-decay2*abs(t)**stretch2)
    return weight1*strexp1 + (1-weight1)*strexp2
# Create model
bg_sumstrexp = Model(sumstrexp,axis='t')
bg_sumstrexp.description = 'Sum of two stretched exponentials background model'
# Parameters
bg_sumstrexp.decay1.set(description='Decay rate of 1st component', lb=0, ub=200, par0=0.25, units='μs⁻¹')
bg_sumstrexp.decay2.set(description='Decay rate of 2nd component', lb=0, ub=200, par0=0.25, units='μs⁻¹')
bg_sumstrexp.weight1.set(description='Weight of the 1st component', lb=0, ub=1, par0=0.5, units='')
bg_sumstrexp.stretch1.set(description='Stretch factor of 1st component', lb=0, ub=6, par0=1, units='')
bg_sumstrexp.stretch2.set(description='Stretch factor of 2nd component', lb=0, ub=6, par0=1, units='')
# Add documentation
bg_sumstrexp.__doc__ = _docstring(bg_sumstrexp,notes)


#=======================================================================================
#                                     bg_poly1
#=======================================================================================
notes =  r"""
**Model:**

:math:`B(t) = p_0 + p_1 t`
"""  
def poly1(t,p0,p1):
    return np.polyval([p1,p0],abs(t))
# Create model
bg_poly1 = Model(poly1,axis='t')
bg_poly1.description = 'Polynomial 1st-order background model'
# Parameters
bg_poly1.p0.set(description='Intercept', lb=0, ub=200, par0=1, units='')
bg_poly1.p1.set(description='1st order weight', lb=-200, ub=200, par0=-1, units='μs⁻¹')
# Add documentation
bg_poly1.__doc__ = _docstring(bg_poly1,notes)



#=======================================================================================
#                                     bg_poly2
#=======================================================================================
notes =  r"""
**Model:**

:math:`B(t) = p_0 + p_1 t + p_2 t^2`
"""  
def poly2(t,p0,p1,p2):
    return np.polyval([p2,p1,p0],abs(t))
# Create model
bg_poly2 = Model(poly2,axis='t')
bg_poly2.description = 'Polynomial 2nd-order background model'
# Parameters
bg_poly2.p0.set(description='Intercept', lb=0, ub=200, par0=1, units='')
bg_poly2.p1.set(description='1st order weight', lb=-200, ub=200, par0=-1, units=r'μs\ :sup:`-1`')
bg_poly2.p2.set(description='2nd order weight', lb=-200, ub=200, par0=-1, units=r'μs\ :sup:`-2`')
# Add documentation
bg_poly2.__doc__ = _docstring(bg_poly2,notes)


#=======================================================================================
#                                     bg_poly2
#=======================================================================================
notes =  r"""
**Model:**

:math:`B(t) = p_0 + p_1 t + p_2 t^2 + p_3 t^3`
"""  
def poly3(t,p0,p1,p2,p3):
    return np.polyval([p3,p2,p1,p0],abs(t))
# Create model
bg_poly3 = Model(poly3,axis='t')
bg_poly3.description = 'Polynomial 3rd-order background model'
# Parameters
bg_poly3.p0.set(description='Intercept', lb=0, ub=200, par0=1, units='')
bg_poly3.p1.set(description='1st order weight', lb=-200, ub=200, par0=-1, units=r'μs\ :sup:`-1`')
bg_poly3.p2.set(description='2nd order weight', lb=-200, ub=200, par0=-1, units=r'μs\ :sup:`-2`')
bg_poly3.p3.set(description='3rd order weight', lb=-200, ub=200, par0=-1, units=r'μs\ :sup:`-3`')
# Add documentation
bg_poly3.__doc__ = _docstring(bg_poly3,notes)

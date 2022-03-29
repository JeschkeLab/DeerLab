# bg_models.py - Background parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.
 
import numpy as np
import math as m
from numpy import pi
import inspect
from deerlab.dipolarkernel import dipolarkernel
from deerlab.utils import formatted_table
from deerlab.model import Model
from scipy.special import gamma, hyp2f1, sici
from deerlab.constants import *

#---------------------------------------------------------------------------------------
def hyp2f1_repro(a,b,c,z): 
    """
    Gauss Hypergeometric function 2F1 for |z|>1 and non-integer (a-b) based on its "reciprocation" form
    Reference: https://functions.wolfram.com/07.23.17.0057.01
    """
    return gamma(b - a)*gamma(c)/(gamma(b)*gamma(c - a)*(-z)**a)*hyp2f1(a, a - c + 1, a - b + 1, 1/z) + \
         (gamma(a - b)*gamma(c))/(gamma(a)*gamma(c - b)*(-z)**b)*hyp2f1(b, b - c + 1, b - a + 1, 1/z)
#---------------------------------------------------------------------------------------


def _docstring(model,notes):
#---------------------------------------------------------------------------------------
    args = model._parameter_list(order='vector')
    args.insert(model._constantsInfo[0]['argidx'],model._constantsInfo[0]['argkey'])

    parameters = ''
    for arg in args:
        if arg==model._constantsInfo[0]['argkey']:
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

    **Parameter Table**
    """) 
    string += '\n'
    string += '\n'
    table = []
    table.append(['Name','Lower','Upper','Type','Frozen','Unit','Description'])  
    for n, paramname in enumerate(model._parameter_list(order='vector')): 
        param_str = f'``{paramname}``'
        lb_str = f'{np.atleast_1d(getattr(model,paramname).lb)[0]:5.3g}'
        ub_str = f'{np.atleast_1d(getattr(model,paramname).ub)[0]:5.3g}'
        linear_str = "linear" if np.all(getattr(model,paramname).linear) else "nonlin"
        frozen_str = "Yes" if np.all(getattr(model,paramname).frozen) else "No"
        unit_str = str(getattr(model,paramname).unit)
        desc_str = str(getattr(model,paramname).description)
        table.append([param_str,lb_str,ub_str,linear_str,frozen_str,unit_str,desc_str])
    string += formatted_table(table)
    string += f'\n{notes}'

    return string
#---------------------------------------------------------------------------------------

#=======================================================================================
#                                     bg_hom3d
#=======================================================================================
notes = r"""
**Model**

This model describes the inter-molecular interaction of one observer spin with a 3D
homogenous distribution of spins of concentration `c_s`

.. image:: ../images/model_scheme_bg_hom3d.png
   :width: 350px

The expression for this model is

.. math::
   B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\lambda c_s D |t|\right)`

where `c_s` is the spin concentration (entered in spins/m\ :sup:`3` into this expression) and D is the dipolar constant

.. math::
   D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}
"""  
def _hom3d(t,conc,lam):
    # Unit conversion    
    conc = conc*1e-6*1e3*Nav # umol/L -> mol/L -> mol/m^3 -> spins/m^3
    # Compute background function
    κ = 8*pi**2/9/m.sqrt(3)
    B = np.exp(-κ*lam*conc*D*np.abs(t*1e-6))
    return B
# Create model
bg_hom3d = Model(_hom3d,constants='t')
bg_hom3d.description = 'Background from a homogeneous distribution of spins in a 3D medium'
# Add parameters
bg_hom3d.conc.set(description='Spin concentration', lb=0.01, ub=5000, par0=50, unit='μM')
bg_hom3d.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, unit='')
# Add documentation
bg_hom3d.__doc__ = _docstring(bg_hom3d,notes)


#=======================================================================================
#                                     bg_hom3d_phase
#=======================================================================================
notes = r"""
**Model**

This model describes the phase shift due to inter-molecular interactions between one observer spin with a 3D homogenous distribution of spins of concentration `c_s`

The expression for this model is

.. math::
   B(t) = \mathrm{exp}\left(-\ii\frac{8\pi}{9\sqrt{3}}(\sqrt{3} + \mathrm{ln}(2-\sqrt{3}))\lambda c_s D |t|\right)`

where `c_s` is the spin concentration (entered in spins/m\ :sup:`3` into this expression) and D is the dipolar constant

.. math::
   D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

"""  
def _hom3dphase(t,conc,lam):
    # Unit conversion    
    conc = conc*1e-6*1e3*Nav # umol/L -> mol/L -> mol/m^3 -> spins/m^3
    # Compute background function
    ξ = 8*pi/9/np.sqrt(3)*(np.sqrt(3)+np.log(2-np.sqrt(3)))*D
    B = np.exp(-1j*ξ*lam*conc*(t*1e-6))

    return B
# Create model
bg_hom3d_phase = Model(_hom3dphase,constants='t')
bg_hom3d_phase.description = 'Phase shift from a homogeneous distribution of spins in a 3D medium'
# Add parameters
bg_hom3d_phase.conc.set(description='Spin concentration', lb=0.01, ub=5000, par0=50, unit='μM')
bg_hom3d_phase.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, unit='')
# Add documentation
bg_hom3d_phase.__doc__ = _docstring(bg_hom3d_phase,notes)

#=======================================================================================
#                                     bg_hom3dex
#=======================================================================================
notes = r"""
**Model**

.. image:: ../images/model_scheme_bg_hom3dex.png
   :width: 350px

This implements a hard-shell excluded-volume model, with spin concentration `c_s` (in μM) and the radius of the spherical excluded volume `R_\mathrm{ex}` (in nm).

The expression for this model is

.. math:: B(t) = \exp \Bigg(- c_\mathrm{s}\lambda_k \bigg( V_\mathrm{ex} K_0(t, R_\mathrm{ex}) + \mathcal{I}_\mathrm{S}(t) \bigg) 

where `\mathcal{I}_\mathrm{S}(t)` is an integral without analytical form given by 

.. math:: \mathcal{I}_\mathrm{S}(t) = \frac{4\pi}{3} D\,t \int_0^1 \mathrm{d}z~(1 - 3z^2) ~ \mathrm{S_i}\left( \frac{D\,t (1 - 3z^2)}{R_\mathrm{ex}^3 } \right)  

where `\mathrm{S_i}` is the sine integral function and `D` is the dipolar constant

.. math:: D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

"""  
def _hom3dex(t,conc,rex,lam):
    # Conversion: µmol/L -> mol/L -> mol/m^3 -> spins/m^3
    conc = conc*1e-6*1e3*Nav 

    # Excluded volume
    Vex = 4*np.pi/3*(rex*1e-9)**3

    # Averaging integral
    z = np.linspace(0,1,1000)[np.newaxis,:]
    Dt = D*t[:,np.newaxis]*1e-6
    Is =  4*np.pi/3*np.trapz(Dt*(1-3*z**2)*sici((Dt*(1-3*z**2))/((rex*1e-9)**3))[0],z,axis=1)
    
    # Background function
    C_k = -Vex + Is + np.squeeze(Vex*(dipolarkernel(t,rex,integralop=False)))
    B = np.exp(-lam*conc*C_k)

    return B
# Create model
bg_hom3dex = Model(_hom3dex,constants='t')
bg_hom3dex.description = 'Background from a homogeneous distribution of spins with excluded volume'
# Add parameters
bg_hom3dex.conc.set(description='Spin concentration', lb=0.01, ub=5000, par0=50, unit='μM')
bg_hom3dex.rex.set(description='Exclusion radius', lb=0.01, ub=20, par0=1, unit='nm')
bg_hom3dex.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, unit='')
# Add documentation
bg_hom3dex.__doc__ = _docstring(bg_hom3dex,notes)



#=======================================================================================
#                                     bg_hom3dex_phase
#=======================================================================================
notes = r"""
**Model**

.. image:: ../images/model_scheme_bg_hom3dex.png
   :width: 350px

This implements the phase-shift arising from a hard-shell excluded-volume model, with spin concentration `c_s` (in μM) and the radius of the spherical excluded volume `R_\mathrm{ex}` (in nm).

The expression for this model is

.. math:: B(t) = \exp \Bigg(- \ii c_\mathrm{s}\lambda_k \bigg( V_\mathrm{ex} \mathrm{Im}\{\mathcal{K}_0(t, R_\mathrm{ex})\} + \mathcal{I}_\mathrm{C}(t) \bigg) 

where `\mathcal{I}_\mathrm{C}(t)` is an integral without analytical form given by 

.. math:: \mathcal{I}_\mathrm{C}(t) = \frac{4\pi}{3} D\,t \int_0^1 \mathrm{d}z~(1 - 3z^2) ~ \mathrm{C_i}\left( \frac{D\,t (1 - 3z^2)}{R_\mathrm{ex}^3 } \right)  

where `\mathrm{C_i}` is the cosine integral function and `D` is the dipolar constant

.. math:: D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

"""  
def _hom3dex_phase(t,conc,rex,lam):    
    # Conversion: µmol/L -> mol/L -> mol/m^3 -> spins/m^3
    conc = conc*1e-6*1e3*Nav 

    # Excluded volume
    Vex = 4*np.pi/3*(rex*1e-9)**3

    # Averaging integral
    ξ = 8*pi**2/9/np.sqrt(3)*(np.sqrt(3)+np.log(2-np.sqrt(3)))/np.pi*D
    z = np.linspace(0,1,1000)[np.newaxis,:]
    Dt = D*t[:,np.newaxis]*1e-6
    Ic =  -ξ*(t*1e-6) + 4*np.pi/3*np.trapz(Dt*(1-3*z**2)*sici((Dt*np.abs(1-3*z**2))/((rex*1e-9)**3))[1],z,axis=1)
    
    # Background function
    C_k = - Ic - np.squeeze(Vex*(dipolarkernel(t,rex,integralop=False,complex=True)).imag)
    B = np.exp(-1j*lam*conc*C_k)

    return B
# Create model
bg_hom3dex_phase = Model(_hom3dex_phase,constants='t')
bg_hom3dex_phase.description = 'Phase shift from a homogeneous distribution of spins with excluded volume'
# Add parameters
bg_hom3dex_phase.conc.set(description='Spin concentration', lb=0.01, ub=5000, par0=50, unit='μM')
bg_hom3dex_phase.rex.set(description='Exclusion radius', lb=0.01, ub=20, par0=1, unit='nm')
bg_hom3dex_phase.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, unit='')
# Add documentation
bg_hom3dex_phase.__doc__ = _docstring(bg_hom3dex_phase,notes)


#=======================================================================================
#                                     bg_homfractal
#=======================================================================================
notes = r"""
**Model**

This implements the background due to a homogeneous distribution of spins in a `fdim`-dimensional
space, with the `fdim`-dimensional spin concentration ``fconc``.
"""  
def _homfractal(t,fconc,fdim,lam):
    # Fractal dimension (not defined for d=[0, 1.5, 3, 4.5, 6])
    d = float(fdim)

    # Unit conversion of concentration
    conc = fconc*1e-6*(np.power(10,d))*Nav # µmol/dm^d -> mol/m^d -> spins/m^d

    # Compute prefactor
    if d==3:
        κd = 8*np.pi**2/9/np.sqrt(3) # d->3 limit of general expression
    elif d==1.5: 
        κd = 8.71135  # d->1.5 limit of general expression
    elif d==4.5:
        κd = 5.35506  # d->4.5 limit of general expression
    else:
        κd = 2/9*(-1)**(-d/3)*pi*np.cos(d*pi/6)*gamma(-d/3)*(
                (-1 + (-1)**(d/3))*np.sqrt(3*np.pi)*gamma(1+d/3)/gamma(3/2+d/3)
                + 6*hyp2f1_repro(1/2, -d/3, 3/2, 3)
              )
        κd = κd.real  # Imaginary part is always negligible 
    
    # Compute background function
    B = np.exp(-κd*lam*conc*abs(D*t*1e-6)**(d/3))
    return B
# ======================================================================
# Create model
bg_homfractal = Model(_homfractal, constants='t')
bg_homfractal.description = 'Background from homogeneous spin distribution in a space of fractal dimension'
# Add parameters
bg_homfractal.fconc.set(description='Fractal concentration of spins', lb=1e-20, ub=1e20, par0=1.0e-6, unit='μmol/dmᵈ')
bg_homfractal.fdim.set(description='Fractal dimensionality', lb=0.01, ub=5.99, par0=2.2, unit='')
bg_homfractal.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, unit='')
# Add documentation
bg_homfractal.__doc__ = _docstring(bg_homfractal, notes)




#=======================================================================================
#                                     bg_homfractal_phase
#=======================================================================================
notes = r"""
**Model**

This implements the phase shift due to a homogeneous distribution of spins in a `d`-dimensional space, with `d`-dimensional spin concentration ``c_d``.
"""  
def _homfractal_phase(t,fconc,fdim,lam):
    # Fractal dimension (not defined for d=[0, 1.5, 3, 4.5, 6])
    d = float(fdim) 

    # Unit conversion of concentration    
    fconc = fconc*1e-6*(np.power(10,d))*Nav # umol/dm^d -> mol/m^d -> spins/m^d

    # Compute constant
    if d==3: 
        ξd = 0.33462*D**(d/3) # Limit of d->3 of equation below 
    elif d==1.5: 
        ξd = 1j*np.inf # Limit of d->1.5 of equation below 
    elif d==4.5:
        ξd = 1j*np.inf # Limit of d->4.5 of equation below 
    else:
        ξd = D**(d/3)*pi**(3/2)/9/gamma(3/2 + d/3) * (
                np.sqrt(3)*pi*np.cos(d*pi/6)/np.cos(d*pi/3) 
                - 2**(2+2*d/3)*3**(1 + d/3)*gamma(-1-2*d/3)*np.sin(d*pi/6)*gamma(3/2+d/3)*hyp2f1((-3-2*d)/6, -d/3, (3-2*d)/6, 1/3)/gamma((3-2*d)/6) 
            ) 

    # Compute background function
    B = np.exp(-1j*ξd*fconc*lam*np.sign(t)*abs(t*1e-6)**(d/3))
    return B
# ======================================================================
# Create model
bg_homfractal_phase = Model(_homfractal_phase,constants='t')
bg_homfractal_phase.description = 'Phase shift from a homogeneous distribution of spins in a fractal medium'
# Add parameters
bg_homfractal_phase.fconc.set(description='Fractal concentration of spins', lb=1e-20, ub=1e20, par0=1.0e-6, unit='μmol/dmᵈ')
bg_homfractal_phase.fdim.set(description='Fractal dimensionality', lb=0.01, ub=5.99, par0=2.2, unit='')
bg_homfractal_phase.lam.set(description='Pathway amplitude', lb=0, ub=1, par0=1, unit='')
# Add documentation
bg_homfractal_phase.__doc__ = _docstring(bg_homfractal_phase,notes)


#=======================================================================================
#                                     bg_exp
#=======================================================================================
notes= r"""
**Model**

.. math::

   B(t) = \exp\left(-\kappa \vert t \vert\right)

Although the ``bg_exp`` model has the same functional form as ``bg_hom3d``, it is distinct since its 
parameter is a decay rate constant and not a spin concentration like for ``bg_hom3d``. 
    """  
def _exp(t,decay):
    return np.exp(-decay*np.abs(t))
# Create model
bg_exp = Model(_exp,constants='t')
bg_exp.description = 'Exponential background model'
# Add parameters
bg_exp.decay.set(description='Decay rate', lb=0, ub=200, par0=0.35, unit='μs⁻¹')
# Add documentation
bg_exp.__doc__ = _docstring(bg_exp,notes)



#=======================================================================================
#                                     bg_strexp
#=======================================================================================
notes = r"""
**Model**

.. math::

    B(t) = \exp\left(-\kappa \vert t\vert^{d}\right)

Although the ``bg_strexp`` model has the same functional form as ``bg_homfractal``, it is distinct since its
first parameter is a decay rate constant and not a spin concentration like for ``bg_homfractal``.
"""  
def _strexp(t,decay,stretch):
    return np.exp(-decay*abs(t)**stretch)
# Create model
bg_strexp = Model(_strexp,constants='t')
bg_strexp.description = 'Stretched exponential background model'
# Add parameters
bg_strexp.decay.set(description='Decay rate', lb=0, ub=200, par0=0.25, unit='μs⁻¹')
bg_strexp.stretch.set(description='Stretch factor', lb=0, ub=6, par0=1, unit='')
# Add documentation
bg_strexp.__doc__ = _docstring(bg_strexp,notes)



#=======================================================================================
#                                     bg_prodstrexp
#=======================================================================================
notes = r"""
**Model**

:math:`B(t) = \exp\left(-\kappa_1 \vert t \vert^{d_1}\right) \exp\left(-\kappa_2 \vert t\vert^{d_2}\right)`
"""  
def _prodstrexp(t,decay1,stretch1,decay2,stretch2):
    strexp1 = np.exp(-decay1*abs(t)**stretch1)
    strexp2 = np.exp(-decay2*abs(t)**stretch2)
    return strexp1*strexp2
# Create model
bg_prodstrexp = Model(_prodstrexp,constants='t')
bg_prodstrexp.description = 'Product of two stretched exponentials background model'
# Add parameters
bg_prodstrexp.decay1.set(description='Decay rate of 1st component', lb=0, ub=200, par0=0.25, unit='μs⁻¹')
bg_prodstrexp.decay2.set(description='Decay rate of 2nd component', lb=0, ub=200, par0=0.25, unit='μs⁻¹')
bg_prodstrexp.stretch1.set(description='Stretch factor of 1st component', lb=0, ub=6, par0=1, unit='')
bg_prodstrexp.stretch2.set(description='Stretch factor of 2nd component', lb=0, ub=6, par0=1, unit='')
# Add documentation
bg_prodstrexp.__doc__ = _docstring(bg_prodstrexp,notes)



#=======================================================================================
#                                     bg_sumstrexp
#=======================================================================================
notes = r"""
**Model**

:math:`B(t) = A_1\exp \left(-\kappa_1 \vert t \vert^{d_1}\right) + (1-A_1)\exp\left(-\kappa_2 \vert t \vert^{d_2}\right)`
""" 
def _sumstrexp(t,decay1,stretch1,weight1,decay2,stretch2):
    strexp1 = np.exp(-decay1*abs(t)**stretch1)
    strexp2 = np.exp(-decay2*abs(t)**stretch2)
    return weight1*strexp1 + (1-weight1)*strexp2
# Create model
bg_sumstrexp = Model(_sumstrexp,constants='t')
bg_sumstrexp.description = 'Sum of two stretched exponentials background model'
# Add parameters
bg_sumstrexp.decay1.set(description='Decay rate of 1st component', lb=0, ub=200, par0=0.25, unit='μs⁻¹')
bg_sumstrexp.decay2.set(description='Decay rate of 2nd component', lb=0, ub=200, par0=0.25, unit='μs⁻¹')
bg_sumstrexp.weight1.set(description='Weight of the 1st component', lb=0, ub=1, par0=0.5, unit='')
bg_sumstrexp.stretch1.set(description='Stretch factor of 1st component', lb=0, ub=6, par0=1, unit='')
bg_sumstrexp.stretch2.set(description='Stretch factor of 2nd component', lb=0, ub=6, par0=1, unit='')
# Add documentation
bg_sumstrexp.__doc__ = _docstring(bg_sumstrexp,notes)


#=======================================================================================
#                                     bg_poly1
#=======================================================================================
notes =  r"""
**Model**

:math:`B(t) = p_0 + p_1 t`
"""  
def _poly1(t,p0,p1):
    return np.polyval([p1,p0],abs(t))
# Create model
bg_poly1 = Model(_poly1,constants='t')
bg_poly1.description = 'Polynomial 1st-order background model'
# Add parameters
bg_poly1.p0.set(description='Intercept', lb=0, ub=200, par0=1, unit='')
bg_poly1.p1.set(description='1st order weight', lb=-200, ub=200, par0=-1, unit='μs⁻¹')
# Add documentation
bg_poly1.__doc__ = _docstring(bg_poly1,notes)



#=======================================================================================
#                                     bg_poly2
#=======================================================================================
notes =  r"""
**Model**

:math:`B(t) = p_0 + p_1 t + p_2 t^2`
"""
def _poly2(t,p0,p1,p2):
    return np.polyval([p2,p1,p0],abs(t))
# Create model
bg_poly2 = Model(_poly2,constants='t')
bg_poly2.description = 'Polynomial 2nd-order background model'
# Add parameters
bg_poly2.p0.set(description='Intercept', lb=0, ub=200, par0=1, unit='')
bg_poly2.p1.set(description='1st order weight', lb=-200, ub=200, par0=-1, unit=r'μs\ :sup:`-1`')
bg_poly2.p2.set(description='2nd order weight', lb=-200, ub=200, par0=-1, unit=r'μs\ :sup:`-2`')
# Add documentation
bg_poly2.__doc__ = _docstring(bg_poly2,notes)


#=======================================================================================
#                                     bg_poly2
#=======================================================================================
notes =  r"""
**Model**

:math:`B(t) = p_0 + p_1 t + p_2 t^2 + p_3 t^3`
"""
def _poly3(t,p0,p1,p2,p3):
    return np.polyval([p3,p2,p1,p0],abs(t))
# Create model
bg_poly3 = Model(_poly3,constants='t')
bg_poly3.description = 'Polynomial 3rd-order background model'
# Add parameters
bg_poly3.p0.set(description='Intercept', lb=0, ub=200, par0=1, unit='')
bg_poly3.p1.set(description='1st order weight', lb=-200, ub=200, par0=-1, unit=r'μs\ :sup:`-1`')
bg_poly3.p2.set(description='2nd order weight', lb=-200, ub=200, par0=-1, unit=r'μs\ :sup:`-2`')
bg_poly3.p3.set(description='3rd order weight', lb=-200, ub=200, par0=-1, unit=r'μs\ :sup:`-3`')
# Add documentation
bg_poly3.__doc__ = _docstring(bg_poly3,notes)

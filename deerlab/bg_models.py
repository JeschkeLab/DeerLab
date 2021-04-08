# bg_models.py - Background parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.
 
import numpy as np
import math as m
import scipy as scp
from numpy import pi
import inspect
from deerlab.utils import load_exvolume_redfactor, metadata


# =================================================================
def docstr_header1(title,fcnstr): 
    """ Definition of the header for all experiment models taking lambda as well"""
    return f"""
{title}

The function takes a list or array of parameters and returns the calculated background model::

        B = {fcnstr}(t,param)
        B = {fcnstr}(t,param,lam)

The built-in information on the model can be accessed via its attributes::

        {fcnstr}.parameters  # String list of parameter names
        {fcnstr}.units       # String list of metric units of parameters
        {fcnstr}.start       # List of values used as start values during optimization 
        {fcnstr}.lower       # List of values used as lower bounds during optimization
        {fcnstr}.upper       # List of values used as upper bounds during optimization 


Parameters
----------
t : array_like
    Time axis, in microseconds.
param : array_like
    List of model parameter values.
lam : float scalar
    Pathway amplitude. If not specified it is set to 1.

Returns
-------
    
B : ndarray
    Dipolar background. 
"""
# =================================================================

# =================================================================
def docstr_header2(title,fcnstr): 
    """Definition of the header for all experiment models"""
    return f"""
{title}

The function takes a list or array of parameters and returns the calculated background model::

        B = {fcnstr}(t,param)

The built-in information on the model can be accessed via its attributes::

        {fcnstr}.parameters  # String list of parameter names
        {fcnstr}.units       # String list of metric units of parameters
        {fcnstr}.start       # List of values used as start values during optimization 
        {fcnstr}.lower       # List of values used as lower bounds during optimization
        {fcnstr}.upper       # List of values used as upper bounds during optimization 


Parameters
----------
t : array_like
    Time axis, in microseconds.
param : array_like
    List of model parameter values.

Returns
-------
    
B : ndarray
    Dipolar background. 
"""
# =================================================================

# =================================================================
def docstring(takes_lambda=False):
    """
    Decorator: Insert docstring header to a pre-existing docstring
    """
    sep="\n"
    def _decorator(func):
        docstr = func.__doc__
        title = docstr.split("Model parameters:",1)[0]
        docstr = docstr.replace(title,"")
        if takes_lambda:
            func.__doc__ = sep.join([docstr_header1(title,func.__name__),docstr])
        else:
            func.__doc__ = sep.join([docstr_header2(title,func.__name__),docstr])
        return func
    return _decorator
# =================================================================

# =================================================================
def _parsargs(t,p,npar):
    t,p = np.atleast_1d(t,p)
    # Check that the correct number of parmameters have been specified
    if len(p)!=npar:
        raise ValueError(f'The model function requires {npar} parameters, but {len(p)} are provided.')
 
    return t,p
# =================================================================


# =================================================================
@metadata(
parameters = ['Concentration of pumped spins'],
units = ['μM'],
start = np.asarray([50]),
lower = np.asarray([0.01]),
upper = np.asarray([5000]))
@docstring(takes_lambda=True)
def bg_hom3d(t,param,lam=1):
    r"""
Background from homogeneous distribution of spins in a 3D medium
    
Notes
-----

**Model:**

This model describes the inter-molecular interaction of one observer spin with a 3D homogenous distribution of pump-spins of concentration `c_p`

.. image:: ../images/model_scheme_bg_hom3d.png
   :width: 350px

The expression for this model is

.. math::
   B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\lambda c_p D |t|\right)`

where `c_p` is the pumped-spin concentration (entered in spins/m\ :sup:`3` into this expression) and D is the dipolar constant

.. math::
   D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

============== =============== ============= ============= ============= =================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= =================================
``param[0]``   :math:`c_p`        50            0.01           5000       Pumped spin concentration (μM)
============== =============== ============= ============= ============= =================================
    """  
    t,param = _parsargs(t,param,npar=1) 

    conc = param            # concentration, µM
    Nav = 6.02214076e23      # Avogadro constant, mol^-1
    muB = 9.2740100783e-24  # Bohr magneton, J/T (CODATA 2018 value)
    mu0 = 1.25663706212e-6  # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
    h = 6.62607015e-34      # Planck constant, J/Hz (CODATA 2018)
    ge = 2.00231930436256   # free-electron g factor (CODATA 2018 value)
    hbar = h/2/pi         # reduced Planck constant, J/(rad/s)

    D = (mu0/4/pi)*(muB*ge)**2/hbar   # dipolar constant, m^3 s^-1
    
    conc = conc*1e-6*1e3*Nav # umol/L -> mol/L -> mol/m^3 -> spins/m^3
    
    # Compute background function
    B = np.exp(-8*pi**2/9/m.sqrt(3)*lam*conc*D*np.abs(t*1e-6))
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Spin concentration','Exclusion distance'],
units = ['μM','nm'],
start = np.asarray([50,   1]),
lower = np.asarray([0.01, 0.01]),
upper = np.asarray([5000, 20]))
@docstring(takes_lambda=True)
def bg_hom3dex(t,param,lam=1):
    r"""
Background from homogeneous distribution of spins with excluded-volume effects
    
Notes
-----

**Model:**

.. image:: ../images/model_scheme_bg_hom3dex.png
   :width: 350px

This implements a hard-shell excluded-volume model, with pumped spin concentration ``c`` (first parameter, in μM) and distance of closest approach ``R`` (second parameter, in nm).

The expression for this model is

.. math:: B(t) = \mathrm{exp}\left(-\frac{8\pi^2}{9\sqrt{3}}\alpha(R) \lambda c D |t|\right)`

where :math:`c` is the spin concentration (entered in spins/m\ :sup:`3` into this expression) and :math:`D` is the dipolar constant

.. math:: D = \frac{\mu_0}{4\pi}\frac{(g_\mathrm{e}\mu_\mathrm{B})^2}{\hbar}

The function :math:`\alpha(R)` of the exclusion distance :math:`R` captures the excluded-volume effect. It is a smooth function, but doesn't have an analytical representation. For details, see `Kattnig et al, J.Phys.Chem.B 2013, 117, 16542 <https://pubs.acs.org/doi/abs/10.1021/jp408338q>`_.

============== =============== ============= ============= ============= =================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= =================================
``param[0]``    :math:`c`              50         0.01          5000          Spin concentration (μM)
``param[1]``    :math:`R`              1          0.1            20           Exclusion distance (nm)
============== =============== ============= ============= ============= =================================
    """  
    t,param = _parsargs(t,param, npar=2) 
    
    # Load precalculated reduction factor look-up table (Kattnig Eq.(18))
    dR_tab,alphas_tab = load_exvolume_redfactor()

    # Get parameters
    conc = param[0] # µM
    R = param[1]    # nm

    NAv = 6.02214076e23      # Avogadro constant, mol^-1
    conc = conc*1e-6*1e3*NAv # umol/L -> mol/L -> mol/m^3 -> spins/m^3
    ge = 2.00231930436256   # free-electron g factor (CODATA 2018 value)
    mu0 = 1.25663706212e-6  # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
    muB = 9.2740100783e-24  # Bohr magneton, J/T (CODATA 2018 value)
    h = 6.62607015e-34      # Planck constant, J/Hz (CODATA 2018)
    hbar = h/2/pi

    A = (mu0/4/pi)*(ge*muB)**2/hbar # Eq.(6) m^3 s^-1

    # Calculate reduction factor (Eq.(18))
    if R==0:
        alpha = 1
    else:
        dR = A*abs(t*1e-6)/(R*1e-9)**3 # unitless
        
        # Use interpolation of look-up table for small dR
        small = dR < max(dR_tab)
        alpha = np.zeros(np.shape(dR))
        alpha[small] = np.interp(dR[small], dR_tab, alphas_tab)
        
        # For large dR, use limiting dR->inf expression
        alpha[~small] = 1 - (3/2/pi)*np.sqrt(3)/dR[~small]

    K = 8*pi**2/9/np.sqrt(3)*A*abs(t*1e-6)*alpha # Eq.(17)
    B = np.exp(-lam*conc*K) # Eq.(13)
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Fractal Concentration of pumped spins','Fractal dimensionality'],
units = ['μmol/dmᵈ',''],
start = np.asarray([50,   3]),
lower = np.asarray([0.01, 0+np.finfo(float).eps]),
upper = np.asarray([5000, 6-np.finfo(float).eps]))
@docstring(takes_lambda=True)
def bg_homfractal(t,param,lam=1):
    r"""
Background from homogeneous distribution of spins in a fractal medium
    
Notes
-----

**Model:**

This implements the background due to a homogeneous distribution of spins in a d-dimensional space, with d-dimensional spin concentration ``c_d``.

============= ============= ============= ============= ============= ==========================================================
 Variable       Symbol       Start Value   Lower bound   Upper bound      Description
============= ============= ============= ============= ============= ==========================================================
``param[0]``   :math:`c_d`     50          0.01          5000          Pumped spin fractal concentration (μmol/dm\ :sup:`d`)
``param[1]``   :math:`d`       3           0                6          Fractal dimension
============= ============= ============= ============= ============= ==========================================================
    """  
    t,param = _parsargs(t,param,npar=2) 

    # Unpack model paramters
    conc = param[0]         # concentration, umol/dm^d
    d = float(param[1])     # fractal dimension    

    # Natural constants
    NA = 6.02214076e23      # Avogadro constant, mol^-1
    muB = 9.2740100783e-24  # Bohr magneton, J/T (CODATA 2018 value)
    mu0 = 1.25663706212e-6  # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
    h = 6.62607015e-34      # Planck constant, J/Hz (CODATA 2018)
    ge = 2.00231930436256   # free-electron g factor (CODATA 2018 value)
    hbar = h/2/pi         # reduced Planck constant, J/(rad/s)
    D = (mu0/4/pi)*(muB*ge)**2/hbar   # dipolar constant, m^3 s^-1
    # Units conversion of concentration    
    conc = conc*1e-6*(np.power(10,d))*NA # umol/dm^d -> mol/m^d -> spins/m^d
    

    # Compute constants
    if d==3:
        c = -pi/2
        Lam = 4/3/np.sqrt(3)
    else:
        c = np.cos(d*pi/6)*scp.special.gamma(-d/3)
        integrand = lambda z: abs(1-3*z**2)**(d/3)
        Lam,_ = scp.integrate.quad(integrand,0,1,limit=1000)

    # Compute background function
    B = np.exp(4*pi/3*c*Lam*lam*conc*D**(d/3)*abs(t*1e-6)**(d/3))

    return B
 # ======================================================================


# =================================================================
@metadata(
parameters = ['Decay Rate'],
units = ['μs⁻¹'],
start = np.asarray([0.35]),
lower = np.asarray([0]),
upper = np.asarray([200]))
@docstring()
def bg_exp(t,param):
    r"""
Exponential background model
    
Notes
-----

**Model:**

.. math::

   B(t) = \exp\left(-\kappa \vert t \vert\right)

============== =============== ============= ============= ============= ================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= ================================
``param[0]``   :math:`\kappa`      0.35         0            200          Decay rate (μs\ :sup:`-1`)
============== =============== ============= ============= ============= ================================


Although the ``bg_exp`` model has the same functional form as ``bg_hom3d``, it is distinct since its 
parameter is a decay rate constant and not a spin concentration like for ``bg_hom3d``. 
    """  
    t,param = _parsargs(t,param,npar=1) 
    t = np.atleast_1d(t)
    param = np.atleast_1d(param)
    kappa = param[0]
    B = np.exp(-kappa*np.abs(t))
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Decay Rate','Stretch factor'],
units = ['μs⁻¹',''],
start = np.asarray([0.25, 1]),
lower = np.asarray([0,    0]),
upper = np.asarray([200,  6]))
@docstring()
def bg_strexp(t,param):
    r"""
Stretched exponential background model
    
Notes
-----

**Model:**

.. math::

    B(t) = \exp\left(-\kappa \vert t\vert^{d}\right)

============== ================= ============= ============= ============= =================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`\kappa`      0.25              0              200        Decay rate (μs\ :sup:`-d`)
``param[1]``   :math:`d`           1                 0              6          Stretch factor
============== ================= ============= ============= ============= =================================

Although the ``bg_strexp`` model has the same functional form as ``bg_homfractal``, it is distinct since its
first parameter is a decay rate constant and not a spin concentration like for ``bg_homfractal``.
    """  
    t,param = _parsargs(t,param,npar=2) 
    kappa = param[0]         # decay rate, µs^-1
    d = param[1]            # fractal dimension    
    B = np.exp(-kappa*abs(t)**d)
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Decay Rate of 1st component','Stretch factor of 1st component',
              'Decay Rate of 2nd component','Stretch factor of 2nd component'],
units = ['µs^-1','','µs^-1',''],
start = np.asarray([0.25, 1, 0.25, 1]),
lower = np.asarray([ 0,   0,  0,   0]),
upper = np.asarray([200,  6, 200,  6]))
@docstring()
def bg_prodstrexp(t,param):
    r"""
Product of two stretched exponentials background model
    
Notes
-----

**Model:**

:math:`B(t) = \exp\left(-\kappa_1 \vert t \vert^{d_1}\right) \exp\left(-\kappa_2 \vert t\vert^{d_2}\right)`

============== ================= ============= ============= ============= =================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`\kappa_1`    0.25          0              200         1st strexp decay rate
``param[1]``   :math:`d_1`         1             0              6           1st strexp stretch factor
``param[2]``   :math:`\kappa_2`    0.25          0              200         2nd strexp decay rate
``param[3]``   :math:`d_2`         1             0              6           2nd strexp stretch factor
============== ================= ============= ============= ============= =================================
    """  
    t,param = _parsargs(t,param,npar=4) 
    kappa1 = param[0]
    d1 = param[1]
    kappa2 = param[2]
    d2 = param[3]
    strexp1 = np.exp(-kappa1*abs(t)**d1)
    strexp2 = np.exp(-kappa2*abs(t)**d2)
    B = strexp1*strexp2
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Decay Rate of 1st component','Stretch factor of 1st component',
              'Amplitude of 1st component','Decay Rate of 2nd component','Stretch factor of 2nd component'],
units = ['μs⁻¹','','','μs⁻¹',''],
start = np.asarray([0.25, 1, 0.5, 0.25, 1]),
lower = np.asarray([ 0,   0,  0,   0,   0]),
upper = np.asarray([200,  6,  1,  200,  6]))
@docstring()
def bg_sumstrexp(t,param):
    r"""
Sum of two stretched exponentials background model
    
Notes
-----

**Model:**

:math:`B(t) = A_1\exp \left(-\kappa_1 \vert t \vert^{d_1}\right) + (1-A_1)\exp\left(-\kappa_2 \vert t \vert^{d_2}\right)`

============== ================= ============= ============= ============= ========================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= ========================================
``param[0]``   :math:`\kappa_1`     0.25            0            200         1st strexp decay rate (μs\ :sup:`-d`)
``param[1]``   :math:`d_1`          1               0            6           1st strexp stretch factor
``param[2]``   :math:`\kappa_2`     0.25            0            200         2nd strexp decay rate (μs\ :sup:`-d`)
``param[3]``   :math:`d_2`          1               0            6           2nd strexp stretch factor
``param[4]``   :math:`A_1`          0.50            0            1           Relative amplitude
============== ================= ============= ============= ============= ========================================
    """ 
    t,param = _parsargs(t,param,npar=5) 
    kappa1 = param[0]
    d1 = param[1]
    w1 = param[2]
    kappa2 = param[3]
    d2 = param[4]
    strexp1 = np.exp(-kappa1*abs(t)**d1)
    strexp2 = np.exp(-kappa2*abs(t)**d2)
    B = w1*strexp1 + (1-w1)*strexp2
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Intercept','1st-order coefficient'],
units = ['','μs⁻¹'],
start = np.asarray([ 1,   -1 ]),
lower = np.asarray([ 0,  -200]),
upper = np.asarray([200,  200]))
@docstring()
def bg_poly1(t,param):
    r"""
Polynomial 1st-order background model
    
Notes
-----

**Model:**

:math:`B(t) = p_0 + p_1 t`

============== =============== ============= ============= ============= ====================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= ====================================
``param[0]``    :math:`p_0`      1              0            200          Intercept
``param[1]``    :math:`p_1`     -1              -200         200          1st order weight (μs\ :sup:`-1`)
============== =============== ============= ============= ============= ====================================
    """  
    t,param = _parsargs(t,param, npar=2) 
    p = np.copy(np.flip(param))
    p[:-1] = p[:-1]
    B = np.polyval(p,abs(t))
    return B 
# ======================================================================


# =================================================================
@metadata(
parameters = ['Intercept','1st-order coefficient','2nd-order coefficient'],
units = ['','μs⁻¹','μs⁻²'],
start = np.asarray([ 1,   -1 , -1]),
lower = np.asarray([ 0,  -200, -200]),
upper = np.asarray([200,  200,  200]))
@docstring()
def bg_poly2(t,param):
    r"""
Polynomial 2nd-order background model
     
Notes
-----

**Model:**

:math:`B(t) = p_0 + p_1 t + p_2t^2`

============== =============== ============= ============= ============= ===================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= ===================================
``param[0]``   :math:`p_0`       1               0            200          Intercept
``param[1]``   :math:`p_1`      -1               -200         200          1st order weight (μs\ :sup:`-1`)
``param[2]``   :math:`p_2`      -1               -200         200          2nd order weight (μs\ :sup:`-2`)
============== =============== ============= ============= ============= ===================================
    """  
    t,param = _parsargs(t,param,npar=3) 
    p = np.copy(np.flip(param))
    p[:-1] = p[:-1]
    B = np.polyval(p,abs(t))
    return B
# ======================================================================


# =================================================================
@metadata(
parameters = ['Intercept','1st-order coefficient','2nd-order coefficient','3rd-order coefficient'],
units = ['','μs⁻¹','μs⁻²','μs⁻³'],
start = np.asarray([ 1,  -1  , -1,   -1  ]),
lower = np.asarray([ 0,  -200, -200, -200]),
upper = np.asarray([200,  200,  200,  200]))
@docstring()
def bg_poly3(t,param):
    r"""
Polynomial 3rd-order background model
     
Notes
-----

**Model:**

:math:`B(t) = p_0 + p_1t + p_2t^2 + p_3t^3`

============== =============== ============= ============= ============= ===================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= ===================================
``param[0]``   :math:`p_0`        1             0             200          Intercept
``param[1]``   :math:`p_1`       -1          -200             200          1st order weight (μs\ :sup:`-1`)
``param[2]``   :math:`p_2`       -1          -200             200          2nd order weight (μs\ :sup:`-2`)
``param[3]``   :math:`p_3`       -1          -200             200          3rd order weight (μs\ :sup:`-3`)
============== =============== ============= ============= ============= ===================================
    """  
    t,param = _parsargs(t,param,npar=4) 
    p = np.copy(np.flip(param))
    p[:-1] = p[:-1]
    B = np.polyval(p,abs(t))
    return B
# ======================================================================
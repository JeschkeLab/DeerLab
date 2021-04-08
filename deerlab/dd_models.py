# dd_models.py - Distance distribtion parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import math as m
import numpy as np
import scipy.special as spc
import inspect
from deerlab.utils import metadata 

# =================================================================
def docstr_header(title,fcnstr):
    "Definition of the header for all distribution models"
    return f"""
{title}

The function takes a list or array of parameters and returns the calculated distance distribution::

        P = {fcnstr}(r,param)

The built-in information on the model can be accessed via its attributes::

        {fcnstr}.parameters  # String list of parameter names
        {fcnstr}.units       # String list of metric units of parameters
        {fcnstr}.start       # List of values used as start values during optimization 
        {fcnstr}.lower       # List of values used as lower bounds during optimization
        {fcnstr}.upper       # List of values used as upper bounds during optimization 

Parameters
----------
r : array_like
    Distance axis, in nanometers.
    
param : array_like
    List of model parameter values.

Returns
-------

P : ndarray
    Distance distribution.
"""
# =================================================================


# =================================================================
def docstr_example(fcnstr): 
    return f""" 
Examples
--------

Example of the model evaluated at the start values of the parameters:

.. plot::

    import deerlab as dl
    import matplotlib.pyplot as plt 
    import numpy as np 
    model = dl.{fcnstr}
    r = np.linspace(2,5,400)
    info = model() 
    par0 = info['Start']
    P = model(r,par0)
    plt.figure(figsize=[6,3])
    plt.plot(r,P)
    plt.xlabel('r (nm)',fontsize=13)
    plt.ylabel('P (nm⁻¹)',fontsize=13)
    plt.grid(alpha=0.4)
    plt.tick_params(labelsize=12)
    plt.tick_params(labelsize=12)
    plt.tight_layout()
"""
# =================================================================


# =================================================================
def docstring():
    """
    Decorator: Insert docstring header to a pre-existing docstring
    """
    sep="\n"
    def _decorator(func):
        docstr = func.__doc__
        title = docstr.split("Notes",1)[0]
        docstr = docstr.replace(title,"")
        func.__doc__ = sep.join([docstr_header(title,func.__name__),docstr])
        func.__doc__ = sep.join([func.__doc__,docstr_example(func.__name__)])
        return func
    return _decorator
# =================================================================

# =================================================================
def _parsparam(r,p,npar):
    r,p = np.atleast_1d(r,p)
    if len(p)!=npar:
        raise ValueError(f'The model function requires {npar} parameters, but {len(p)} are provided.')
    return r,p
# =================================================================

# =================================================================
def _normalize(r,P):
    if not all(P==0):
        P = P/np.trapz(P,r)
    return P
# =================================================================

# =================================================================
def _multigaussfun(r,r0,sig,a):
    "Compute a distribution with multiple Gaussians"    
    n = len(r0)
    P = np.zeros_like(r)
    for k in range(n):
        P += a[k]*m.sqrt(1/(2*m.pi))*1/sig[k]*np.exp(-0.5*((r-r0[k])/sig[k])**2)
    P = _normalize(r,P)
    return P
# =================================================================

def _multirice3dfun(r,nu,sig,a):
# =================================================================
    "Compute a distribution with multiple Gaussians"    
    N = len(nu)
    nu = np.maximum(nu,0) # to avoid invalid values
    n = 3 # degrees of freedom
    P = np.zeros_like(r)
    for k in range(N):
        s2 = sig[k]**2
        I_scaled = spc.ive(n/2-1, nu[k]*r/s2)
        P =+ a[k]*nu[k]**(n/2-1)/s2*r**(n/2)*np.exp(-(r**2+nu[k]**2)/(2*s2)+nu[k]*r/s2)*I_scaled
    P[P<0] = 0
    P = _normalize(r,P)
    return P
# =================================================================

# =================================================================
@metadata(
parameters = ('Mean','Standard deviation'),
units = ('nm','nm'),
start = np.asarray([3.5, 0.2]),
lower = np.asarray([1, 0.05]),
upper = np.asarray([20, 2.5]))
@docstring()
def dd_gauss(r,param):    
    r"""
Gaussian distribution
    
Notes
-----

**Model:**

:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sigma^2}\right)`

==============  ========================  =============  =============  =============  ===========================
Variable        Symbol                    Start value    Lower bound    Upper bound    Description
==============  ========================  =============  =============  =============  ===========================
``param[0]``    :math:`\left<r\right>`    3.5            1.0            20             Mean (nm)
``param[1]``    :math:`\sigma`            0.2            0.05           2.5            Standard deviation (nm)
==============  ========================  =============  =============  =============  ===========================
    """  
    r,param = _parsparam(r,param,npar=2)
    r0    = [param[0]]
    sigma = [param[1]]
    a = [1.0]
    P = _multigaussfun(r,r0,sigma,a)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Mean of 1st Gaussian', 'Standard deviation of 1st Gaussian', 'Amplitude of 1st Gaussian',
                'Mean of 2nd Gaussian', 'Standard deviation of 2nd Gaussian', 'Amplitude of 2nd Gaussian'),
units = ('nm','nm','','nm','nm',''),
start = np.asarray([2.5, 0.2, 0.5, 3.5, 0.2, 0.5]),
lower = np.asarray([1, 0.05, 0, 1, 0.05, 0]),
upper = np.asarray([20, 2.5, 1, 20, 2.5, 1]))
@docstring()
def dd_gauss2(r,param):
    r"""
Sum of two Gaussian distributions
        
Notes
-----

**Model:**

:math:`P(r) = a_1\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_1}\exp\left(-\frac{(r-\left<r_1\right>)^2}{\sigma_1^2}\right) + a_2\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_2}\exp\left(-\frac{(r-\left<r_2\right>)^2}{\sigma_2^2}\right)`


============== ========================= ============= ============= ============= ======================================
Variable         Symbol                  Start Value   Lower bound   Upper bound      Description
============== ========================= ============= ============= ============= ======================================
``param[0]``   :math:`\left<r_1\right>`     2.5            1.0            20        1st Gaussian mean distance (nm)
``param[1]``   :math:`\sigma_1`             0.2            0.05          2.5        1st Gaussian standard deviation (nm)
``param[2]``   :math:`a_1`                  0.5            0              1         1st Gaussian amplitude
``param[3]``   :math:`\left<r_2\right>`     3.5            1.0            20        2nd Gaussian mean distance (nm)
``param[4]``   :math:`\sigma_2`             0.2            0.05          2.5        2nd Gaussian standard deviation (nm)
``param[5]``   :math:`a_2`                  0.5            0              1         2nd Gaussian amplitude
============== ========================= ============= ============= ============= ======================================
    """
    r,param = _parsparam(r,param,npar=6)
    r0    = [param[0], param[3]]
    sigma = [param[1], param[4]]
    a     = [param[2], param[5]]
    P = _multigaussfun(r,r0,sigma,a)
    return P
# =================================================================
    

# =================================================================
@metadata(
parameters = ('Mean of 1st Gaussian', 'Standard deviation of 1st Gaussian', 'Amplitude of 1st Gaussian',
                'Mean of 2nd Gaussian', 'Standard deviation of 2nd Gaussian', 'Amplitude of 2nd Gaussian',
                'Mean of 3rd Gaussian', 'Standard deviation of 3rd Gaussian', 'Amplitude of 3rd Gaussian'),
units = ('nm','nm','','nm','nm','','nm','nm',''),
start = np.asarray([2.5, 0.2, 0.3, 3.5, 0.2, 0.3, 5, 0.2, 0.3]),
lower = np.asarray([1, 0.05, 0, 1, 0.05, 0, 1, 0.05, 0]),
upper = np.asarray([20, 2.5, 1, 20, 2.5, 1,  20, 2.5, 1]))
@docstring()
def dd_gauss3(r,param):
    r"""
Sum of three Gaussian distributions
        
Notes
-----

**Model:**

:math:`P(r) = a_1\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_1}\exp\left(-\frac{(r-\left<r_1\right>)^2}{\sigma_1^2}\right) + a_2\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_2}\exp\left(-\frac{(r-\left<r_2\right>)^2}{\sigma_2^2}\right) + a_3\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_3}\exp\left(-\frac{(r-\left<r_3\right>)^2}{\sigma_3^2}\right)`

============== ========================== ============= ============= ============= =======================================
Variable         Symbol                   Start Value   Lower bound   Upper bound      Description
============== ========================== ============= ============= ============= =======================================
``param[0]``     :math:`\left<r_1\right>`     2.5           1.0           20         1st Gaussian mean distance (nm)
``param[1]``     :math:`\sigma_1`             0.2           0.05          2.5        1st Gaussian standard deviation (nm)
``param[2]``     :math:`a_1`                  0.3           0             1          1st Gaussian amplitude
``param[3]``     :math:`\left<r_2\right>`     3.5           1.0           20         2nd Gaussian mean distance (nm)
``param[4]``     :math:`\sigma_2`             0.2           0.05          2.5        2nd Gaussian standard deviation (nm)
``param[5]``     :math:`a_2`                  0.3           0             1          2nd Gaussian amplitude
``param[6]``     :math:`\left<r_3\right>`     5.0           1.0           20         3rd Gaussian mean distance (nm)
``param[7]``     :math:`\sigma_3`             0.2           0.05          2.5        3rd Gaussian standard deviation (nm)
``param[8]``     :math:`a_3`                  0.3           0             1          3rd Gaussian amplitude
============== ========================== ============= ============= ============= =======================================
    """
    r,param = _parsparam(r,param,npar=9)
    r0    = [param[0], param[3], param[6]]
    sigma = [param[1], param[4], param[7]]
    a     = [param[2], param[5], param[8]]
    P = _multigaussfun(r,r0,sigma,a)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Mean','Spread','Kurtosis'),
units = ('nm','nm',''),
start = np.asarray([3.5, 0.5,0.5]),
lower = np.asarray([1, 0.05, 0.25]),
upper = np.asarray([20, 5, 15]))
@docstring()
def dd_gengauss(r,param):    
    r"""
Generalized Gaussian distribution model

Notes
-----

**Model:**

.. image:: ../images/model_scheme_dd_gengauss.png
    :width: 450px

:math:`P(r) = \frac{\beta}{2\sigma\Gamma(1/\beta)}\exp\left(-\left(\frac{(r-\left<r\right>)}{\sigma}\right)^\beta \right)`

============== ========================== ============= ============= ============= =======================================
Variable         Symbol                   Start Value   Lower bound   Upper bound      Description
============== ========================== ============= ============= ============= =======================================
``param[0]``   :math:`r_0`                     3.5          1.0              20         Mean (nm)
``param[1]``   :math:`\sigma`                  0.2          0.05             2.5        Spread (nm)
``param[2]``   :math:`\beta`                   5.0          0.25             15         kurtosis
============== ========================== ============= ============= ============= =======================================

    """  
    r,param = _parsparam(r,param,npar=3)
    r0, sigma, beta = param
    x = abs(r-r0)/sigma
    P = beta/(2*sigma*spc.gamma(1/beta))*np.exp(-x**beta)
    P = _normalize(r,P)
    return P
# =================================================================
    

# =================================================================
@metadata(
parameters = ('Center','Spread','Kurtosis'),
units = ('nm','nm',''),
start = np.asarray([3.5, 0.2, 5]),
lower = np.asarray([1, 0.05, -25]),
upper = np.asarray([20, 5, 25]))
@docstring()
def dd_skewgauss(r,param):    
    r"""
Skew Gaussian distribution model


Notes
-----

**Model:**

.. image:: ../images/model_scheme_dd_skewgauss.png
    :width: 650px

:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sqrt(2)\sigma^2}\right)\frac{1}{2}\left(1 + erf\left(\frac{(r-\left<r\right>)}{\sqrt{2}\sigma}\right) \right)`

============== ============== ============= ============= ============= =========================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``   :math:`r_0`          3.5         1.0            20         Location (nm)
``param[1]``   :math:`\sigma`       0.2         0.05           2.5        Spread (nm)
``param[2]``   :math:`\alpha`       5.0         -25            25         Skewness
============== ============== ============= ============= ============= =========================

    """  
    r,param = _parsparam(r,param,npar=3)
    r0, sigma, alpha = param
    x = (r-r0)/sigma/np.sqrt(2)
    P = 1/np.sqrt(2*np.pi)*np.exp(-x**2)*(1 + spc.erf(alpha*x))
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Location','Spread'),
units = ('nm','nm'),
start = np.asarray([3.5, 0.7]),
lower = np.asarray([1, 0.1]),
upper = np.asarray([10, 5])) 
@docstring()
def dd_rice(r,param):    
    r"""
3D-Rice distribution

Notes
-----

**Model:** 

:math:`P(r) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`. This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ============= ============= ============= =======================================
Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =======================================
``param[0]``     :math:`\nu`                3.5           1.0              10          Location (nm)
``param[1]``     :math:`\sigma`             0.7           0.1              5           Spread (nm)
============== ======================== ============= ============= ============= =======================================
    """  
    r,param = _parsparam(r,param,npar=2)
    nu = [param[0]]
    sig = [param[1]]
    a = [1.0]
    P = _multirice3dfun(r,nu,sig,a)
    return P
# =================================================================
    

# =================================================================
@metadata(
parameters = ('Location of 1st Rician', 'Spread of 1st Rician', 'Amplitude of 1st Rician',
                'Location of 2nd Rician', 'Spread of 2nd Rician', 'Amplitude of 2nd Rician'),
units = ('nm','nm','','nm','nm',''),
start = np.asarray([2.5, 0.7, 0.5, 4.0, 0.7, 0.5]),
lower = np.asarray([1, 0.1, 0, 1, 0.1, 0]),
upper = np.asarray([10, 5, 1, 10, 5, 1])) 
@docstring()
def dd_rice2(r,param):
    r"""
Sum of two 3D-Rice distributions

Notes
-----

**Model:**

:math:`P(r) = a_1 R(r,\nu_1,\sigma_1) + a_2 R(r,\nu_2,\sigma_2)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ============= ============= ============= =======================================
Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =======================================
``param[0]``   :math:`\nu_1`                2.5              1         10           1st Rician location (nm)
``param[1]``   :math:`\sigma_1`             0.7             0.1         5           1st Rician spread (nm)
``param[2]``   :math:`a_1`                  0.5              0          1           1st Rician amplitude
``param[3]``   :math:`\nu_2`                 4               1         10           2nd Rician location (nm)
``param[4]``   :math:`\sigma_2`             0.7             0.1         5           2nd Rician spread (nm)
``param[5]``   :math:`a_2`                  0.5              0          1           2nd Rician amplitude
============== ======================== ============= ============= ============= =======================================
    """
    r,param = _parsparam(r,param,npar=6)
    nu  = [param[0], param[3]]
    sig = [param[1], param[4]]
    a   = [param[2], param[5]]
    P = _multirice3dfun(r,nu,sig,a)
    return P
# =================================================================
    

# =================================================================
@metadata(
parameters = ('Location of 1st Rician', 'Spread of 1st Rician', 'Amplitude of 1st Rician',
              'Location of 2nd Rician', 'Spread of 2nd Rician', 'Amplitude of 2nd Rician',
              'Location of 3rd Rician', 'Spread of 3rd Rician', 'Amplitude of 3rd Rician'),
units = ('nm','nm','','nm','nm','','nm','nm',''),
start = np.asarray([2.5, 0.7, 0.3, 3.5, 0.7, 0.3, 5, 0.7, 0.3]),
lower = np.asarray([1, 0.1, 0, 1, 0.1, 0, 1, 0.1, 0]),
upper = np.asarray([10, 5, 1, 10, 5, 1,  10, 5, 1])) 
@docstring()
def dd_rice3(r,param):
    r"""
Sum of three 3D-Rice distributions

Notes
-----

**Model:**
    
:math:`P(r) = a_1 R(r,\nu_1,\sigma_1) + a_2 R(r,\nu_2,\sigma_2) + a_3 R(r,\nu_3,\sigma_3)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ============= ============= ============= =======================================
Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =======================================
``param[0]``   :math:`\nu_1`                2.5           1.0           10          1st Rician location (nm)
``param[1]``   :math:`\sigma_1`             0.7           0.1           5           1st Rician spread (nm)
``param[2]``   :math:`a_1`                  0.3             0           1           1st Rician amplitude
``param[3]``   :math:`\nu_2`                4.0           1.0           10          2nd Rician location (nm)
``param[4]``   :math:`\sigma_2`             0.7           0.1           5           2nd Rician spread (nm)
``param[5]``   :math:`a_2`                  0.3           0             1           2nd Rician amplitude
``param[6]``   :math:`\nu_3`                5.0           1.0           10          3rd Rician location (nm)
``param[7]``   :math:`\sigma_3`             0.7           0.1           5           3rd Rician spread (nm)
``param[8]``   :math:`a_3`                  0.3           0             1           3rd Rician amplitude
============== ======================== ============= ============= ============= =======================================
    """
    r,param = _parsparam(r,param,npar=9)
    nu  = [param[0], param[3], param[6]]
    sig = [param[1], param[4], param[7]]
    a   = [param[2], param[5], param[8]]
    P = _multirice3dfun(r,nu,sig,a)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Number of residues','Segment length','Scaling exponent'),
units = ('','nm',''),
start = np.asarray([50,   0.2, 0.602]),
lower = np.asarray([2,    0.1, 0.33 ]),
upper = np.asarray([1000, 0.4, 1    ])) 
@docstring()
def dd_randcoil(r,param):    
    r"""
Random-coil model for an unfolded peptide/protein

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_randcoil.png
    :width: 25%

:math:`P(r) = \frac{3}{(2\pi\nu_0)^{3/2}}4\pi r^2\exp(-\frac{3 r^2}{\nu_0})`

where :math:`\nu_0 = 3/(12\pi r_0 N \nu)^{3/2}`

============== ============= ============= ============= ============= =======================================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============= ============= ============= ============= =======================================
``param[0]``   :math:`N`         50                2         1000           Number of residues
``param[1]``   :math:`R_0`       0.20           0.10         0.40           Segment length (nm)
``param[2]``   :math:`\nu`       0.602          0.33         1.00           Scaling exponent
============== ============= ============= ============= ============= =======================================

    """  
    r,param = _parsparam(r,param,npar=3)
    N  = param[0] # number of residues
    nu = param[1] # scaling exponent
    R0 = param[2] # residue length

    rsq = 6*(R0*N**nu)**2 # mean square end-to-end distance from radius of gyration
    normFact = 3/(2*np.pi*rsq)**(3/2) # normalization prefactor
    ShellSurf = 4*np.pi*r**2 # spherical shell surface
    Gaussian = np.exp(-3*r**2/(2*rsq))
    P = normFact*ShellSurf*Gaussian
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Center','Radius'),
units = ('nm','nm'),
start = np.asarray([3, 0.5]),
lower = np.asarray([1, 0.1]),
upper = np.asarray([20, 5 ])) 
@docstring()
def dd_circle(r,param):    
    r"""
Semicircle distribution model

Notes
-----

**Model:**

This provides a `semi-circle distribution <https://en.wikipedia.org/wiki/Wigner_semicircle_distribution>`_, defined by
:math:`P(r) = 2\pi\sqrt{(r-r_0)^2/R^2+1}` for :math:`r_0-R\le r\le r_0+R` and zero otherwise.

============== ================= ============= ============= ============= =================================
Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`r_0`            3.0           1              20          Center (nm)
``param[1]``   :math:`R`              0.5          0.1              5          Radius (nm)
============== ================= ============= ============= ============= =================================
    """  
    r,param = _parsparam(r,param,npar=2)
    r0 = param[0]
    R = abs(param[1])
    dr = r - r0
    idx = abs(dr)<R
    P = np.zeros(len(r))
    P[idx] = 2/np.pi/R**2*np.sqrt(R**2 - dr[idx]**2)
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Center','FWHM'),
units = ('nm','nm'),
start = np.asarray([3, 0.5]),
lower = np.asarray([1, 0.1]),
upper = np.asarray([20, 5 ])) 
@docstring()
def dd_cos(r,param):    
    r"""
Raised-cosine parametric model

Notes
-----

**Model:**

This provides a `raised-cosine distribution <https://en.wikipedia.org/wiki/Raised_cosine_distribution>`_, defined by 
:math:`P(r) = \frac{1}{2w}\cos\left(\frac{r-r_0}{w}\pi\right)` for :math:`r_0-w \le r \le r_0+w`, and zero otherwise.

============== ================= ============= ============= ============= =================================
Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`r_0`            3.0           1              20          Center (nm)
``param[1]``   :math:`w`              0.5          0.1             5           FWHM (nm)
============== ================= ============= ============= ============= =================================
    """  
    r,param = _parsparam(r,param,npar=2)
    r0 = param[0]
    fwhm = param[1]

    phi = (r-r0)/fwhm*np.pi
    P = (1 + np.cos(phi))/2/fwhm
    P[(r<(r0-fwhm)) | (r>(r0+fwhm))] = 0
    P = _normalize(r,P)
    return P
# =================================================================


def _pb(r,R):
# =================================================================
    P = np.zeros(len(r))
    idx = (r >= 0) & (r <= 2*R)
    P[idx] = 3*r[idx]**5/(16*R**6) - 9*r[idx]**3/(4*R**4) + 3*r[idx]**2/(R**3)

    return P
# =================================================================


def _pbs(r,R1,R2):
# =================================================================
    P = np.zeros(len(r))
    # Case1
    idx = (r >= 0) & (r < np.minimum(2*R1,R2 - R1)) 
    P[idx] = 12*r[idx]**3*R1**2 - r[idx]**5

    # Case2
    idx = (r >= R2 - R1) & (r < 2*R1)
    P[idx] = 8*r[idx]**2*(R2**3 - R1**3) - 3*r[idx]*(R2**2 - R1**2)**2 - 6*r[idx]**3*(R2 - R1)*(R2 + R1)

    # Case3
    idx = (r >= 2*R1) & (r < R2 - R1)
    P[idx] = 16*r[idx]**2*R1**3

    # Case4
    idx = (r >= np.maximum(R2 - R1,2*R1)) & (r < R1 + R2)
    P[idx] = r[idx]**5 - 6*r[idx]**3*(R2**2 + R1**2) + 8*r[idx]**2*(R2**3 + R1**3) - 3*r[idx]*(R2**2 - R1**2)**2

    P = P*3/(16*R1**3*(R2**3 - R1**3))
    return P
# =================================================================

# =================================================================
@metadata(
parameters = ('Inner shell radius','Shell thickness'),
units = ('nm','nm'),
lower = np.asarray([0.1, 0.1]),
upper = np.asarray([20,  20 ]),
start = np.asarray([1.5, 0.5])) 
@docstring()
def dd_shell(r,param):    
    r"""
Uniform spherical shell

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_shell.png
    :width: 25%

:math:`P(r) = \left(R_2^6 P_\mathrm{B}(r|R_2) - R_1^6 P_\mathrm{B}(r|R_1) - 2(r_2^3 - r_1^3)P_\mathrm{BS}(r|R_1,R_2)\right)/(R_2^3 - R_1^3)^2`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

:math:`P_\mathrm{B}(r|R_i) = \begin{cases} \frac{3r^5}{16R_i^6} - \frac{9r^3}{4R_i^4} + \frac{3r^2}{R_i^3} \quad \text{for} \quad 0 \leq r < 2R_i \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w`

============== ============== ============= ============= ============= =======================================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =======================================
``param[0]``     :math:`R`       1.5            0.1           20         Inner shell radius (nm)
``param[1]``     :math:`w`       0.5            0.1           20         Shell thickness (nm)
============== ============== ============= ============= ============= =======================================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 

    """  
    r,param = _parsparam(r,param,npar=2)
    R1 = float(param[0])
    w = float(param[1])
    R2 = R1 + w

    P = np.zeros(len(r))
    P = R2**6*_pb(r,R2) - R1**6*_pb(r,R1) - 2*(R2**3 - R1**3)*_pbs(r,R1,R2)

    P = P/(R2**3 - R1**3)**2

    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Sphere radius','Distance to point'),
units = ('nm','nm'),
lower = np.asarray([0.1, 0.1]),
upper = np.asarray([20,  20 ]),
start = np.asarray([1.5, 3.5])) 
@docstring()
def dd_spherepoint(r,param):    
    r"""
One particle distanced from particles distributed on a sphere

Notes
-----

**Model:**

.. image:: ../images/model_scheme_dd_spherepoint.png
    :width: 25%

:math:`P(r) = \begin{cases} \frac{3r(R^2-(d-r)^2)}{4dR^3} \quad \text{for} \quad d-R \leq r < d+R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

============== ============== ============= ============= ============= =========================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``     :math:`R`       1.5            0.1            20        Sphere radius (nm)
``param[1]``     :math:`d`       3.5            0.1            20        Distance to point (nm)
============== ============== ============= ============= ============= =========================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """ 
    r,param = _parsparam(r,param,npar=2)
    R = float(param[0])
    d = float(param[1])
    P = np.zeros(len(r))
    idx = (r >= d - R) & (r<= d + R)
    P[idx] = 3*r[idx]*(R**2 - (d - r[idx])**2)/(4*d*R**3)
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Sphere radius',),
units = ('nm',),
lower = np.asarray([0.1]),
upper = np.asarray([20]),
start = np.asarray([2.5])) 
@docstring()
def dd_spheresurf(r,param):    
    r"""
Particles distributed on a sphere's surface

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_spheresurf.png
    :width: 25%

:math:`P(r) = \begin{cases} \frac{r}{2R^2} \quad \text{for} \quad 0 \leq r < 2R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

============== ============== ============= ============= ============= =========================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``     :math:`R`       2.5              0.1          20        Sphere radius (nm)
============== ============== ============= ============= ============= =========================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """ 
    r,param = _parsparam(r,param,npar=1)
    R = float(param[0])
    P = np.zeros(len(r))
    idx = (r >= 0) & (r<= 2*R)
    P[idx] = r[idx]/R**2
    P = _normalize(r,P)
    return P
#=================================================================


# =================================================================
@metadata(
parameters = ('Inner shell radius','1st Shell thickness','2nd Shell thickness'),
units = ('nm','nm','nm'),
lower = np.asarray([0.1, 0.1, 0.1]),
upper = np.asarray([20,  20,  20 ]),
start = np.asarray([1.5, 0.5, 0.5])) 
@docstring()
def dd_shellshell(r,param):    
    r"""
Uniform spherical shell inside another spherical shell

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_shellshell.png
    :width: 25%

:math:`P(r) = (R_1^3(R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) - R_1^3(R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - R_2^3(R_3^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_3))/((R_3^3 - R_2^3)(R_2^3 - R_1^3))`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w_1`

:math:`R_3 = R + w_1 + w_2`

============== ============== ============= ============= ============= =======================================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =======================================
``param[0]``     :math:`R`       1.5             0.1           20         Inner shell radius (nm)
``param[1]``     :math:`w_1`     0.5             0.1           20         1st Shell thickness (nm)
``param[2]``     :math:`w_2`     0.5             0.1           20         2nd Shell thickness (nm)
============== ============== ============= ============= ============= =======================================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
    r,param = _parsparam(r,param,npar=3)
    R1 = float(param[0])
    w1 = float(param[1])
    w2 = float(param[2])

    R2 = R1 + w1
    R3 = R2 + w2

    delta21 = R2**3 - R1**3
    q21 = delta21*_pbs(r,R1,R2)
    delta31 = R3**3 - R1**3
    q31 = delta31*_pbs(r,R1,R3)
    delta32 = R3**3 - R2**3
    q32 = delta32*_pbs(r,R2,R3)

    P = R1**3*q21 - R1**3*q31 + R2**3*q32
    P = P/(delta21*delta32)

    P = _normalize(r,P)

    return P
#=================================================================


# =================================================================
@metadata(
parameters = ('Sphere radius','Shell thickness'),
units = ('nm','nm'),
lower = np.asarray([0.1, 0.1]),
upper = np.asarray([20,  20]),
start = np.asarray([1.5, 0.5])) 
@docstring()
def dd_shellsphere(r,param):    
    r"""
Particles distributed on a sphere inside a spherical shell

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_sphereshell.png
    :width: 25%

:math:`P(r) = \frac{3}{16R_1^3(R_2^3 - R_1^3)}\begin{cases} 12r^3R_1^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_1,R_2 - R_1) \\ 8r^2(R_2^3 - R_1^3) - 3r(R_2^2 - R_1^2)^2 - 6r^3(R_2 - R_1)(R_2 + R_1) \quad \text{for} \quad R_2-R_1 \leq r < 2R_1 \\ 16r^2R_1^3 \quad \text{for} \quad 2R_1\leq r < R_2 - R_1  \\  r^5 - 6r^3(R_2^2 + R_1^2) + 8r^2(R_2^3 + R_1^3) - 3r(R_2^2 - R1_2)^2 \quad \text{for} \quad \max(R_2-R_1,2R_1) \leq r < R_1+R_2 \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

with 

:math:`R_1 = R`

:math:`R_2 = R + w_1`

============== ============== ============= ============= ============= =========================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``     :math:`R`       1.5             0.1            20         Sphere radius (nm)
``param[1]``     :math:`w`       0.5             0.1            20         Shell thickness (nm)
============== ============== ============= ============= ============= =========================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
    r,param = _parsparam(r,param,npar=2)
    R1 = float(param[0])
    w = float(param[1])
    R2 = R1 + w
    P = _pbs(r,R1,R2)

    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Sphere radius','1st Shell thickness','2nd Shell thickness','Shell-Shell separation'),
units = ('nm','nm','nm','nm'),
lower = np.asarray([0.1, 0.1, 0.1, 0.1]),
upper = np.asarray([20,  20, 20, 2]),
start = np.asarray([0.75, 1, 1, 0.5])) 
@docstring()
def dd_shellvoidshell(r,param):    
    r"""
Particles distributed on a spherical shell inside another spherical shell separated by a void 

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_shellvoidshell.png
    :width: 25%

:math:`P(r) = \left(R_1^3((R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - (R_4^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_4)) + R_2^3((R_4^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_4) - (R_3^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_3)) \right)/((R_4^3 - R_3^3)(R_2^3 - R_1^3))`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w_1`

:math:`R_3 = R + w_1 + d`

:math:`R_4 = R + w_1 + d + w_2`

============== ============== ============= ============= ============= =============================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =============================
``param[0]``     :math:`R`       0.75             0.1           20        Sphere radius (nm)
``param[1]``     :math:`w_1`     1.00             0.1           20        1st Shell thickness (nm)
``param[2]``     :math:`w_2`     1.00             0.1           20        2nd Shell thickness (nm)
``param[3]``     :math:`d`       0.50             0.1           20        Shell-Shell separation (nm)
============== ============== ============= ============= ============= =============================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
    r,param = _parsparam(r,param,npar=4)
    R1 = float(param[0])
    w1 = float(param[1])
    w2 = float(param[2])
    d  = float(param[3])

    R2 = R1 + w1
    R3 = R1 + w1 + w2
    R4 = R1 + w1 + w2 + d

    delta21 = R2**3 - R1**3
    delta31 = R3**3 - R1**3
    q31 = delta31*_pbs(r,R1,R3)
    delta32 = R3**3 - R2**3
    q32 = delta32*_pbs(r,R2,R3)
    delta41 = R4**3 - R1**3
    q41 = delta41*_pbs(r,R1,R4)
    delta42 = R4**3 - R2**3
    q42 = delta42*_pbs(r,R2,R4)
    delta43 = R4**3 - R3**3

    P = (R1**3*(q31 - q41) + R2**3*(q42 - q32))/(delta43*delta21)
    P = np.round(P,15)

    P = _normalize(r,P)

    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Sphere radius','Shell thickness','Shell-Sphere separation'),
units = ('nm','nm','nm'),
lower = np.asarray([0.1, 0.1, 0.1]),
upper = np.asarray([20, 20, 20]),
start = np.asarray([1.5, 1, 0.5])) 
@docstring()
def dd_shellvoidsphere(r,param):    
    r"""
Particles distributed on a sphere inside a spherical shell separated by a void 


Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_shellvoidsphere.png
    :width: 25%

:math:`P(r) = ((R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - (R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) )/(R_3^3 - R_2^3)`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + d`

:math:`R_3 = R + d + w`

============== ============== ============= ============= ============= ==============================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= ==============================
``param[0]``     :math:`R`       1.5            0.1            20        Sphere radius (nm)
``param[1]``     :math:`w`       1.0            0.1            20        Shell thickness (nm)
``param[2]``     :math:`d`       0.5            0.1            20        Shell-Sphere separation (nm)
============== ============== ============= ============= ============= ==============================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
    r,param = _parsparam(r,param,npar=3)
    R1 = float(param[0])
    w  = float(param[1])
    d  = float(param[2])

    R2 = R1 + w
    R3 = R1 + w + d

    delta21 = R2**3 - R1**3
    q21 = delta21*_pbs(r,R1,R2)
    delta31 = R3**3 - R1**3
    q31 = delta31*_pbs(r,R1,R3)
    delta32 = R3**3 - R2**3
    
    P = (q31 - q21)/delta32
    P = np.round(P,15)
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Sphere radius',),
units = ('nm',),
lower = np.asarray([0.1]),
upper = np.asarray([20]),
start = np.asarray([2.5])) 
@docstring()
def dd_sphere(r,param):    
    r"""
Particles distributed on a sphere

Notes
-----

**Model:**
    
.. image:: ../images/model_scheme_dd_sphere.png
    :width: 25%

:math:`P(r) = \begin{cases} \frac{3r^5}{16R^6} - \frac{9r^3}{4R^4} + \frac{3r^2}{R^3} \quad \text{for} \quad 0 \leq r < 2R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

============== ============== ============= ============= ============= =========================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``     :math:`R`       2.5            0.1            20         Sphere radius (nm)
============== ============== ============= ============= ============= =========================

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
    r,param = _parsparam(r,param,npar=1)
    R = float(param[0])
    P = _pb(r,R)
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Mode','Left width','Right width'),
units = ('nm','nm','nm'),
lower = np.asarray([1, 0.1, 0.1]),
upper = np.asarray([20,  20, 20]),
start = np.asarray([3.5, 1, 0.5]))
@docstring()
def dd_triangle(r,param):    
    r"""
Triangle distribution model

Notes
-----

**Model:**

This provides a simple triangular distribution.

============== ======================== ============= ============= ============= =========================
Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =========================
``param[0]``   :math:`r_0`                 3.5            1.0              20         Mode (nm)
``param[1]``   :math:`w_\mathrm{L}`        0.3            0.1              5          Left width (nm)
``param[2]``   :math:`w_\mathrm{R}`        0.3            0.1              5          Right width (nm)
============== ======================== ============= ============= ============= =========================
    """  
    r,param = _parsparam(r,param,npar=3)
    r0 = param[0]
    wL = abs(param[1])
    wR = abs(param[2])
    rL = r0 - wL
    rR = r0 + wR
    idxL = (r >= r0-wL) & (r <= r0)
    idxR = (r <= r0+wR) & (r >= r0)
    P = np.zeros(len(r))
    if wL>0:
        P[idxL] = (r[idxL]-rL)/wL/(wL+wR)
    if wR>0:
        P[idxR] = -(r[idxR]-rR)/wR/(wL+wR)
    P = _normalize(r,P)
    return P
# =================================================================


# =================================================================
@metadata(
parameters = ('Left edge','Right edge'),
units = ('nm','nm'),
lower = np.asarray([0.1, 0.2]),
upper = np.asarray([6, 20]),
start = np.asarray([2.5, 3]))
@docstring()
def dd_uniform(r,param):    
    r"""
Uniform distribution model

Notes
-----

**Model:**
    
This provides a simple uniform distribution.

============== ======================== ============= ============= ============= =========================
Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =========================
``param[0]``   :math:`r_\mathrm{L}`         2.5             0.1            6           Left edge (nm)
``param[1]``   :math:`r_\mathrm{R}`         3.0             0.2            20          Right edge (nm)
============== ======================== ============= ============= ============= =========================
    """  
    r,param = _parsparam(r,param,npar=2)
    rL = min(abs(param))
    rR = max(abs(param))
    P = np.zeros(len(r))
    P[(r>=rL) & (r<=rR)] = 1
    P = _normalize(r,P)
    return P
# =================================================================

# -----------------------------------------------------------------
def wlc(r,L,Lp):
    
    P = np.zeros(len(r))
    
    kappa = Lp/L
    rnorm = r/L
    crit  = kappa*(1 - rnorm)

    idx = crit>0.2
    rcrit = rnorm[idx]
    P[idx] = 2*kappa/(4*np.pi)*(np.pi**2*(-1)**2*np.exp(-kappa*np.pi**2*(1-rcrit)) + np.pi**2*4*(-1)**3*np.exp(-kappa*np.pi**2*4*(1-rcrit)))

    idx = (crit>0) & (crit<0.2)
    rcrit = rnorm[idx]
    P[idx] = kappa/(4*np.pi*2*np.sqrt(np.pi))*((1/(kappa*(1 - rcrit))**(3/2)*np.exp(-(1 - 1/2)**2/(kappa*(1 - rcrit)))*(4*((1 - 1/2)/np.sqrt(kappa*(1-rcrit)))**2-2)) + 1/(kappa*(1 - rcrit))**(3/2)*np.exp(-(2 - 1/2)**2/(kappa*(1 - rcrit)))*(4*((2 - 1/2)/np.sqrt(kappa*(1-rcrit)))**2-2))

    return P
# -----------------------------------------------------------------


# =================================================================
@metadata(
parameters = ('Contour length','Persistence length'),
units = ('nm','nm'),
lower = np.asarray([1.5, 2]),
upper = np.asarray([20, 100]),
start = np.asarray([3.7, 10]))
@docstring()
def dd_wormchain(r,param):    
    r"""
Worm-like chain model near the rigid limit

Notes
-----

**Model:**

============== ============ ============= ============= ============= =========================
Variable         Symbol     Start Value   Lower bound   Upper bound      Description
============== ============ ============= ============= ============= =========================
``param[0]``   :math:`L`      3.7            1.5            10         Contour length (nm)
``param[1]``   :math:`L_p`    10             2              100        Persistence length (nm)
============== ============ ============= ============= ============= =========================

References
----------
.. [1] J. Wilhelm, E. Frey,
    Radial Distribution Function of Semiflexible Polymers
    Phys. Rev. Lett. 77(12), 2581-2584, 1996
    """  
    r,param = _parsparam(r,param,npar=2)
    L = param[0]
    Lp = param[1]
    P = wlc(r,L,Lp)
    P = _normalize(r,P)
    return P
#=================================================================


# =================================================================
@metadata(
parameters = ('Contour length','Persistence length','Gaussian standard deviation'),
units = ('nm','nm','nm'),
lower = np.asarray([1.5, 2, 0.001]),
upper = np.asarray([20, 100, 2]),
start = np.asarray([3.7, 10, 0.2]))
@docstring()
def dd_wormgauss(r,param):    
    r"""
Worm-like chain model near the rigid limit with Gaussian convolution

Notes
-----

**Model:**

============== ============== ============= ============= ============= ==================================
Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= ==================================
``param[0]``   :math:`L`         3.7            1.5           10         Contour length (nm)
``param[1]``   :math:`L_p`       10              2            100        Persistence length (nm)
``param[2]``   :math:`\sigma`    0.2            0.01           2         Gaussian standard deviation (nm)
============== ============== ============= ============= ============= ==================================

References
----------
.. [1] J. Wilhelm, E. Frey,
    Radial Distribution Function of Semiflexible Polymers
    Phys. Rev. Lett. 77(12), 2581-2584, 1996
    
    """  
    r,param = _parsparam(r,param,npar=3)
    L = param[0]
    Lp = param[1]
    sigma = param[2]
    P = wlc(r,L,Lp)

    # Compute Gaussian convolution window
    idx = np.argmax(P)
    gauss = np.exp(-((r - r[idx])/sigma)**2)
    
    # Convolution with size retention
    P = np.convolve(gauss,P,mode='full')

    # Adjust new convoluted axis
    idxconv = np.argmax(P)
    rconv = np.linspace(min(r),max(r)*2,len(P))
    rconv = rconv - abs(r[idx] - rconv[idxconv])

    #Interpolate down to original axis
    P = np.interp(r,rconv,P)
    P = _normalize(r,P)

    return P
#=================================================================
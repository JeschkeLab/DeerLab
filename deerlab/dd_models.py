# dd_models.py - Distance distribtion parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.


import inspect
import numpy as np
import scipy.special as spc
from deerlab.model import Model
from deerlab.utils import formatted_table

def _dd_docstring(model,notes):
#---------------------------------------------------------------------------------------
    args = model._parameter_list(order='vector')
    args.insert(model._constantsInfo[0]['argidx'],model._constantsInfo[0]['argkey'])

    parameters = ''
    for arg in args:
        if arg==model._constantsInfo[0]['argkey']:
            type = 'array_like'
            parameters += f'\n    {arg} : {type} \n        Distance axis, in nanometers.'
        elif len(np.atleast_1d(getattr(model,arg).idx))>1:
            type = 'array_like'
            parameters += f'\n    {arg} : {type} \n        {str(getattr(model,arg).description):s}.'
        else: 
            type = 'scalar'
            parameters += f'\n    {arg} : {type} \n        {str(getattr(model,arg).description):s}.'

        
    string = inspect.cleandoc(f"""

    {model.description}

    Parameters
    ----------
    {parameters}
    
    Returns
    -------

    P : ndarray
        Distance distribution.

    Notes
    -----

    **Parameter List**
    """) 
    string += '\n'
    string += '\n'
    table = []
    table.append(['Name','Lower','Upper','Type','Frozen','Unit','Description'])  
    for n,paramname in enumerate(model._parameter_list(order='vector')): 
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
    r = np.linspace(1,6,400)
    par0 = model.getmetadata()['par0']
    P = model(r,*par0)
    plt.figure(figsize=[6,3])
    plt.plot(r,P,color='#4550e6')
    plt.xlabel('r (nm)',fontsize=13)
    plt.ylabel('P (nm⁻¹)',fontsize=13)
    plt.grid(alpha=0.4)
    plt.tick_params(labelsize=12)
    plt.tick_params(labelsize=12)
    plt.tight_layout()
    plt.show()
"""
# =================================================================

# =================================================================
def _normalize(r,P):
    if not all(P==0):
        P = P/np.trapz(P,r)
    return P
# =================================================================

# =================================================================
def _multigaussfun(r,r0,sig):
    "Compute a distribution with multiple Gaussians"    
    r,r0,sig = np.atleast_2d(r,r0,sig)
    P = np.sqrt(1/(2*np.pi))*1/sig*np.exp(-0.5*((r.T-r0)/sig)**2)
    if not np.all(P==0):
        # Normalization
        P = np.squeeze(P)/np.sum([np.trapz(c,r) for c in P.T])
    else: 
        P = np.squeeze(P)
    return P
# =================================================================

def _multirice3dfun(r,nu,sig):
# =================================================================
    "Compute a distribution with multiple 3D-Rice distributions"    
    r,r0,sig = np.atleast_2d(r,nu,sig)
    nu = np.maximum(nu,0) # to avoid invalid values
    n = 3 # degrees of freedom
    r = r.T
    s2 = sig**2
    I_scaled = spc.ive(n/2-1, nu*r/s2)
    P = nu**(n/2-1)/s2*r**(n/2)*np.exp(-(r**2+nu**2)/(2*s2)+nu*r/s2)*I_scaled
    P[P<0] = 0
    
    # Normalization
    P = np.squeeze(P)/np.sum([np.trapz(c,np.squeeze(r)) for c in P.T])
    return P
# =================================================================

def freedist(r):
    def _nonparametric():
        return np.eye(len(r))
    # Create model
    dd_nonparametric = Model(_nonparametric,constants='r')
    dd_nonparametric.description = 'Non-parametric distribution model'
    # Parameters
    dd_nonparametric.addlinear('P',vec=len(r),lb=0,par0=0,description='Non-parametric distance distribution',unit="nm⁻¹", normalization=lambda P: _normalize(r,P))
    return dd_nonparametric


#=======================================================================================
#                                     dd_gauss
#=======================================================================================
notes = r"""
**Model**

:math:`P(r) = \frac{1}{\sigma\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r\right>)^2}{2\sigma^2}\right)`

where `\left<r\right>` is the mean distance and `\sigma` the standard deviation.
"""  
def _gauss(r,mean,std):
    return _multigaussfun(r,mean,std)
# Create model
dd_gauss = Model(_gauss,constants='r')
dd_gauss.description = 'Gaussian distribution model'
# Parameters
dd_gauss.mean.set(description='Mean', lb=1.0, ub=20, par0=3.5, unit='nm')
dd_gauss.std.set(description='Standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
# Add documentation
dd_gauss.__doc__ = _dd_docstring(dd_gauss,notes) + docstr_example('dd_gauss')



#=======================================================================================
#                                     dd_gauss2
#=======================================================================================
notes = r"""
**Model**

:math:`P(r) = a_1\frac{1}{\sigma_1\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r_1\right>)^2}{2\sigma_1^2}\right) + a_2\frac{1}{\sigma_2\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r_2\right>)^2}{2\sigma_2^2}\right)`

where `\left<r\right>_i` are the mean distances, `\sigma_i` the standard deviations, and `a_i` are the amplitudes of the Gaussians.
"""
def _gauss2(r,mean1,std1,mean2,std2):
    return _multigaussfun(r,[mean1,mean2],[std1,std2])
# Create model
dd_gauss2 = Model(_gauss2,constants='r')
dd_gauss2.description = 'Sum of two Gaussian distributions model'
# Parameters
dd_gauss2.mean1.set(description='1st Gaussian mean', lb=1.0, ub=20, par0=2.5, unit='nm')
dd_gauss2.mean2.set(description='2nd Gaussian mean', lb=1.0, ub=20, par0=4.5, unit='nm')
dd_gauss2.std1.set(description='1st Gaussian standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_gauss2.std2.set(description='2nd Gaussian standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_gauss2.addlinear('amp1',description='1st Gaussian amplitude', lb=0, par0=1, unit='')
dd_gauss2.addlinear('amp2',description='2nd Gaussian amplitude', lb=0, par0=1, unit='')
# Add documentation
dd_gauss2.__doc__ = _dd_docstring(dd_gauss2,notes) + docstr_example('dd_gauss2')



#=======================================================================================
#                                     dd_gauss3
#=======================================================================================
ntoes = r"""
**Model**

:math:`P(r) = a_1\frac{1}{\sigma_1\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r_1\right>)^2}{2\sigma_1^2}\right) + a_2\frac{1}{\sigma_2\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r_2\right>)^2}{2\sigma_2^2}\right) + a_3\frac{1}{\sigma_3\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r_3\right>)^2}{2\sigma_3^2}\right)`

where `\left<r\right>_i` are the mean distances, `\sigma_i` the standard deviations, and `a_i` are the amplitudes of the Gaussians.
"""
def _gauss3(r,mean1,std1,mean2,std2,mean3,std3):
    return _multigaussfun(r,[mean1,mean2,mean3],[std1,std2,std3])
# Create model
dd_gauss3 = Model(_gauss3,constants='r')
dd_gauss3.description = 'Sum of three Gaussian distributions model'
# Parameters
dd_gauss3.mean1.set(description='1st Gaussian mean', lb=1.0, ub=20, par0=2.5, unit='nm')
dd_gauss3.mean2.set(description='2nd Gaussian mean', lb=1.0, ub=20, par0=3.5, unit='nm')
dd_gauss3.mean3.set(description='3rd Gaussian mean', lb=1.0, ub=20, par0=5.0, unit='nm')
dd_gauss3.std1.set(description='1st Gaussian standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_gauss3.std2.set(description='2nd Gaussian standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_gauss3.std3.set(description='3rd Gaussian standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_gauss3.addlinear('amp1',description='1st Gaussian amplitude', lb=0, par0=1, unit='')
dd_gauss3.addlinear('amp2',description='2nd Gaussian amplitude', lb=0, par0=1, unit='')
dd_gauss3.addlinear('amp3',description='3rd Gaussian amplitude', lb=0, par0=1, unit='')
# Add documentation
dd_gauss3.__doc__ = _dd_docstring(dd_gauss3,notes) +  docstr_example('dd_gauss3')



#=======================================================================================
#                                     dd_gengauss
#=======================================================================================
notes =  r"""
**Model**

.. raw:: html 
    
    <img src="../_images/model_scheme_dd_gengauss.png", style="width: 50%">
    <br><br><br>    

:math:`P(r) = \frac{\beta}{2\sigma\Gamma(1/\beta)}\exp\left(-\left(\frac{(r-\left<r\right>)}{\sigma}\right)^\beta \right)`

where `\left<r\right>` is the mean distance,`\sigma` is the standard deviation, and `\beta` is the kurtosis of the distribution.
"""  
def _gengauss(r,mean,std,kurt):
    P = kurt/(2*std*spc.gamma(1/kurt))*np.exp(-abs(r-mean)/std**kurt)
    return _normalize(r,P) 
# Create model
dd_gengauss = Model(_gengauss,constants='r')
dd_gengauss.description = 'Generalized Gaussian distribution model'
# Parameters
dd_gengauss.mean.set(description='Mean', lb=1.0, ub=20, par0=3.5, unit='nm')
dd_gengauss.std.set(description='Standard deviation', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_gengauss.kurt.set(description='Kurtosis', lb=0.25, ub=15, par0=5.0, unit='')
# Add documentation
dd_gengauss.__doc__ = _dd_docstring(dd_gengauss,notes) +  docstr_example('dd_gengauss')

    

#=======================================================================================
#                                     dd_skewgauss
#=======================================================================================
notes = r"""
**Model**

.. raw:: html 
    
    <img src="../_images/model_scheme_dd_skewgauss.png", style="width: 50%">
    <br><br><br>  

:math:`P(r) = \frac{1}{\sqrt{2\pi}}\exp\left(-\frac{(r-\left<r\right>)^2}{2\sigma^2}\right)\left(1 + \mathrm{erf}\left(\frac{(r-\left<r\right>)}{\sqrt{2}\sigma}\right) \right)`

where `\left<r\right>` is the center distance,`\sigma` is the spread, and `\alpha` is the skewness of the distribution.
"""
def _skewgauss(r,center,std,skew):
    x = (r-center)/std
    P = 1/np.sqrt(2*np.pi)*np.exp(-x**2/2)*(1 + spc.erf(skew*x/np.sqrt(2)))
    return _normalize(r,P) 
# Create model
dd_skewgauss = Model(_skewgauss,constants='r')
dd_skewgauss.description = 'Skew Gaussian distribution model'
# Parameters
dd_skewgauss.center.set(description='Center', lb=1.0, ub=20, par0=3.5, unit='nm')
dd_skewgauss.std.set(description='Spread', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_skewgauss.skew.set(description='Skewness', lb=-25, ub=25, par0=5, unit='')
# Add documentation
dd_skewgauss.__doc__ = _dd_docstring(dd_skewgauss,notes) +  docstr_example('dd_skewgauss')



#=======================================================================================
#                                     dd_rice
#=======================================================================================
notes = r"""
**Model** 

:math:`P(r) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`\nu` is the location of the distribution, :math:`\sigma` its spread, :math:`n=3`, and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`. This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.
"""  
def _rice(r,location,spread): 
    return _multirice3dfun(r,[location],[spread])
# Create model
dd_rice = Model(_rice,constants='r')
dd_rice.description = '3D-Rice distribution model'
# Parameters
dd_rice.location.set(description='Location', lb=1.0, ub=20, par0=3.5, unit='nm')
dd_rice.spread.set(description='Spread', lb=0.1, ub=5, par0=0.7, unit='nm')
# Add documentation
dd_rice.__doc__ = _dd_docstring(dd_rice,notes) +  docstr_example('dd_rice')
    


#=======================================================================================
#                                     dd_rice2
#=======================================================================================
notes = r"""
**Model**

:math:`P(r) = a_1 R(r,\nu_1,\sigma_1) + a_2 R(r,\nu_2,\sigma_2)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.
"""
def _rice2(r,location1,spread1,location2,spread2): 
    return _multirice3dfun(r,[location1,location2],[spread1,spread2])
# Create model
dd_rice2 = Model(_rice2,constants='r')
dd_rice2.description = 'Sum of two 3D-Rice distributions model'
# Parameters
dd_rice2.location1.set(description='1st Rician location', lb=1.0, ub=20, par0=2.5, unit='nm')
dd_rice2.location2.set(description='2nd Rician location', lb=1.0, ub=20, par0=4.5, unit='nm')
dd_rice2.spread1.set(description='1st Rician spread', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_rice2.spread2.set(description='2nd Rician spread', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_rice2.addlinear('amp1',description='1st Rician amplitude', lb=0, par0=1, unit='')
dd_rice2.addlinear('amp2',description='2nd Rician amplitude', lb=0, par0=1, unit='')
# Add documentation
dd_rice2.__doc__ = _dd_docstring(dd_rice2,notes) +  docstr_example('dd_rice2')
    


#=======================================================================================
#                                     dd_rice3
#=======================================================================================
notes =  r"""
**Model**
    
:math:`P(r) = a_1 R(r,\nu_1,\sigma_1) + a_2 R(r,\nu_2,\sigma_2) + a_3 R(r,\nu_3,\sigma_3)`

:math:`R(r,\nu,\sigma) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`.
This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.
"""
def _rice3(r,location1,spread1,location2,spread2,location3,spread3): 
    return _multirice3dfun(r,[location1,location2,location3],[spread1,spread2,spread3])
# Create model
dd_rice3 = Model(_rice3,constants='r')
dd_rice3.description = 'Sum of two 3D-Rice distributions model'
# Parameters
dd_rice3.location1.set(description='1st Rician location', lb=1.0, ub=20, par0=2.5, unit='nm')
dd_rice3.location2.set(description='2nd Rician location', lb=1.0, ub=20, par0=4.5, unit='nm')
dd_rice3.location3.set(description='3rd Rician location', lb=1.0, ub=20, par0=5.5, unit='nm')
dd_rice3.spread1.set(description='1st Rician spread', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_rice3.spread2.set(description='2nd Rician spread', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_rice3.spread3.set(description='3rd Rician spread', lb=0.05, ub=2.5, par0=0.2, unit='nm')
dd_rice3.addlinear('amp2',description='2nd Rician amplitude', lb=0, par0=1,  unit='')
dd_rice3.addlinear('amp1',description='1st Rician amplitude', lb=0, par0=1, unit='')
dd_rice3.addlinear('amp3',description='3rd Rician amplitude', lb=0, par0=1, unit='')
# Add documentation
dd_rice3.__doc__ = _dd_docstring(dd_rice3,notes) +  docstr_example('dd_rice3')



#=======================================================================================
#                                     dd_randcoil
#=======================================================================================
notes = r"""
**Model**

.. raw:: html 
    
    <img src="../_images/model_scheme_dd_randcoil.png", style="width: 50%">
    <br><br><br>

:math:`P(r) = \frac{3}{(2\pi\nu_0)^{3/2}}4\pi r^2\exp(-\frac{3 r^2}{\nu_0})`

where `\nu_0 = 3/(12\pi r_0 N \nu)^{3/2}`, and `N` is the number of residues, `R_0` 
is the segment length, and `\nu` is the scaling exponent.
"""  
def _randcoil(r,Nres,scaling,length):
    rsq = 6*(length*Nres**scaling)**2 # mean square end-to-end distance from radius of gyration
    normFact = 3/(2*np.pi*rsq)**(3/2) # normalization prefactor
    ShellSurf = 4*np.pi*r**2 # spherical shell surface
    Gaussian = np.exp(-3*r**2/(2*rsq))
    P = normFact*ShellSurf*Gaussian
    P = _normalize(r,P)
    return P
# Create model
dd_randcoil = Model(_randcoil,constants='r')
dd_randcoil.description = 'Random-coil model for an unfolded peptide/protein'
# Parameters
dd_randcoil.Nres.set(description='Number of residues', lb=2.0, ub=1000, par0=50, unit='')
dd_randcoil.scaling.set(description='Segment length', lb=0.1, ub=0.4, par0=0.2, unit='nm')
dd_randcoil.length.set(description='Scaling exponent', lb=0.33, ub=1.00, par0=0.602, unit='')
# Add documentation
dd_randcoil.__doc__ = _dd_docstring(dd_randcoil,notes) +  docstr_example('dd_randcoil')




#=======================================================================================
#                                     dd_circle
#=======================================================================================
notes = r"""
**Model**

This provides a `semi-circle distribution <https://en.wikipedia.org/wiki/Wigner_semicircle_distribution>`_, defined by

:math:`P(r) = 2\pi\sqrt{(r-r_0)^2/R^2+1}` for :math:`r_0-R\le r\le r_0+R` and zero otherwise.
"""  
def _circle(r,center,radius):
    dr = r - center
    idx = abs(dr)<radius
    P = np.zeros(len(r))
    P[idx] = 2/np.pi/radius**2*np.sqrt(radius**2 - dr[idx]**2)
    P = _normalize(r,P)
    return P
# Create model
dd_circle = Model(_circle,constants='r')
dd_circle.description = 'Semicircle distribution model'
# Parameters
dd_circle.center.set(description='Center', lb=1, ub=20, par0=3, unit='nm')
dd_circle.radius.set(description='Radius', lb=0.1, ub=5, par0=0.5, unit='nm')
# Add documentation
dd_circle.__doc__ = _dd_docstring(dd_circle,notes) +  docstr_example('dd_circle')




#=======================================================================================
#                                     dd_cos
#=======================================================================================
notes = r"""
**Model**

This provides a `raised-cosine distribution <https://en.wikipedia.org/wiki/Raised_cosine_distribution>`_, defined by 
:math:`P(r) = \frac{1}{2w}\cos\left(\frac{r-r_0}{w}\pi\right)` for :math:`r_0-w \le r \le r_0+w`, and zero otherwise.
"""  
def _rcos(r,center,fwhm):

    phi = (r-center)/fwhm*np.pi
    P = (1 + np.cos(phi))/2/fwhm
    P[(r<(center-fwhm)) | (r>(center+fwhm))] = 0
    P = _normalize(r,P)
    return P
# Create model
dd_cos = Model(_rcos,constants='r')
dd_cos.description = 'Raised-cosine parametric model'
# Parameters
dd_cos.center.set(description='Center', lb=1, ub=20, par0=3, unit='nm')
dd_cos.fwhm.set(description='FWHM', lb=0.1, ub=5, par0=0.5, unit='nm')
# Add documentation
dd_cos.__doc__ = _dd_docstring(dd_cos,notes) +  docstr_example('dd_cos')


#------------------------------------------------------------------
def _pb(r,R):
    P = np.zeros(len(r))
    idx = (r >= 0) & (r <= 2*R)
    P[idx] = 3*r[idx]**5/(16*R**6) - 9*r[idx]**3/(4*R**4) + 3*r[idx]**2/(R**3)

    return P
#------------------------------------------------------------------
def _pbs(r,R1,R2):
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
#------------------------------------------------------------------


#=======================================================================================
#                                     dd_shell
#=======================================================================================
notes = r"""
**Model**
    
.. raw:: html 
    
    <img src="../_images/model_scheme_dd_shell.png", style="width: 50%">
    <br><br><br>

:math:`P(r) = \left(R_2^6 P_\mathrm{B}(r|R_2) - R_1^6 P_\mathrm{B}(r|R_1) - 2(r_2^3 - r_1^3)P_\mathrm{BS}(r|R_1,R_2)\right)/(R_2^3 - R_1^3)^2`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

:math:`P_\mathrm{B}(r|R_i) = \begin{cases} \frac{3r^5}{16R_i^6} - \frac{9r^3}{4R_i^4} + \frac{3r^2}{R_i^3} \quad \text{for} \quad 0 \leq r < 2R_i \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w`

where :math:`R` is the inner shell radius, and :math:`w` is the shell thickness.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
"""  
def _shell(r,radius,thickness):
    R1 = float(radius)
    w = float(thickness)
    R2 = R1 + w
    P = np.zeros(len(r))
    P = R2**6*_pb(r,R2) - R1**6*_pb(r,R1) - 2*(R2**3 - R1**3)*_pbs(r,R1,R2)
    P = P/(R2**3 - R1**3)**2
    P = _normalize(r,P)
    return P
# Create model
dd_shell = Model(_shell,constants='r')
dd_shell.description = 'Uniform distribution of particles on a spherical shell'
# Parameters
dd_shell.radius.set(description='Inner shell radius', lb=0.1, ub=20, par0=1.5, unit='nm')
dd_shell.thickness.set(description='Shell thickness', lb=0.1, ub=20, par0=0.5, unit='nm')
# Add documentation
dd_shell.__doc__ = _dd_docstring(dd_shell,notes) +  docstr_example('dd_shell')




#=======================================================================================
#                                     dd_spherepoint
#=======================================================================================
notes = r"""
**Model**

.. raw:: html 
    
    <img src="../_images/model_scheme_dd_spherepoint.png", style="width: 50%">
    <br><br><br>

:math:`P(r) = \begin{cases} \frac{3r(R^2-(d-r)^2)}{4dR^3} \quad \text{for} \quad d-R \leq r < d+R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

where :math:`R` is the sphere's radius, and :math:`d` is the distance from the sphere to the point.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
""" 
def _spherepoint(r,radius,dist):
    R = float(radius)
    d = float(dist)
    P = np.zeros(len(r))
    idx = (r >= d - R) & (r<= d + R)
    P[idx] = 3*r[idx]*(R**2 - (d - r[idx])**2)/(4*d*R**3)
    P = _normalize(r,P)
    return P
# Create model
dd_spherepoint = Model(_spherepoint,constants='r')
dd_spherepoint.description = 'One particle distanced from particles uniformly distributed on a sphere'
# Parameters
dd_spherepoint.radius.set(description='Sphere radius', lb=0.1, ub=20, par0=1.5, unit='nm')
dd_spherepoint.dist.set(description='Distance to point', lb=0.1, ub=20, par0=3.5, unit='nm')
# Add documentation
dd_spherepoint.__doc__ = _dd_docstring(dd_spherepoint,notes) +  docstr_example('dd_spherepoint')




#=======================================================================================
#                                     dd_spheresurf
#=======================================================================================
notes = r"""
**Model**

.. raw:: html 
    
    <img src="../_images/model_scheme_dd_spheresurf.png", style="width: 50%">
    <br><br><br>

:math:`P(r) = \begin{cases} \frac{r}{2R^2} \quad \text{for} \quad 0 \leq r < 2R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

where :math:`R` is the sphere's radius.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
""" 
def _spheresurf(r,radius):
    R = float(radius)
    P = np.zeros(len(r))
    idx = (r >= 0) & (r<= 2*R)
    P[idx] = r[idx]/R**2
    P = _normalize(r,P)
    return P
# Create model
dd_spheresurf = Model(_spheresurf,constants='r')
dd_spheresurf.description = "Particles uniformly distributed on a sphere's surface."
# Parameters
dd_spheresurf.radius.set(description='Sphere radius', lb=0.1, ub=20, par0=2.5, unit='nm')
# Add documentation
dd_spheresurf.__doc__ = _dd_docstring(dd_spheresurf,notes) +  docstr_example('dd_spheresurf')



#=======================================================================================
#                                     dd_shellshell
#=======================================================================================
notes = r"""
**Model**

.. raw:: html 
    
    <img src="../_images/model_scheme_dd_shellshell.png", style="width: 50%">
    <br><br><br>    

:math:`P(r) = (R_1^3(R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) - R_1^3(R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - R_2^3(R_3^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_3))/((R_3^3 - R_2^3)(R_2^3 - R_1^3))`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w_1`

:math:`R_3 = R + w_1 + w_2`

where :math:`R` is the inner shell's radius, and :math:`w1` and :math:`w2` denote the thickness of the inner and outer shells, respectively.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
def _shellshell(r,radius,thickness1,thickness2):
    R1 = float(radius)
    w1 = float(thickness1)
    w2 = float(thickness2)
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
# Create model
dd_shellshell = Model(_shellshell,constants='r')
dd_shellshell.description = 'Particles uniformly distributed on a spherical shell and on another concentric spherical shell.'
# Parameters
dd_shellshell.radius.set(description='Inner shell radius', lb=0.1, ub=20, par0=1.5, unit='nm')
dd_shellshell.thickness1.set(description='Inner shell thickness', lb=0.1, ub=20, par0=0.5, unit='nm')
dd_shellshell.thickness2.set(description='Outer shell thickness', lb=0.1, ub=20, par0=0.5, unit='nm')
# Add documentation
dd_shellshell.__doc__ = _dd_docstring(dd_shellshell,notes) +  docstr_example('dd_shellshell')


#=======================================================================================
#                                     dd_shellsphere
#=======================================================================================
notes = r"""
**Model**
    
.. raw:: html 
    
    <img src="../_images/model_scheme_dd_sphereshell.png", style="width: 50%">
    <br><br><br>     

:math:`P(r) = \frac{3}{16R_1^3(R_2^3 - R_1^3)}\begin{cases} 12r^3R_1^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_1,R_2 - R_1) \\ 8r^2(R_2^3 - R_1^3) - 3r(R_2^2 - R_1^2)^2 - 6r^3(R_2 - R_1)(R_2 + R_1) \quad \text{for} \quad R_2-R_1 \leq r < 2R_1 \\ 16r^2R_1^3 \quad \text{for} \quad 2R_1\leq r < R_2 - R_1  \\  r^5 - 6r^3(R_2^2 + R_1^2) + 8r^2(R_2^3 + R_1^3) - 3r(R_2^2 - R1_2)^2 \quad \text{for} \quad \max(R_2-R_1,2R_1) \leq r < R_1+R_2 \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

with 

:math:`R_1 = R`

:math:`R_2 = R + w_1`

where :math:`R` is the inner sphere's radius, and :math:`w` is the outer shell's thickness.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
def _shellsphere(r,radius,thickness):
    R1 = float(radius)
    w = float(thickness)
    R2 = R1 + w
    P = _pbs(r,R1,R2)
    P = _normalize(r,P)
    return P
# Create model
dd_shellsphere = Model(_shellsphere,constants='r')
dd_shellsphere.description = 'Particles uniformly distributed on a sphere and on an outer spherical shell.'
# Parameters
dd_shellsphere.radius.set(description='Inner shell radius', lb=0.1, ub=20, par0=1.5, unit='nm')
dd_shellsphere.thickness.set(description='Inner shell thickness', lb=0.1, ub=20, par0=0.5, unit='nm')
# Add documentation
dd_shellsphere.__doc__ = _dd_docstring(dd_shellsphere,notes) +  docstr_example('dd_shellsphere')



#=======================================================================================
#                                     dd_shellvoidshell
#=======================================================================================
notes = r"""
**Model**
    
.. raw:: html 
    
    <img src="../_images/model_scheme_dd_shellvoidshell.png", style="width: 50%">
    <br><br><br>    

:math:`P(r) = \left(R_1^3((R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - (R_4^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_4)) + R_2^3((R_4^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_4) - (R_3^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_3)) \right)/((R_4^3 - R_3^3)(R_2^3 - R_1^3))`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w_1`

:math:`R_3 = R + w_1 + d`

:math:`R_4 = R + w_1 + d + w_2`

where :math:`R` is the inner shell's radius, :math:`w_1` and :math:`w_2` denote the thickness of the inner and outer shells, respectively, and :math:`d` is
the shell-shell separation.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
def _shellvoidshell(r,radius,thickness1,thickness2,separation):
    R1 = float(radius)
    w1 = float(thickness1)
    w2 = float(thickness2)
    d  = float(separation)
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
# Create model
dd_shellvoidshell = Model(_shellvoidshell,constants='r')
dd_shellvoidshell.description = 'Particles uniformly distributed on a spherical shell and on another concentric spherical shell separated by a void.'
# Parameters
dd_shellvoidshell.radius.set(description='Inner shell radius', lb=0.1, ub=20, par0=0.75, unit='nm')
dd_shellvoidshell.thickness1.set(description='Inner shell thickness', lb=0.1, ub=20, par0=1.0, unit='nm')
dd_shellvoidshell.thickness2.set(description='Outer shell thickness', lb=0.1, ub=20, par0=1.0, unit='nm')
dd_shellvoidshell.separation.set(description='Shell-shell separation', lb=0.1, ub=20, par0=0.5, unit='nm')
# Add documentation
dd_shellvoidshell.__doc__ = _dd_docstring(dd_shellvoidshell,notes) +  docstr_example('dd_shellvoidshell')




#=======================================================================================
#                                     dd_shellvoidsphere
#=======================================================================================
notes = r"""
**Model**
    
.. raw:: html 
    
    <img src="../_images/model_scheme_dd_shellvoidsphere.png", style="width: 50%">
    <br><br><br>    


:math:`P(r) = ((R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - (R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) )/(R_3^3 - R_2^3)`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + d`

:math:`R_3 = R + d + w`

where :math:`R` is the inner sphere's radius, :math:`w` is the thickness of the outer shell, and :math:`d` is the shell-sphere separation.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
def _shellvoidsphere(r,radius,thickness,separation):
    R1 = float(radius)
    w  = float(thickness)
    d  = float(separation)

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
# Create model
dd_shellvoidsphere = Model(_shellvoidsphere,constants='r')
dd_shellvoidsphere.description = 'Particles uniformly distributed on a sphere and on a concentric outer spherical shell separated by a void.'
# Parameters
dd_shellvoidsphere.radius.set(description='Sphere radius', lb=0.1, ub=20, par0=1.5, unit='nm')
dd_shellvoidsphere.thickness.set(description='Outer shell thickness', lb=0.1, ub=20, par0=1.0, unit='nm')
dd_shellvoidsphere.separation.set(description='Shell-sphere separation', lb=0.1, ub=20, par0=0.5, unit='nm')
# Add documentation
dd_shellvoidsphere.__doc__ = _dd_docstring(dd_shellvoidsphere,notes) +  docstr_example('dd_shellvoidsphere')





#=======================================================================================
#                                     dd_sphere
#=======================================================================================
notes =  r"""
**Model**
    
.. raw:: html 
    
    <img src="../_images/model_scheme_dd_sphere.png", style="width: 50%">
    <br><br><br>    

:math:`P(r) = \begin{cases} \frac{3r^5}{16R^6} - \frac{9r^3}{4R^4} + \frac{3r^2}{R^3} \quad \text{for} \quad 0 \leq r < 2R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

where :math:`R` is the sphere's radius.

References
----------
.. [1] D.R. Kattnig, D. Hinderberger,
    Analytical distance distributions in systems of spherical symmetry with applications to double electron-electron resonance, JMR, 230, 50-63, 2013 
    """  
def _sphere(r,radius):
    P = _pb(r,radius)
    P = _normalize(r,P)
    return P
# Create model
dd_sphere = Model(_sphere,constants='r')
dd_sphere.description = 'Particles uniformly distributed on a sphere.'
# Parameters
dd_sphere.radius.set(description='Sphere radius', lb=0.1, ub=20, par0=2.5, unit='nm')
# Add documentation
dd_sphere.__doc__ = _dd_docstring(dd_sphere,notes) +  docstr_example('dd_sphere')




#=======================================================================================
#                                     dd_triangle
#=======================================================================================
notes = r"""
**Model**

This provides a simple `triangular distribution <https://en.wikipedia.org/wiki/Triangular_distribution>`_.

"""  
def _triangle(r,mode,left,right):
    r0 = mode
    wL = abs(left)
    wR = abs(right)
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
# Create model
dd_triangle = Model(_triangle,constants='r')
dd_triangle.description = 'Triangular distribution model.'
# Parameters
dd_triangle.mode.set(description='Mode', lb=1, ub=20, par0=3.5, unit='nm')
dd_triangle.left.set(description='Left width', lb=0.1, ub=5, par0=0.3, unit='nm')
dd_triangle.right.set(description='Right width', lb=0.1, ub=5, par0=0.3, unit='nm')
# Add documentation
dd_triangle.__doc__ = _dd_docstring(dd_triangle,notes) +  docstr_example('dd_triangle')




#=======================================================================================
#                                     dd_uniform
#=======================================================================================
notes = r"""
**Model**
    
This provides a simple uniform distribution.
"""  
def _uniform(r,left,right):
    rL = np.min(np.abs([left,right]))
    rR = np.max(np.abs([left,right]))
    P = np.zeros(len(r))
    P[(r>=rL) & (r<=rR)] = 1
    P = _normalize(r,P)
    return P
# Create model
dd_uniform = Model(_uniform,constants='r')
dd_uniform.description = 'Uniform distribution model.'
# Parameters
dd_uniform.left.set(description='Left edge', lb=0.1, ub=6, par0=2.5, unit='nm')
dd_uniform.right.set(description='Right edge', lb=0.2, ub=20, par0=3.5, unit='nm')
# Add documentation
dd_uniform.__doc__ = _dd_docstring(dd_uniform,notes) +  docstr_example('dd_uniform')



# ----------------------------------------------------------------------------------
def _wlc(r,L,Lp):
    
    P = np.zeros(len(r))
    
    kappa = Lp/np.maximum(L,1e-16)
    rnorm = r/np.maximum(L,1e-16)
    crit  = kappa*(1 - rnorm)

    idx = crit>0.2
    rcrit = rnorm[idx]
    P[idx] = 2*kappa/(4*np.pi)*(np.pi**2*(-1)**2*np.exp(-kappa*np.pi**2*(1-rcrit)) + np.pi**2*4*(-1)**3*np.exp(-kappa*np.pi**2*4*(1-rcrit)))

    idx = (crit>0) & (crit<0.2)
    rcrit = rnorm[idx]
    P[idx] = kappa/(4*np.pi*2*np.sqrt(np.pi))*((1/(kappa*(1 - rcrit))**(3/2)*np.exp(-(1 - 1/2)**2/(kappa*(1 - rcrit)))*(4*((1 - 1/2)/np.sqrt(kappa*(1-rcrit)))**2-2)) + 1/(kappa*(1 - rcrit))**(3/2)*np.exp(-(2 - 1/2)**2/(kappa*(1 - rcrit)))*(4*((2 - 1/2)/np.sqrt(kappa*(1-rcrit)))**2-2))

    return P
# ----------------------------------------------------------------------------------


#=======================================================================================
#                                     dd_wormchain
#=======================================================================================
notes = r"""

References
----------
.. [1] J. Wilhelm, E. Frey,
    Radial Distribution Function of Semiflexible Polymers
    Phys. Rev. Lett. 77(12), 2581-2584, 1996
"""  
def _wormchain(r,contour,persistence):
    P = _wlc(r,contour,persistence)
    P = _normalize(r,P)
    return P
# Create model
dd_wormchain = Model(_wormchain,constants='r')
dd_wormchain.description = 'Worm-like chain model near the rigid limit.'
# Parameters
dd_wormchain.contour.set(description='Contour length', lb=1.5, ub=10, par0=3.7, unit='nm')
dd_wormchain.persistence.set(description='Persistence length', lb=2, ub=100, par0=10, unit='nm')
# Add documentation
dd_wormchain.__doc__ = _dd_docstring(dd_wormchain,notes) +  docstr_example('dd_wormchain')



#=======================================================================================
#                                     dd_wormgauss
#=======================================================================================
notes = r"""
References
----------
.. [1] J. Wilhelm, E. Frey,
    Radial Distribution Function of Semiflexible Polymers
    Phys. Rev. Lett. 77(12), 2581-2584, 1996
"""  
def _wormgauss(r,contour,persistence,std):
    P = _wlc(r,contour,persistence)
    # Compute Gaussian convolution window
    idx = np.argmax(P)
    gauss = np.exp(-((r - r[idx])/std)**2)
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
# Create model
dd_wormgauss = Model(_wormgauss,constants='r')
dd_wormgauss.description = 'Worm-like chain model near the rigid limit with Gaussian convolution.'
# Parameters
dd_wormgauss.contour.set(description='Contour length', lb=1.5, ub=10, par0=3.7, unit='nm')
dd_wormgauss.persistence.set(description='Persistence length', lb=2, ub=100, par0=10, unit='nm')
dd_wormgauss.std.set(description='Gaussian standard deviation', lb=0.01, ub=5, par0=0.2, unit='nm')
# Add documentation
dd_wormgauss.__doc__ = _dd_docstring(dd_wormgauss,notes) +  docstr_example('dd_wormgauss')

# dd_models.py - Distance distribtion parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import math as m
import numpy as np
import scipy.special as spc
import inspect

def _parsargs(args,npar):
#=================================================================
    name = inspect.stack()[1][3]
    if len(args)!=2:
        raise KeyError('The model function {} requires two input arguments: {}(r,params).'.format(name,name))
    r,p = args
    if len(p)!=npar:
        raise ValueError('The model function {} requires {} parameters, but {} are provided.'.format(name,npar,len(p)))
    return r,p
#=================================================================

def _normalize(r,P):
#=================================================================
    if not all(P==0):
        P = P/np.trapz(P,r)
    return P
#=================================================================

def _multigaussfun(r,r0,sig,a):
#=================================================================
    "Compute a distribution with multiple Gaussians"    
    n = len(r0)
    P = np.zeros_like(r)
    for k in range(n):
        P += a[k]*m.sqrt(1/(2*m.pi))*1/sig[k]*np.exp(-0.5*((r-r0[k])/sig[k])**2)
    P = _normalize(r,P)
    return P
#=================================================================

def _multirice3dfun(r,nu,sig,a):
#=================================================================
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
#=================================================================


def dd_gauss(*args):    
#=================================================================
    r"""
    Gaussian distribution
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_gauss()


    Otherwise the function returns to calculated distance distribution::

        P = dd_gauss(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -------------------------------------------------------------
      Parameter                 Units     Lower    Upper    Start
     -------------------------------------------------------------
      Mean                       nm        1        20       3.5 
      Standard deviation         nm       0.05      2.5      0.2 
     -------------------------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Mean','Standard deviation'),
            Units = ('nm','nm'),
            Start = np.asarray([3.5, 0.2]),
            Lower = np.asarray([1, 0.05]),
            Upper = np.asarray([20, 2.5])
        )
        return info
    r,p = _parsargs(args,npar=2)

    r0 = [p[0]]
    sigma = [p[1]]
    a = [1.0]
    P = _multigaussfun(r,r0,sigma,a)
    return P
#=================================================================
    

def dd_gauss2(*args):
#=================================================================
    r"""
    Sum of two Gaussian distributions
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_gauss2()


    Otherwise the function returns to calculated distance distribution::

        P = dd_gauss2(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
     -------------------------------------------------------------------------
      Parameter                             Units     Lower    Upper    Start
     -------------------------------------------------------------------------
      Mean of 1st Gaussian                   nm        1        20       2.5 
      Standard deviation of 1st Gaussian     nm       0.05      2.5      0.2 
      Amplitude of 1st Gaussian                        0        1        0.5 
      Mean of 2nd Gaussian                   nm        1        20       3.5 
      Standard deviation of 2nd Gaussian     nm       0.05      2.5      0.2 
      Amplitude of 2nd Gaussian                        0        1        0.5 
     -------------------------------------------------------------------------
    """

    if not args:
        info = dict(
            Parameters = ('Mean of 1st Gaussian', 'Standard deviation of 1st Gaussian', 'Amplitude of 1st Gaussian',
                          'Mean of 2nd Gaussian', 'Standard deviation of 2nd Gaussian', 'Amplitude of 2nd Gaussian'),
            Units = ('nm','nm','','nm','nm',''),
            Start = np.asarray([2.5, 0.2, 0.5, 3.5, 0.2, 0.5]),
            Lower = np.asarray([1, 0.05, 0, 1, 0.05, 0]),
            Upper = np.asarray([20, 2.5, 1, 20, 2.5, 1])
        )
        return info
    r,p = _parsargs(args,npar=6)

    r0 = [p[0], p[3]]
    sigma = [p[1], p[4]]
    a = [p[2], p[5]]
    P = _multigaussfun(r,r0,sigma,a)

    return P
#=================================================================
    

def dd_gauss3(*args):
#=================================================================
    r"""
    Sum of three Gaussian distributions
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_gauss3()


    Otherwise the function returns to calculated distance distribution::

        P = dd_gauss3(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -------------------------------------------------------------------------
      Parameter                            Units     Lower    Upper    Start
     -------------------------------------------------------------------------
      Mean of 1st Gaussian                  nm        1        20       2.5 
      Standard deviation of 1st Gaussian    nm       0.05      2.5      0.2 
      Amplitude of 1sd Gaussian                       0        1        0.3 
      Mean of 2nd Gaussian                  nm        1        20       3.5 
      Standard deviation of 2nd Gaussian    nm       0.05      2.5      0.2 
      Amplitude of 2nd Gaussian                       0        1        0.3 
      Mean of 3rd Gaussian                  nm        1        20       5.0 
      Standard deviation of 3rd Gaussian    nm       0.05      2.5      0.2 
      Amplitude of 3rd Gaussian                       0        1        0.3 
     -------------------------------------------------------------------------
    """

    if not args:
        info = dict(
            Parameters = ('Mean of 1st Gaussian', 'Standard deviation of 1st Gaussian', 'Amplitude of 1st Gaussian',
                          'Mean of 2nd Gaussian', 'Standard deviation of 2nd Gaussian', 'Amplitude of 2nd Gaussian',
                          'Mean of 3rd Gaussian', 'Standard deviation of 3rd Gaussian', 'Amplitude of 3rd Gaussian'),
            Units = ('nm','nm','','nm','nm','','nm','nm',''),
            Start = np.asarray([2.5, 0.2, 0.3, 3.5, 0.2, 0.3, 5, 0.2, 0.3]),
            Lower = np.asarray([1, 0.05, 0, 1, 0.05, 0, 1, 0.05, 0]),
            Upper = np.asarray([20, 2.5, 1, 20, 2.5, 1,  20, 2.5, 1])
        )
        return info
    r,p = _parsargs(args,npar=9)

    r0 = [p[0], p[3], p[6]]
    sigma = [p[1], p[4], p[7]]
    a = [p[2], p[5], p[8]]
    P = _multigaussfun(r,r0,sigma,a)

    return P
#=================================================================

def dd_gengauss(*args):    
#=================================================================
    r"""
    Generalized Gaussian distribution model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_gengauss()

    Otherwise the function returns to calculated distance distribution::

        P = dd_gengauss(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -------------------------------------------
      Parameter   Units   Lower   Upper   Start
     -------------------------------------------
      Mean         nm       1      20      3.5 
      Spread       nm      0.05   2.5      0.2 
      Kurtosis             0.25    15       5 
     --------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Mean','Spread','Kurtosis'),
            Units = ('nm','nm',''),
            Start = np.asarray([3.5, 0.5,0.5]),
            Lower = np.asarray([1, 0.05, 0.25]),
            Upper = np.asarray([20, 5, 15])
        )
        return info
    r,p = _parsargs(args,npar=3)
    
    # Compute the model distance distribution
    r0 = p[0]
    sigma = p[1]
    beta = p[2]
    x = abs(r-r0)/sigma
    P = beta/(2*sigma*spc.gamma(1/beta))*np.exp(-x**beta)
    P = _normalize(r,P)

    return P
#=================================================================
    

def dd_skewgauss(*args):    
#=================================================================
    r"""
    Skew Gaussian distribution model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_skewgauss()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_skewgauss(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -------------------------------------------
      Parameter   Units   Lower   Upper   Start
     -------------------------------------------
      Center       nm       1      20      3.5 
      Spread       nm      0.05    2.5     0.2 
      Skewness             -25     25      5 
     --------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Center','Spread','Kurtosis'),
            Units = ('nm','nm',''),
            Start = np.asarray([3.5, 0.2, 5]),
            Lower = np.asarray([1, 0.05, -25]),
            Upper = np.asarray([20, 5, 25])
        )
        return info
    r,p = _parsargs(args,npar=3)
    
    # Compute the model distance distribution
    r0 = p[0]
    sigma = p[1]
    alpha = p[2]
    x = (r-r0)/sigma/np.sqrt(2)
    P = 1/np.sqrt(2*np.pi)*np.exp(-x**2)*(1 + spc.erf(alpha*x))
    P = _normalize(r,P)

    return P
#=================================================================


def dd_rice(*args):    
#=================================================================
    r"""
    3D-Rice distribution
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_rice()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_rice(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ------------------------------------------------
      Parameter  Units     Lower    Upper    Start
     ------------------------------------------------
      Location       nm        1        10       3.5 
      Spread         nm       0.1       5        0.7 
     ------------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Location','Spread'),
            Units = ('nm','nm'),
            Start = np.asarray([3.5, 0.7]),
            Lower = np.asarray([1, 0.1]),
            Upper = np.asarray([10, 5])
        )
        return info
    r,p = _parsargs(args,npar=2)

    nu = [p[0]]
    sig = [p[1]]
    a = [1.0]
    P = _multirice3dfun(r,nu,sig,a)
    return P
#=================================================================
    

def dd_rice2(*args):
#=================================================================
    r"""
    Sum of two 3D-Rice distributions
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_rice2()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_rice2(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------
      Parameter                    Units     Lower    Upper    Start
     ---------------------------------------------------------------
      Location of 1st Rician        nm        1        10       2.5 
      Spread of 1st Rician          nm       0.1       5        0.7 
      Amplitude of 1sd Rician                 0        1        0.5 
      Location of 2nd Rician        nm        1        10       4.0 
      Spread of 2nd Rician          nm       0.1       5        0.7 
      Amplitude of 2nd Rician                 0        1        0.5 
     --------------------------------------------------------------
    """

    if not args:
        info = dict(
            Parameters = ('Location of 1st Rician', 'Spread of 1st Rician', 'Amplitude of 1st Rician',
                          'Location of 2nd Rician', 'Spread of 2nd Rician', 'Amplitude of 2nd Rician'),
            Units = ('nm','nm','','nm','nm',''),
            Start = np.asarray([2.5, 0.7, 0.5, 4.0, 0.7, 0.5]),
            Lower = np.asarray([1, 0.1, 0, 1, 0.1, 0]),
            Upper = np.asarray([10, 5, 1, 10, 5, 1])
        )
        return info
    r,p = _parsargs(args,npar=6)

    nu = [p[0], p[3]]
    sig = [p[1], p[4]]
    a = [p[2], p[5]]
    P = _multirice3dfun(r,nu,sig,a)

    return P
#=================================================================
    

def dd_rice3(*args):
#=================================================================
    r"""
    Sum of three 3D-Rice distributions
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_rice3()

    Otherwise the function returns to calculated distance distribution::
    
        P = dd_rice3(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------
      Parameter                    Units     Lower    Upper    Start
     ---------------------------------------------------------------
      Location of 1st Rician        nm        1        10       2.5 
      Spread of 1st Rician         nm       0.1       5        0.7 
      Amplitude of 1sd Rician               0        1        0.3 
      Location of 2nd Rician        nm        1        10       3.5 
      Spread of 2nd Rician         nm       0.1       5        0.7 
      Amplitude of 2nd Rician               0        1        0.3 
      Location of 3rd Rician        nm        1        10       5.0 
      Spread of 3rd Rician         nm       0.1       5        0.7 
      Amplitude of 3rd Rician               0        1        0.3 
     --------------------------------------------------------------
    """

    if not args:
        info = dict(
            Parameters = ('Location of 1st Rician', 'Spread of 1st Rician', 'Amplitude of 1st Rician',
                          'Location of 2nd Rician', 'Spread of 2nd Rician', 'Amplitude of 2nd Rician',
                          'Location of 3rd Rician', 'Spread of 3rd Rician', 'Amplitude of 3rd Rician'),
            Units = ('nm','nm','','nm','nm','','nm','nm',''),
            Start = np.asarray([2.5, 0.7, 0.3, 3.5, 0.7, 0.3, 5, 0.7, 0.3]),
            Lower = np.asarray([1, 0.1, 0, 1, 0.1, 0, 1, 0.1, 0]),
            Upper = np.asarray([10, 5, 1, 10, 5, 1,  10, 5, 1])
        )
        return info
    r,p = _parsargs(args,npar=9)

    nu = [p[0], p[3], p[6]]
    sig = [p[1], p[4], p[7]]
    a = [p[2], p[5], p[8]]
    P = _multirice3dfun(r,nu,sig,a)

    return P
#=================================================================

def dd_randcoil(*args):    
#=================================================================
    r"""
    Random-coil model for an unfolded peptide/protein
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_randcoil()


    Otherwise the function returns to calculated distance distribution::

        P = dd_randcoil(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ------------------------------------------------------
      Parameter           Units   Lower    Upper    Start
     ------------------------------------------------------
      Number of residues           2       1000      50
      Segment length       nm      0.1     0.4       0.2
      Scaling exponent             0.33    1         0.602
     ------------------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('','nm',''),
            Start = np.asarray([50,   0.2, 0.602]),
            Lower = np.asarray([2,    0.1, 0.33 ]),
            Upper = np.asarray([1000, 0.4, 1    ])
        )
        return info
    r,p = _parsargs(args,npar=3)
    
    N  = p[0]  # number of residues
    nu = p[1] # scaling exponent
    R0 = p[2] # residue length

    rsq = 6*(R0*N**nu)**2 # mean square end-to-end distance from radius of gyration
    normFact = 3/(2*np.pi*rsq)**(3/2) # normalization prefactor
    ShellSurf = 4*np.pi*r**2 # spherical shell surface
    Gaussian = np.exp(-3*r**2/(2*rsq))
    P = normFact*ShellSurf*Gaussian
    P = _normalize(r,P)

    return P
#=================================================================


def dd_circle(*args):    
#=================================================================
    r"""
    Semicircle distribution model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_circle()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_circle(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ----------------------------------------------
      Parameter   Units   Lower    Upper    Start
     ----------------------------------------------
      Center       nm        1      20       3 
      Radius       nm      0.1      5      0.5 
     ----------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Start = np.asarray([3, 0.5]),
            Lower = np.asarray([1, 0.1]),
            Upper = np.asarray([20, 5 ])
        )
        return info
    r,p = _parsargs(args,npar=2)
    
    # Compute the model distance distribution
    r0 = p[0]
    R = abs(p[1])

    dr = r - r0
    idx = abs(dr)<R

    P = np.zeros(len(r))
    P[idx] = 2/np.pi/R**2*np.sqrt(R**2 - dr[idx]**2)
    P = _normalize(r,P)

    return P
#=================================================================


def dd_cos(*args):    
#=================================================================
    r"""
    Raised-cosine parametric model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_cos()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_cos(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ----------------------------------------------
      Parameter   Units   Lower    Upper    Start
     ----------------------------------------------
      Center       nm        1      20       3 
      FWHM         nm      0.1      5      0.5 
     ----------------------------------------------
    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Start = np.asarray([3, 0.5]),
            Lower = np.asarray([1, 0.1]),
            Upper = np.asarray([20, 5 ])
        )
        return info
    r,p = _parsargs(args,npar=2)
    
    # Compute the model distance distribution
    r0 = p[0]
    fwhm = p[1]

    phi = (r-r0)/fwhm*np.pi
    P = (1 + np.cos(phi))/2/fwhm
    P[(r<(r0-fwhm)) | (r>(r0+fwhm))] = 0
    P = _normalize(r,P)

    return P
#=================================================================


def _pb(r,R):
#=================================================================
    P = np.zeros(len(r))
    idx = (r >= 0) & (r <= 2*R)
    P[idx] = 3*r[idx]**5/(16*R**6) - 9*r[idx]**3/(4*R**4) + 3*r[idx]**2/(R**3)

    return P
#=================================================================


def _pbs(r,R1,R2):
#=================================================================
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
#=================================================================



def dd_shell(*args):    
#=================================================================
    r"""
    Uniform spherical shell
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_shell()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_shell(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -------------------------------------------------
      Parameter       Units   Lower    Upper    Start
     -------------------------------------------------
      Shell radius     nm     0.1      20       1.5 
      Shell thickness  nm     0.1      20       0.5 
     -------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1]),
            Upper = np.asarray([20,  20 ]),
            Start = np.asarray([1.5, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=2)
        
    # Compute the model distance distribution
    R1 = float(p[0])
    w = float(p[1])
    R2 = R1 + w

    P = np.zeros(len(r))
    P = R2**6*_pb(r,R2) - R1**6*_pb(r,R1) - 2*(R2**3 - R1**3)*_pbs(r,R1,R2)

    P = P/(R2**3 - R1**3)**2

    P = _normalize(r,P)

    return P
#=================================================================


def dd_spherepoint(*args):    
#=================================================================
    r"""
    One particle distanced from particles distributed on a sphere
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_spherepoint()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_spherepoint(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ---------------------------------------------------
      Parameter         Units   Lower    Upper    Start
     ---------------------------------------------------
      Sphere radius      nm     0.1     20      1.5 
      Distance to point  nm     0.1     20      3.5 
     ---------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1]),
            Upper = np.asarray([20,  20 ]),
            Start = np.asarray([1.5, 3.5]),
        )
        return info
    r,p = _parsargs(args,npar=2)
        
    # Compute the model distance distribution
    R = float(p[0])
    d = float(p[1])
    P = np.zeros(len(r))
    idx = (r >= d - R) & (r<= d + R)
    P[idx] = 3*r[idx]*(R**2 - (d - r[idx])**2)/(4*d*R**3)

    P = _normalize(r,P)

    return P
#=================================================================


def dd_spheresurf(*args):    
#=================================================================
    r"""
    Particles distributed on a sphere's surface
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_spheresurf()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_spheresurf(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ---------------------------------------------------
      Parameter         Units   Lower    Upper    Start
     ---------------------------------------------------
      Sphere radius      nm     0.1     20      2.5 
     ---------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1]),
            Upper = np.asarray([20]),
            Start = np.asarray([2.5]),
        )
        return info
    r,p = _parsargs(args,npar=1)
        
    # Compute the model distance distribution
    R = float(p[0])
    P = np.zeros(len(r))
    idx = (r >= 0) & (r<= 2*R)
    P[idx] = r[idx]/R**2

    P = _normalize(r,P)

    return P
#=================================================================


def dd_shellshell(*args):    
#=================================================================
    r"""
    Uniform spherical shell inside another spherical shell
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_shellshell()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_shellshell(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -------------------------------------------------------
      Parameter             Units   Lower    Upper    Start
     -------------------------------------------------------
      Inner shell radius     nm     0.1      20       1.5 
      Inner shell thickness  nm     0.1      20       0.5 
      Outer shell thickness  nm     0.1      20       0.5 
     -------------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1, 0.1]),
            Upper = np.asarray([20,  20,  20 ]),
            Start = np.asarray([1.5, 0.5, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=3)
        
    # Compute the model distance distribution
    R1 = float(p[0])
    w1 = float(p[1])
    w2 = float(p[2])

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



def dd_shellsphere(*args):    
#=================================================================
    r"""
    Particles distributed on a sphere inside a spherical shell
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_shellsphere()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_shellsphere(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     --------------------------------------------------
      Parameter       Units   Lower    Upper    Start
     --------------------------------------------------
      Sphere radius    nm     0.1      20       1.5 
      Shell thickness  nm     0.1      20       0.5 
     --------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1]),
            Upper = np.asarray([20,  20]),
            Start = np.asarray([1.5, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=2)
        
    # Compute the model distance distribution
    R1 = float(p[0])
    w = float(p[1])
    R2 = R1 + w
    P = _pbs(r,R1,R2)

    P = _normalize(r,P)

    return P
#=================================================================

def dd_shellvoidshell(*args):    
#=================================================================
    r"""
    Particles distributed on a spherical shell inside another spherical shell separated by a void 
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_shellvoidshell()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_shellvoidshell(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     -----------------------------------------------------
      Parameter              Units   Lower   Upper  Start
     -----------------------------------------------------
      Inner sphere radius     nm      0.1     20     0.75 
      Inner shell thickness   nm      0.1     20      1 
      Outer shell thickness   nm      0.1     20      1 
      Shell-shell separation  nm      0.1      2      0.5 
     -----------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1, 0.1, 0.1]),
            Upper = np.asarray([20,  20, 20, 2]),
            Start = np.asarray([0.75, 1, 1, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=4)
        
    # Compute the model distance distribution
    R1 = float(p[0])
    w1 = float(p[1])
    w2 = float(p[2])
    d  = float(p[3])

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
#=================================================================


def dd_shellvoidsphere(*args):    
#=================================================================
    r"""
    Particles distributed on a sphere inside a spherical shell separated by a void 
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_shellvoidsphere()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_shellvoidsphere(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     --------------------------------------------------------
      Parameter               Units   Lower   Upper   Start
     --------------------------------------------------------
      Sphere radius            nm     0.1     20      1.5 
      Shell thickness          nm     0.1     20       1 
      Shell-sphere separation  nm     0.1     20      0.5 
     ---------------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1, 0.1]),
            Upper = np.asarray([20, 20, 20]),
            Start = np.asarray([1.5, 1, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=3)
        
    # Compute the model distance distribution
    R1 = float(p[0])
    w  = float(p[1])
    d  = float(p[2])

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
#=================================================================


def dd_sphere(*args):    
#=================================================================
    r"""
    Particles distributed on a sphere
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_sphere()


    Otherwise the function returns to calculated distance distribution::

        P = dd_sphere(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ------------------------------------------------------
      Parameter              Units   Lower   Upper   Start
     ------------------------------------------------------
      Sphere radius            nm     0.1     20      1.5 
      Shell thickness          nm     0.1     20       1 
      Shell-sphere separation  nm     0.1     20      0.5 
     ------------------------------------------------------

    See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
    http://doi.org/10.1016/j.jmr.2013.01.007

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.1, 0.1]),
            Upper = np.asarray([20,  20, 20]),
            Start = np.asarray([1.5, 1, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=3)
        
    # Compute the model distance distribution
    R = float(p[0])
    P = _pb(r,R)

    P = _normalize(r,P)

    return P
#=================================================================


def dd_triangle(*args):    
#=================================================================
    r"""
    Triangle distribution model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_triangle()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_triangle(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ---------------------------------------------
      Parameter   Units     Lower   Upper   Start
     ---------------------------------------------
      Center       nm       1       20       3.5 
      Width left   nm      0.1       5       0.3 
      Width right  nm      0.1       5       0.3 
     ---------------------------------------------

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([1, 0.1, 0.1]),
            Upper = np.asarray([20,  20, 20]),
            Start = np.asarray([3.5, 1, 0.5]),
        )
        return info
    r,p = _parsargs(args,npar=3)
        
    # Compute the model distance distribution
    r0 = p[0]
    wL = abs(p[1])
    wR = abs(p[2])
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
#=================================================================


def dd_uniform(*args):    
#=================================================================
    r"""
    Uniform distribution model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_uniform()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_uniform(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     --------------------------------------------
      Parameter   Units   Lower   Upper   Start
     --------------------------------------------
      Left edge    nm      0.1     6       2.5 
      Right edge   nm      0.2     20       3 
     --------------------------------------------

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([0.1, 0.2]),
            Upper = np.asarray([6, 20]),
            Start = np.asarray([2.5, 3]),
        )
        return info
    r,p = _parsargs(args,npar=2)
        
    # Compute the model distance distribution
    rL = min(abs(p))
    rR = max(abs(p))
    P = np.zeros(len(r))
    P[(r>=rL) & (r<=rR)] = 1

    P = _normalize(r,P)

    return P
#=================================================================


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

def dd_wormchain(*args):    
#=================================================================
    r"""
    Worm-like chain model near the rigid limit
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_wormchain()


    Otherwise the function returns to calculated distance distribution::
    
        P = dd_wormchain(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ---------------------------------------------------
      Parameter           Units   Lower   Upper   Start
     ---------------------------------------------------
      Chain length         nm      1.5     20      3.7 
      Persistence length   nm       2      100     10 
     ---------------------------------------------------

    See: J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
    https://doi.org/10.1103/PhysRevLett.77.2581

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([1.5, 2]),
            Upper = np.asarray([20, 100]),
            Start = np.asarray([3.7, 10]),
        )
        return info
    r,p = _parsargs(args,npar=2)
    
    # Compute the model distance distribution
    L = p[0]
    Lp = p[1]
    P = wlc(r,L,Lp)

    P = _normalize(r,P)

    return P
#=================================================================


def dd_wormgauss(*args):    
#=================================================================
    r"""
    Worm-like chain model near the rigid limit with Gaussian convolution
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = dd_wormgauss()


    Otherwise the function returns to calculated distance distribution::
    

        P = dd_wormgauss(r,param)
       
 
    Parameters
    ----------
    r : array_like
        TDistance axis, in nanoseconds.
    param : array_like
        List of model parameter values.

    Returns
    -------
    info : dict
        Dictionary containing the built-in information of the model:
        
        * ``info['Parameters']`` - string list of parameter names
        * ``info['Units']`` - string list of metric units of parameters
        * ``info['Start']`` - list of values used as start values during optimization 
        * ``info['Lower']`` - list of values used as lower bounds during optimization 
        * ``info['Upper']`` - list of values used as upper bounds during optimization 
    P : ndarray
        Distance distribution.
 
    Model parameters:
    -------------------

     ------------------------------------------------------------=
      Parameter                    Units   Lower   Upper   Start
     -------------------------------------------------------------
      Chain length                  nm      1.5     20      3.7 
      Persistence length            nm       2      100     10 
      Gaussian standard deviation   nm    0.001      2      0.2 
     -------------------------------------------------------------

    See: J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
    https://doi.org/10.1103/PhysRevLett.77.2581

    """  
    if not args:
        info = dict(
            Parameters = ('Number of residues','Segment length','Scaling exponent'),
            Units = ('nm','nm'),
            Lower = np.asarray([1.5, 2, 0.001]),
            Upper = np.asarray([20, 100, 2]),
            Start = np.asarray([3.7, 10, 0.2]),
        )
        return info
    r,p = _parsargs(args,npar=3)
    
    # Compute the model distance distribution
    L = p[0]
    Lp = p[1]
    sigma = p[2]
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
    print(P.shape)
    P = _normalize(r,P)

    return P
#=================================================================
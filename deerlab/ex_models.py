# ex_models.py - Experiment parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np


def _parsargs(param,npar):
#=================================================================
    param = np.atleast_1d(param)
    if len(param)!=npar:
        raise ValueError('This model requires ',npar,' parameters. A total of ',len(param),' where specified.')
    return param
#=================================================================


def ex_4pdeer(param=None):
# ===================================================================
    r"""
    Single-pathway 4-pulse DEER experiment model 
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_4pdeer()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_4pdeer(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output
        
    pathways : ndarray
        Dipolar pathways of the experiment


    Model parameters:
    -------------------
    
     -------------------------------------------------------
     Parameter            Units   Lower    Upper    Start
     -------------------------------------------------------
     Modulation depth               0        1       0.3 
     -------------------------------------------------------
    """  
    def model(param):     
        # Dipolar pathways
        lam = param[0]
        pathways = [
            [1-lam],
            [lam, 0]
        ]
        return pathways
        
    if param is None:
        info = dict(
            Parameters = ['Modulation depth'],
            Units = [''],
            Start = np.asarray([0.3]),
            Lower = np.asarray([0]),
            Upper = np.asarray([1]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=1) 
        return model(param) 
# ===================================================================


def ex_ovl4pdeer(param=None):
# ===================================================================
    r"""
    4-pulse DEER with band overlap experiment model 
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_ovl4pdeer()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_ovl4pdeer(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output
        
    pathways : ndarray
        Dipolar pathways of the experiment

    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.7
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.1
      Refocusing time of 2nd modulated pathway   μs      0       20        5        
     ---------------------------------------------------------------------------
    """  
    def model(param):   
        # Dipolar pathways
        lam = param[[0,1,2]]
        T0 = param[3]
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0]
        pathways[2] = [lam[2], T0]  
        return pathways  

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway',
                          'Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
            Units = ['','','','μs'],
            Start = np.asarray([0.7, 0.3, 0.1, 5]),
            Lower = np.asarray([0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 20]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=4) 
        return model(param) 
# ===================================================================



def ex_5pdeer(param=None):
# ===================================================================
    r"""
    5-pulse DEER experiment model
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_5pdeer()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_5pdeer(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output
        
    pathways : ndarray
        Dipolar pathways of the experiment
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.4
      Amplitude of 1st modulated pathway                 0       1        0.4
      Amplitude of 2nd modulated pathway                 0       1        0.2
      Refocusing time of 2nd modulated pathway   μs      0       20        5        
     ---------------------------------------------------------------------------
    """  
    def model(param):   
                  # Dipolar pathways
        lam = param[[0,1,2]]
        T0 = param[3]
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0]
        pathways[2] = [lam[2], T0]
        return pathways  

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway','Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
            Units = ['','','','μs'],
            Start = np.asarray([0.4, 0.4, 0.2,5]),
            Lower = np.asarray([0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 20]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=4) 
        return model(param) 
# ===================================================================



def ex_7pdeer(param=None):
# ===================================================================
    r"""
    7-pulse DEER experiment model
    ==============================
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_7pdeer()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_7pdeer(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output
        
    pathways : ndarray
        Dipolar pathways of the experiment
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.3
      Amplitude of 1st modulated pathway                 0       1        0.5
      Amplitude of 2nd modulated pathway                 0       1        0.3
      Amplitude of 3rd modulated pathway                 0       1        0.2 
      Refocusing time of 2nd modulated pathway   μs      0       20       1.5        
      Refocusing time of 3rd modulated pathway   μs      0       20       3.5 
     ---------------------------------------------------------------------------
    """  
    def model(param):   
        # Dipolar pathways
        lam = param[[0,1,2,3]]
        T0 = param[[4,5]]
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0]
        pathways[2] = [lam[2], T0[0]]
        pathways[3] = [lam[3], T0[1]]    
        return pathways  

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway',
                          'Amplitude of 2nd modulated pathway','Amplitude of 3rd modulated pathway',
                          'Refocusing time of 2nd modulated pathway','Refocusing time of 3rd modulated pathway'],
            Units = ['','','','','μs','μs'],
            Start = np.asarray([0.3, 0.5, 0.3, 0.2, 1.5, 3.5]),
            Lower = np.asarray([0, 0, 0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 1, 20, 20]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=6) 
        return model(param) 

# ===================================================================


def ex_ridme1(param=None):
# ===================================================================
    r"""
    RIDME experiment model (spin S=1/2)
    ===================================
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_ridme1()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_ridme1(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output
        
    pathways : ndarray
        Dipolar pathways of the experiment
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated contribution              0       1        0.3
      Amplitude of 1st harmonic pathway                  0       1        0.5
     ---------------------------------------------------------------------------
    """  
    def model(param):   
        # Dipolar pathways
        lam = param.copy()
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0, 1]
        return pathways      

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st modulated pathway'],
            Units = ['',''],
            Start = np.asarray([0.3, 0.5]),
            Lower = np.asarray([0, 0]),
            Upper = np.asarray([1, 1]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=2) 
        return model(param) 

# ===================================================================


def ex_ridme3(param=None):
# ===================================================================
    r"""
    RIDME experiment model (spin S=3/2)
    ====================================
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_ridme3()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_ridme3(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output

    pathways : ndarray
        Dipolar pathways of the experiment
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated contribution              0       1        0.3
      Amplitude of 1st harmonic pathway                  0       1        0.5
      Amplitude of 2nd harmonic pathway                  0       1        0.3
      Amplitude of 3rd harmonic pathway                  0       1        0.2
     ---------------------------------------------------------------------------
    """  
    def model(param):   
        # Dipolar pathways
        lam = param.copy()
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0, 1]
        pathways[2] = [lam[2], 0, 2]
        pathways[3] = [lam[3], 0, 3]
        return pathways

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st harmonic pathway',
                          'Amplitude of 2nd harmonic pathway','Amplitude of 3rd harmonic pathway'],
            Units = ['','','',''],
            Start = np.asarray([0.3, 0.5, 0.3, 0.2]),
            Lower = np.asarray([0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 1]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=4) 
        return model(param) 

# ===================================================================


def ex_ridme5(param=None):
# ===================================================================
    r"""
    RIDME experiment model (spin S=5/2)
    ====================================
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_ridme5()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_ridme5(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output

    pathways : ndarray
        Dipolar pathways of the experiment
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated component                 0       1        0.3
      Amplitude of 1st harmonic pathway                  0       1        0.5
      Amplitude of 2nd harmonic pathway                  0       1        0.3
      Amplitude of 3rd harmonic pathway                  0       1        0.2 
      Amplitude of 4th harmonic pathway                  0       1        0.1 
      Amplitude of 5th harmonic pathway                  0       1        0.05 
     ---------------------------------------------------------------------------
    """  
    def model(param):     
        # Dipolar pathways
        lam = param.copy()
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0, 1]
        pathways[2] = [lam[2], 0, 2]
        pathways[3] = [lam[3], 0, 3]
        pathways[4] = [lam[4], 0, 4]
        pathways[5] = [lam[5], 0, 5]
        return pathways

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st harmonic pathway',
                          'Amplitude of 2nd harmonic pathway','Amplitude of 3rd harmonic pathway',
                          'Amplitude of 4th harmonic pathway','Amplitude of 5th harmonic pathway'],
            Units = ['','','','','',''],
            Start = np.asarray([0.3, 0.5, 0.3, 0.2,0.1,0.05]),
            Lower = np.asarray([0, 0, 0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 1, 1, 1]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=6) 
        return model(param) 

# ===================================================================


def ex_ridme7(param=None):
# ===================================================================
    r"""
    RIDME experiment model (spin S=7/2)
    ====================================
 
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = ex_ridme7()


    Otherwise the function returns to calculated experiment dipolar pathways::
    
        pathways = ex_ridme7(param)
 
 
    Parameters
    ----------
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
        * ``info['ModelFcn']`` - function used to calculate the model output 
    pathways : ndarray
        Dipolar pathways of the experiment
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated component                 0       1        0.3
      Amplitude of 1st harmonic pathway                  0       1        0.5
      Amplitude of 2nd harmonic pathway                  0       1        0.3
      Amplitude of 3rd harmonic pathway                  0       1        0.2 
      Amplitude of 4th harmonic pathway                  0       1        0.1 
      Amplitude of 5th harmonic pathway                  0       1        0.05 
      Amplitude of 6th harmonic pathway                  0       1        0.02 
      Amplitude of 7th harmonic pathway                  0       1        0.01 
     ---------------------------------------------------------------------------
    """  
    def model(param):        
        # Dipolar pathways
        lam = param.copy()
        pathways = [[] for _ in lam]
        pathways[0] = [lam[0]]
        pathways[1] = [lam[1], 0, 1]
        pathways[2] = [lam[2], 0, 2]
        pathways[3] = [lam[3], 0, 3]
        pathways[4] = [lam[4], 0, 4]
        pathways[5] = [lam[5], 0, 5]
        pathways[6] = [lam[6], 0, 6]
        pathways[7] = [lam[7], 0, 7]
        return pathways

    if param is None:
        info = dict(
            Parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st harmonic pathway',
                          'Amplitude of 2nd harmonic pathway','Amplitude of 3rd harmonic pathway',
                          'Amplitude of 4th harmonic pathway','Amplitude of 5th harmonic pathway',
                          'Amplitude of 6th harmonic pathway','Amplitude of 7th harmonic pathway'],
            Units = ['','','','','','','',''],
            Start = np.asarray([0.3, 0.5, 0.3, 0.2,0.1,0.05,0.02,0.01]),
            Lower = np.asarray([0, 0, 0, 0, 0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 1, 1, 1, 1, 1]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=8) 
        return model(param) 

# ===================================================================
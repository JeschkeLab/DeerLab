# ex_models.py - Experiment parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np

# Definition of the header for all experiment models
docstr_header = lambda fcnstr: """
    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries::

        info = {}()


    Otherwise the function returns to calculated experiment dipolar pathways::

        pathways = {}(param)


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
""".format(fcnstr,fcnstr)

# =================================================================
def docstring():
    """
    Decorator: Insert docstring header to a pre-existing docstring
    """
    sep="\n"
    def _decorator(func):
        docstr = func.__doc__
        docstr = docstr.split("Model parameters:",1)[0]
        func.__doc__ = sep.join([docstr,docstr_header(func.__name__)])
        return func
    return _decorator
# =================================================================

def _parsargs(param,npar):
# =================================================================
    r"""Simple validation of input parameter vector."""
    param = np.atleast_1d(param)
    if len(param)!=npar:
        raise ValueError('This model requires ',npar,' parameters. A total of ',len(param),' where specified.')
    return param
# =================================================================

def _simplified_pathway_model(param,par0):
# =================================================================
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
            Start = np.asarray(par0),
            Lower = np.asarray([0]),
            Upper = np.asarray([1]),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=1) 
        return model(param) 
# =================================================================


def _dipolarpathways_model(param,par0,Npathways=1):
# =================================================================
    def model(param):   
        lam_idx = np.arange(2,2+2*(Npathways-1),2)
        T0_idx = np.arange(3,3+2*(Npathways-1),2)

        pathways = [[] for _ in range(Npathways+1)]
        for n in range(Npathways+1):
            if n==0:
                pathways[n] = [param[n]]
            elif n==1:
                pathways[n] = [param[n], 0]
            else:
                pathways[n] = [param[lam_idx[n-2]], param[T0_idx[n-2]] ]
        return pathways  

    if param is None:
        ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
        parameters,lb,ub,units = [],[],[],[]
        for n in range(Npathways+1):
            if n==0:
                parameters.append('Amplitude of unmodulated components')
                units.append('')
                lb.append(0)
                ub.append(1)
            elif n==1:
                parameters.append('Amplitude of 1st modulated pathway')
                units.append('')
                lb.append(0)
                ub.append(1)
            else:
                parameters.append('Amplitude of {} modulated pathway'.format(ordinal(n)))
                units.append('')
                lb.append(0)
                ub.append(1)
                parameters.append('Refocusing time of {} modulated pathway'.format(ordinal(n)))
                units.append('μs')
                lb.append(0)
                ub.append(20)

        info = dict(
            Parameters = parameters,
            Units = units,
            Start = np.asarray(par0),
            Lower = np.asarray(lb),
            Upper = np.asarray(ub),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=2*Npathways) 
        return model(param) 
# =================================================================

def _harmonicpathways_model(param,par0,Npathways=1):
# =================================================================
    def model(param):   
        pathways = [[] for _ in range(Npathways+1)]
        for n in range(Npathways+1):
            if n==0:
                pathways[n] = [param[n]]
            else:
                pathways[n] = [param[n], 0, n]
        return pathways  

    if param is None:
        ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
        parameters,lb,ub,units = [],[],[],[]
        for n in range(Npathways+1):
            if n==0:
                parameters.append('Amplitude of unmodulated components')
                units.append('')
                lb.append(0)
                ub.append(1)
            else:
                parameters.append('Amplitude of {} harmonic pathway'.format(ordinal(n)))
                units.append('')
                lb.append(0)
                ub.append(1)

        info = dict(
            Parameters = parameters,
            Units = units,
            Start = np.asarray(par0),
            Lower = np.asarray(lb),
            Upper = np.asarray(ub),
            ModelFcn = model
        )
        return info
    else:
        param = _parsargs(param, npar=1+Npathways) 
        return model(param) 
# =================================================================

@docstring()
def ex_3pdeer1(param=None):
# ===================================================================
    r"""
    Single dipolar pathways 3-pulse DEER experiment model 

    Model parameters:
    -------------------
    
     -------------------------------------------------------
     Parameter            Units   Lower    Upper    Start
     -------------------------------------------------------
     Modulation depth               0        1       0.3 
     -------------------------------------------------------
    """
    par0 = [0.3]
    return _simplified_pathway_model(param,par0)
# ===================================================================


@docstring()
def ex_3pdeer2(param=None):
# ===================================================================
    r"""
    Two dipolar pathways 3-pulse DEER experiment model 


    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.6
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.1
      Refocusing time of 2nd modulated pathway   μs      0       20        5        
     ---------------------------------------------------------------------------
    """
    par0 = [0.6,0.3,0.1,5]
    return _dipolarpathways_model(param,par0,Npathways=2)
# ===================================================================

@docstring()
def ex_4pdeer1(param=None):
# ===================================================================
    r"""
    Single dipolar pathways 4-pulse DEER experiment model 

    Model parameters:
    -------------------
    
     -------------------------------------------------------
     Parameter            Units   Lower    Upper    Start
     -------------------------------------------------------
     Modulation depth               0        1       0.3 
     -------------------------------------------------------
    """
    par0 = [0.3]
    return _simplified_pathway_model(param,par0)
# ===================================================================

@docstring()
def ex_4pdeer2(param=None):
# ===================================================================
    r"""
    Two dipolar pathways 4-pulse DEER experiment model 


    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.6
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.1
      Refocusing time of 2nd modulated pathway   μs      0       20        5        
     ---------------------------------------------------------------------------
    """
    par0 = [0.6,0.3,0.1,5]
    return _dipolarpathways_model(param,par0,Npathways=2)
# ===================================================================

@docstring()
def ex_4pdeer3(param=None):
# ===================================================================
    r"""
    Three dipolar pathways 4-pulse DEER experiment model 


    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.6
      Amplitude of 1st modulated pathway                 0       1        0.2
      Amplitude of 2nd modulated pathway                 0       1        0.1
      Refocusing time of 2nd modulated pathway   μs      0       20        5 
      Amplitude of 3rd modulated pathway                 0       1        0.1
      Refocusing time of 3rd modulated pathway   μs      0       20       -5          
     ---------------------------------------------------------------------------
    """
    par0 = [0.6, 0.2, 0.1, 5, 0.1, -5]
    return _dipolarpathways_model(param,par0,Npathways=3)
# ===================================================================


@docstring()
def ex_5pdeer1(param=None):
# ===================================================================
    r"""
    Single dipolar pathways 5-pulse DEER experiment model 

    Model parameters:
    -------------------
    
     -------------------------------------------------------
     Parameter            Units   Lower    Upper    Start
     -------------------------------------------------------
     Modulation depth               0        1       0.3 
     -------------------------------------------------------
    """
    par0 = [0.3]
    return _simplified_pathway_model(param,par0)
# ===================================================================


@docstring()
def ex_5pdeer2(param=None):
# ===================================================================
    r"""
    Two dipolar pathways 5-pulse DEER experiment model 
 
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
    par0 = [0.4, 0.4, 0.2, 5]
    return _dipolarpathways_model(param,par0,Npathways=2)
# ===================================================================

@docstring()
def ex_5pdeer3(param=None):
# ===================================================================
    r"""
    Three dipolar pathways 5-pulse DEER experiment model 
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.4
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.2
      Refocusing time of 2nd modulated pathway   μs      0       20        5        
      Amplitude of 3rd modulated pathway                 0       1        0.1
      Refocusing time of 2nd modulated pathway   μs      0       20        2 
     ---------------------------------------------------------------------------
    """  
    par0 = [0.4, 0.3, 0.2, 5, 0.1, 2]
    return _dipolarpathways_model(param,par0,Npathways=3)
# ===================================================================

@docstring()
def ex_5pdeer4(param=None):
# ===================================================================
    r"""
    Four dipolar pathways 5-pulse DEER experiment model 
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.4
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.2
      Refocusing time of 2nd modulated pathway   μs      0       20        5        
      Amplitude of 3rd modulated pathway                 0       1        0.05
      Refocusing time of 2nd modulated pathway   μs      0       20        2        
      Amplitude of 4th modulated pathway                 0       1        0.05
      Refocusing time of 4th modulated pathway   μs      0       20        3 
     ---------------------------------------------------------------------------
    """  
    par0 = [0.4, 0.3, 0.2, 5, 0.05, 2, 0.05, 3]
    return _dipolarpathways_model(param,par0,Npathways=4)
# ===================================================================


@docstring()
def ex_7pdeer1(param=None):
# ===================================================================
    r"""
    Single dipolar pathways 7-pulse DEER experiment model 

    Model parameters:
    -------------------
    
     -------------------------------------------------------
     Parameter            Units   Lower    Upper    Start
     -------------------------------------------------------
     Modulation depth               0        1       0.3 
     -------------------------------------------------------
    """
    par0 = [0.3]
    return _simplified_pathway_model(param,par0)
# ===================================================================


@docstring()
def ex_7pdeer2(param=None):
# ===================================================================
    r"""
    Two dipolar pathways 7-pulse DEER experiment model 


    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.3
      Amplitude of 1st modulated pathway                 0       1        0.5
      Amplitude of 2nd modulated pathway                 0       1        0.3
      Refocusing time of 2nd modulated pathway   μs      0       20       1.5 
     ---------------------------------------------------------------------------
    """  
    par0 = [0.5, 0.3, 0.2, 1.5]
    return _dipolarpathways_model(param,par0,Npathways=2)
# ===================================================================

@docstring()
def ex_7pdeer3(param=None):
# ===================================================================
    r"""
    Three dipolar pathways 7-pulse DEER experiment model 


    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.3
      Amplitude of 1st modulated pathway                 0       1        0.5
      Amplitude of 2nd modulated pathway                 0       1        0.3
      Refocusing time of 2nd modulated pathway   μs      0       20       1.5 
      Amplitude of 3rd modulated pathway                 0       1        0.2        
      Refocusing time of 3rd modulated pathway   μs      0       20       3.5 
     ---------------------------------------------------------------------------
    """  
    par0 = [0.5, 0.3, 0.1, 1.5, 0.1, 3.5]
    return _dipolarpathways_model(param,par0,Npathways=3)
# ===================================================================


@docstring()
def ex_7pdeer4(param=None):
# ===================================================================
    r"""
    Four dipolar pathways 7-pulse DEER experiment model 


    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.4
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.1
      Refocusing time of 2nd modulated pathway   μs      0       20       1.5 
      Amplitude of 3rd modulated pathway                 0       1        0.05        
      Refocusing time of 3rd modulated pathway   μs      0       20       3.5
      Amplitude of 3rd modulated pathway                 0       1        0.05        
      Refocusing time of 3rd modulated pathway   μs      0       20       5.5  
     ---------------------------------------------------------------------------
    """  
    par0 = [0.4, 0.3, 0.1, 1.5, 0.05, 3.5, 5.5, 0.05]
    return _dipolarpathways_model(param,par0,Npathways=4)
# ===================================================================

@docstring()
def ex_ridme1(param=None):
# ===================================================================
    r"""
    Two harmonic pathways RIDME experiment model

 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated contribution              0       1        0.3
      Amplitude of 1st harmonic pathway                  0       1        0.5
     ---------------------------------------------------------------------------
    """  
    par0 = [0.3, 0.5]
    return _harmonicpathways_model(param,par0,Npathways=1)
# ===================================================================

@docstring()
def ex_ridme3(param=None):
# ===================================================================
    r"""
    Three harmonic pathways RIDME experiment model

    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated contribution              0       1        0.3
      Amplitude of 1st harmonic pathway                  0       1        0.4
      Amplitude of 2nd harmonic pathway                  0       1        0.2
      Amplitude of 3rd harmonic pathway                  0       1        0.1
     ---------------------------------------------------------------------------
    """  
    par0 = [0.3,0.4,0.2,0.1]
    return _harmonicpathways_model(param,par0,Npathways=3)
# ===================================================================

@docstring()
def ex_ridme5(param=None):
# ===================================================================
    r"""
    Five harmonic pathways RIDME experiment model
   
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated component                 0       1        0.30
      Amplitude of 1st harmonic pathway                  0       1        0.35
      Amplitude of 2nd harmonic pathway                  0       1        0.15
      Amplitude of 3rd harmonic pathway                  0       1        0.10 
      Amplitude of 4th harmonic pathway                  0       1        0.05 
      Amplitude of 5th harmonic pathway                  0       1        0.05 
     ---------------------------------------------------------------------------
    """  
    par0 = [0.3,0.35,0.15,0.1,0.05,0.05]
    return _harmonicpathways_model(param,par0,Npathways=5)
# ===================================================================

@docstring()
def ex_ridme7(param=None):
# ===================================================================
    r"""

    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated component                 0       1        0.30
      Amplitude of 1st harmonic pathway                  0       1        0.35
      Amplitude of 2nd harmonic pathway                  0       1        0.15
      Amplitude of 3rd harmonic pathway                  0       1        0.10 
      Amplitude of 4th harmonic pathway                  0       1        0.05 
      Amplitude of 5th harmonic pathway                  0       1        0.02 
      Amplitude of 6th harmonic pathway                  0       1        0.02 
      Amplitude of 7th harmonic pathway                  0       1        0.01 
     ---------------------------------------------------------------------------
    """  
    par0 = [0.3,0.35,0.15,0.1,0.05,0.02,0.02,0.01]
    return _harmonicpathways_model(param,par0,Npathways=7)
# ===================================================================
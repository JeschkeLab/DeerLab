

import numpy as np
from deerlab.utils import isempty
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

def _parsargs(param,npar):
#=================================================================
    param = np.atleast_1d(param)
    if len(param)!=npar:
        raise ValueError('This model requires ',npar,' parameters. A total of ',len(param),' where specified.')
    return param
#=================================================================


def ex_4pdeer(param=[]):
# ===================================================================
    """
    Single-pathway 4-pulse DEER experiment model 
    =============================================
   
     info = ex_4pdeer()
     P = ex_4pdeer(param)

    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries. Otherwise, computes the N-point model ``P`` from the N-point distance axis ``r`` according to 
    the paramteres array ``param`. The required parameters can also be found 
    in the ``info`` structure.
 
    Model parameters:
    -------------------

     -------------------------------------------------------
      Parameter            Units   Lower    Upper    Start
     -------------------------------------------------------
      Modulation depth               0        1       0.3 
     -------------------------------------------------------
    """  
    if isempty(param):
        info = dict(
            Parameters = ['Modulation depth'],
            Units = [''],
            Start = np.asarray([0.3]),
            Lower = np.asarray([0]),
            Upper = np.asarray([1])
        )
        return info
    param = _parsargs(param,npar=1)

    # Dipolar pathways
    lam = param[0]
    pathways = np.array([[1-lam, np.NaN], [lam, 0]])
    
    return pathways
# ===================================================================


def ex_ovl4pdeer(param=[]):
# ===================================================================
    """
    4-pulse DEER with band overlap experiment model 
    ================================================
   
     info = ex_ovl4pdeer()
     P = ex_ovl4pdeer(param)

    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries. Otherwise, computes the N-point model ``P`` from the N-point distance axis ``r`` according to 
    the paramteres array ``param`. The required parameters can also be found 
    in the ``info`` structure.
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.3
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.3
      Refocusing time of 2nd modulated pathway   us      0       20        5        
     ---------------------------------------------------------------------------
    """  
    if isempty(param):
        info = dict(
            Parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway','Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
            Units = ['','','','us'],
            Start = np.asarray([0.3, 0.3, 0.3,5]),
            Lower = np.asarray([0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 20])
        )
        return info
    param = _parsargs(param,npar=4)

    # Dipolar pathways
    lam = param[[0,1,2]]
    T02 = param[3]
    pathways = np.array([[lam[0], np.NaN], [lam[1], 0], [lam[2], T02]])
    
    return pathways
# ===================================================================



def ex_5pdeer(param=[]):
# ===================================================================
    """
    5-pulse DEER experiment model
    ==============================
   
     info = ex_5pdeer()
     P = ex_5pdeer(param)

    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries. Otherwise, computes the N-point model ``P`` from the N-point distance axis ``r`` according to 
    the paramteres array ``param`. The required parameters can also be found 
    in the ``info`` structure.
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.3
      Amplitude of 1st modulated pathway                 0       1        0.3
      Amplitude of 2nd modulated pathway                 0       1        0.3
      Refocusing time of 2nd modulated pathway   us      0       20        5        
     ---------------------------------------------------------------------------
    """  
    if isempty(param):
        info = dict(
            Parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway','Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
            Units = ['','','','us'],
            Start = np.asarray([0.3, 0.3, 0.3,5]),
            Lower = np.asarray([0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 20])
        )
        return info
    param = _parsargs(param,npar=4)

    # Dipolar pathways
    lam = param[[0,1,2]]
    T02 = param[3]
    pathways = np.array([[lam[0], np.NaN], [lam[1], 0], [lam[2], T02]])
    
    return pathways
# ===================================================================



def ex_7pdeer(param=[]):
# ===================================================================
    """
    7-pulse DEER experiment model
    ==============================
   
     info = ex_7pdeer()
     P = ex_7pdeer(param)

    If called without arguments, returns an ``info`` dictionary of model parameters and boundaries. Otherwise, computes the N-point model ``P`` from the N-point distance axis ``r`` according to 
    the paramteres array ``param`. The required parameters can also be found 
    in the ``info`` structure.
 
    Model parameters:
    -------------------

     ---------------------------------------------------------------------------
      Parameter                                Units   Lower    Upper    Start
     ---------------------------------------------------------------------------
      Amplitude of unmodulated components                0       1        0.3
      Amplitude of 1st modulated pathway                 0       1        0.5
      Amplitude of 2nd modulated pathway                 0       1        0.3
      Amplitude of 3rd modulated pathway                 0       1        0.2 
      Refocusing time of 2nd modulated pathway   us      0       20       1.5        
      Refocusing time of 3rd modulated pathway   us      0       20       3.5 
     ---------------------------------------------------------------------------
    """  
    if isempty(param):
        info = dict(
            Parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway','Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
            Units = ['','','','us'],
            Start = np.asarray([0.3, 0.5, 0.3, 0.2, 1.5, 3.5]),
            Lower = np.asarray([0, 0, 0, 0, 0, 0]),
            Upper = np.asarray([1, 1, 1, 1, 20, 20])
        )
        return info
    param = _parsargs(param,npar=6)

    # Dipolar pathways
    lam = param[[0,1,2,3]]
    T02 = param[4]
    T03 = param[5]

    pathways = np.array([[lam[0], np.NaN], [lam[1], 0], [lam[2], T02], [lam[3], T03]])
    
    return pathways
# ===================================================================
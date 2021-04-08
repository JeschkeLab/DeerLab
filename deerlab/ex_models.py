# ex_models.py - Experiment parametric models
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.utils import metadata 

# =================================================================
def docstr_header(title,fcnstr):
    "Definition of the header for all experiment models"
    return f"""
{title}

The function takes a list or array of parameters and returns the calculated experiment dipolar pathways::

        pathways = {fcnstr}(param)

The built-in information on the model can be accessed via its attributes::

        {fcnstr}.parameters  # String list of parameter names
        {fcnstr}.units       # String list of metric units of parameters
        {fcnstr}.start       # List of values used as start values during optimization 
        {fcnstr}.lower       # List of values used as lower bounds during optimization
        {fcnstr}.upper       # List of values used as upper bounds during optimization 

Parameters
----------
param : array_like
    List of model parameter values.

Returns
-------
    
pathways : ndarray
    Dipolar pathways of the experiment
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
        return func
    return _decorator
# =================================================================

#=================================================================
def _parsargs(param,npar):
    param = np.atleast_1d(param)
    if len(param)!=npar:
        raise ValueError('This model requires ',npar,' parameters. A total of ',len(param),' where specified.')
    return param
#=================================================================


# =================================================================
@metadata(
parameters = ['Modulation depth'],
units = [''],
start = np.asarray([0.3]),
lower = np.asarray([0]),
upper = np.asarray([1]))
@docstring()
def ex_4pdeer(param):
    r"""
Single-pathway 4-pulse DEER experiment model 
        
Notes
-----

**Model:**

This experiment model has one modulated pathway and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [1-\lambda + \lambda K_0(t-T_0^{(1)},r)]B(t-T_0^{(1)},\lambda)

where :math:`T_0^{(1)}=0` is the refocusing time of the modulated dipolar pathway.


============== ================ ============= ============ ============ ================================================
 Variable        Symbol         Start Values     Lower        Upper                Description
============== ================ ============= ============ ============ ================================================
``param[0]``   :math:`\lambda`     0.3           0            1          Modulated pathway amplitude (modulation depth)
============== ================ ============= ============ ============ ================================================
    """  
    param = _parsargs(param,npar=1) 
    
    # Dipolar pathways
    lam = param[0]
    pathways = [
        [1-lam],
        [lam, 0]
    ]
    return pathways
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway',
                'Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
units = ['','','','μs'],
start = np.asarray([0.7, 0.3, 0.1, 5]),
lower = np.asarray([0, 0, 0, 0]),
upper = np.asarray([1, 1, 1, 20]))
@docstring()
def ex_ovl4pdeer(param):
    r"""
4-pulse DEER with band overlap experiment model 
        
Notes
-----

**Model:**

This experiment model has two modulated pathways and an unmodulated contribution. The second modulated pathway results in the 2+1 contribution at the end of a 4-pulse DEER trace.The kernel is 

.. math::
   K(t,r) =
   [\Lambda_0 + \lambda_1 K_0(t-T_0^{(1)},r) + \lambda_2 K_0(t-T_0^{(2)},r)]
   B(t-T_0^{(1)},\lambda_1)
   B(t-T_0^{(2)},\lambda_2)

where :math:`T_0^{(1)}=0` and :math:`T_0^{(2)}` are the refocusing times of the two modulated dipolar pathways.


============== ======================== ============= ============ ============ ================================================
 Variable        Symbol                  Start Values     Lower        Upper                Description
============== ======================== ============= ============ ============ ================================================
``param[0]``   :math:`\varLambda_0`        0.7              0            1        Unmodulated pathways, amplitude
``param[1]``   :math:`\lambda_1`           0.3              0            1        1st modulated pathway, amplitude
``param[2]``   :math:`\lambda_2`           0.1              0            1        2nd modulated pathway, amplitude
``param[3]``   :math:`T_0^{(2)}`           5.0              0           20        2nd modulated pathway, refocusing time (μs)
============== ======================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param,npar=4) 

    # Dipolar pathways
    lam = param[[0,1,2]]
    T0 = param[3]
    pathways = [[] for _ in lam]
    pathways[0] = [lam[0]]
    pathways[1] = [lam[1], 0]
    pathways[2] = [lam[2], T0]  
    return pathways 
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway','Amplitude of 2nd modulated pathway','Refocusing time of 2nd modulated pathway'],
units = ['','','','μs'],
start = np.asarray([0.4, 0.4, 0.2,5]),
lower = np.asarray([0, 0, 0, 0]),
upper = np.asarray([1, 1, 1, 20]))
@docstring()
def ex_5pdeer(param):
    r"""
5-pulse DEER experiment model
        
Notes
-----

**Model:**

This experiment model has two modulated pathway and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [\varLambda_0 + \lambda_1 K_0(t-T_0^{(1)},r) + \lambda_2 K_0(t-T_0^{(2)},r)]
   B(t-T_0^{(1)},\lambda_1) B(t - T_0^{(2)},\lambda_2)

where :math:`T_0^{(1)}=0` and :math:`T_0^{(2)}` are the refocusing times of the two modulated dipolar pathways.


============== ======================== ============= ============ ============ ================================================
 Variable        Symbol                  Start Values     Lower        Upper                Description
============== ======================== ============= ============ ============ ================================================
``param[0]``   :math:`\varLambda_0`          0.4            0            1       Unmodulated pathways, amplitude
``param[1]``   :math:`\lambda_1`             0.4            0            1       1st modulated pathway, amplitude
``param[2]``   :math:`\lambda_2`             0.2            0            1       2nd modulated pathway, amplitude
``param[3]``   :math:`T_0^{(2)}`             5.0            0            20      2nd modulated pathway, refocusing time (μs)
============== ======================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param,npar=4) 

    # Dipolar pathways
    lam = param[[0,1,2]]
    T0 = param[3]
    pathways = [[] for _ in lam]
    pathways[0] = [lam[0]]
    pathways[1] = [lam[1], 0]
    pathways[2] = [lam[2], T0]
    return pathways
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated components','Amplitude of 1st modulated pathway',
                'Amplitude of 2nd modulated pathway','Amplitude of 3rd modulated pathway',
                'Refocusing time of 2nd modulated pathway','Refocusing time of 3rd modulated pathway'],
units = ['','','','','μs','μs'],
start = np.asarray([0.3, 0.5, 0.3, 0.2, 1.5, 3.5]),
lower = np.asarray([0, 0, 0, 0, 0, 0]),
upper = np.asarray([1, 1, 1, 1, 20, 20]))
@docstring()
def ex_7pdeer(param):
    r"""
7-pulse DEER experiment model
   
Notes
-----

**Model:**

In order to reduce the parameter space, only the dipolar pathways refocusing at positive times (pathways #1-3) are considered in this model:

.. math::

    K(t,r) =
    [\Lambda_0 + \sum^3_{p=1} \lambda_p K_0(t-T_0^{(p)},r)]
    \prod^3_{p=1} B(t - T_0^{(p)},\lambda_p)

where :math:`T_0^{(1)}=0\;\mu s`, :math:`T_0^{(2)}`, and :math:`T_0^{(3)}` are the refocusing times of the three modulated dipolar pathways at positive evolution times.


============== ======================== ============= ============ ============ ================================================
 Variable        Symbol                  Start Values     Lower        Upper                Description
============== ======================== ============= ============ ============ ================================================
``param[0]``   :math:`\varLambda_0`        0.3                0       1          Unmodulated pathways, amplitude
``param[1]``   :math:`\lambda_1`           0.5                0       1          1st modulated pathway, amplitude
``param[2]``   :math:`\lambda_2`           0.3                0       1          2nd modulated pathway, amplitude
``param[3]``   :math:`\lambda_3`           0.2                0       1          3rd modulated pathway, amplitude
``param[4]``   :math:`T_0^{(2)}`           1.5                0       20         2nd modulated pathway, refocusing time (μs)
``param[5]``   :math:`T_0^{(3)}`           3.5                0       20         3rd modulated pathway, refocusing time (μs)
============== ======================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param,npar=6) 

    # Dipolar pathways
    lam = param[[0,1,2,3]]
    T0 = param[[4,5]]
    pathways = [[] for _ in lam]
    pathways[0] = [lam[0]]
    pathways[1] = [lam[1], 0]
    pathways[2] = [lam[2], T0[0]]
    pathways[3] = [lam[3], T0[1]]    
    return pathways  
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st modulated pathway'],
units = ['',''],
start = np.asarray([0.3, 0.5]),
lower = np.asarray([0, 0]),
upper = np.asarray([1, 1]))
@docstring()
def ex_ridme1(param):
    r"""
RIDME experiment model (spin S=1/2)
        
Notes
-----

**Model:**

This experiment model has one harmonic pathway and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [\Lambda_0 + \lambda_1 K_0(t,r)]B(t)

============== =================== ============= ============ ============ ================================================
 Variable        Symbol            Start Values     Lower        Upper                Description
============== =================== ============= ============ ============ ================================================
``param[0]``   :math:`\Lambda_0`     0.5           0            1          Amplitude of unmodulated contribution
``param[1]``   :math:`\lambda_1`     0.3           0            1          Amplitude of 1st harmonic pathway
============== =================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param, npar=2) 
    
    # Dipolar pathways
    lam = param.copy()
    pathways = [[] for _ in lam]
    pathways[0] = [lam[0]]
    pathways[1] = [lam[1], 0, 1]
    return pathways   
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st harmonic pathway',
                'Amplitude of 2nd harmonic pathway','Amplitude of 3rd harmonic pathway'],
units = ['','','',''],
start = np.asarray([0.3, 0.5, 0.3, 0.2]),
lower = np.asarray([0, 0, 0, 0]),
upper = np.asarray([1, 1, 1, 1]))
@docstring()
def ex_ridme3(param):
    r"""
RIDME experiment model (spin S=3/2)
        
Notes
-----

**Model:**
 
This experiment model has three harmonic pathways and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [\Lambda_0 + \lambda_1 K_0(t,r) + \lambda_2 K_0(2t,r) + \lambda_3 K_0(3t,r)]B(t)


============== =================== ============= ============ ============ ================================================
 Variable        Symbol            Start Values     Lower        Upper                Description
============== =================== ============= ============ ============ ================================================
``param[0]``   :math:`\Lambda_0`     0.3           0            1          Amplitude of unmodulated contribution
``param[1]``   :math:`\lambda_1`     0.5           0            1          Amplitude of 1st harmonic pathway
``param[2]``   :math:`\lambda_2`     0.3           0            1          Amplitude of 2nd harmonic pathway
``param[3]``   :math:`\lambda_3`     0.2           0            1          Amplitude of 3rd harmonic pathway
============== =================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param,npar=4) 

    # Dipolar pathways
    lam = param.copy()
    pathways = [[] for _ in lam]
    pathways[0] = [lam[0]]
    pathways[1] = [lam[1], 0, 1]
    pathways[2] = [lam[2], 0, 2]
    pathways[3] = [lam[3], 0, 3]
    return pathways
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st harmonic pathway',
                'Amplitude of 2nd harmonic pathway','Amplitude of 3rd harmonic pathway',
                'Amplitude of 4th harmonic pathway','Amplitude of 5th harmonic pathway'],
units = ['','','','','',''],
start = np.asarray([0.3, 0.5, 0.3, 0.2,0.1,0.05]),
lower = np.asarray([0, 0, 0, 0, 0, 0]),
upper = np.asarray([1, 1, 1, 1, 1, 1]))
@docstring()
def ex_ridme5(param):
    r"""
RIDME experiment model (spin S=5/2)
        
Notes
-----

**Model:**

This experiment model has five harmonic pathways and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [\Lambda_0 + \lambda_1 K_0(t,r) + \lambda_2 K_0(2t,r) + \lambda_3 K_0(3t,r) + \lambda_4 K_0(4t,r) + \lambda_5 K_0(5t,r)]B(t)

============== =================== ============= ============ ============ ================================================
 Variable        Symbol            Start Values     Lower        Upper                Description
============== =================== ============= ============ ============ ================================================
``param[0]``   :math:`\Lambda_0`     0.3           0            1          Amplitude of unmodulated contribution
``param[1]``   :math:`\lambda_1`     0.5           0            1          Amplitude of 1st harmonic pathway
``param[2]``   :math:`\lambda_2`     0.3           0            1          Amplitude of 2nd harmonic pathway
``param[3]``   :math:`\lambda_3`     0.2           0            1          Amplitude of 3rd harmonic pathway
``param[4]``   :math:`\lambda_4`     0.1           0            1          Amplitude of 4th harmonic pathway
``param[5]``   :math:`\lambda_5`     0.05          0            1          Amplitude of 5th harmonic pathway
============== =================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param, npar=6) 

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
# ===================================================================


# =================================================================
@metadata(
parameters = ['Amplitude of unmodulated contribution','Amplitude of 1st harmonic pathway',
                'Amplitude of 2nd harmonic pathway','Amplitude of 3rd harmonic pathway',
                'Amplitude of 4th harmonic pathway','Amplitude of 5th harmonic pathway',
                'Amplitude of 6th harmonic pathway','Amplitude of 7th harmonic pathway'],
units = ['','','','','','','',''],
start = np.asarray([0.3, 0.5, 0.3, 0.2,0.1,0.05,0.02,0.01]),
lower = np.asarray([0, 0, 0, 0, 0, 0, 0, 0]),
upper = np.asarray([1, 1, 1, 1, 1, 1, 1, 1]))
@docstring()
def ex_ridme7(param):
    r"""
RIDME experiment model (spin S=7/2)

Notes
-----

**Model:**

This experiment model has seven harmonic pathways and an unmodulated contribution. The kernel is 

.. math::
   K(t,r) =
   [\Lambda_0 + \lambda_1 K_0(t,r) + \lambda_2 K_0(2t,r) + \lambda_3 K_0(3t,r) + \lambda_4 K_0(4t,r) + \lambda_5 K_0(5t,r) + \lambda_6 K_0(6t,r) + \lambda_7 K_0(7t,r)]B(t)


============== =================== ============= ============ ============ ================================================
 Variable        Symbol            Start Values     Lower        Upper                Description
============== =================== ============= ============ ============ ================================================
``param[0]``   :math:`\Lambda_0`     0.3           0            1          Amplitude of unmodulated contribution
``param[1]``   :math:`\lambda_1`     0.5           0            1          Amplitude of 1st harmonic pathway
``param[2]``   :math:`\lambda_2`     0.3           0            1          Amplitude of 2nd harmonic pathway
``param[3]``   :math:`\lambda_3`     0.2           0            1          Amplitude of 3rd harmonic pathway
``param[4]``   :math:`\lambda_4`     0.1           0            1          Amplitude of 4th harmonic pathway
``param[5]``   :math:`\lambda_5`     0.05          0            1          Amplitude of 5th harmonic pathway
``param[6]``   :math:`\lambda_6`     0.02          0            1          Amplitude of 6th harmonic pathway
``param[7]``   :math:`\lambda_7`     0.01          0            1          Amplitude of 7th harmonic pathway
============== =================== ============= ============ ============ ================================================
    """  
    param = _parsargs(param,npar=8) 

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
# ===================================================================
.. highlight:: python
.. _ex_4pdeer3:


***********************
:mod:`ex_4pdeer3`
***********************

.. autofunction:: deerlab.ex_models.ex_4pdeer3

Model
=========================================

In order to reduce the parameter space, only the dipolar pathways refocusing at positive times (pathways #1-3) are considered in this model:

.. math::

    K(t,r) =
    [\Lambda_0 + \sum^3_{p=1} \lambda_p K_0(t-T_0^{(p)},r)]
    \prod^3_{p=1} B(t - T_0^{(p)},\lambda_p)

where :math:`T_0^{(1)}=0\;\mu s`, :math:`T_0^{(2)}`, and :math:`T_0^{(3)}` are the refocusing times of the three modulated dipolar pathways at positive evolution times.


============== ======================== ============= ============ ============ ================================================
 Variable        Symbol                  Start Values     Lower        Upper                Description
============== ======================== ============= ============ ============ ================================================
``param[0]``   :math:`\varLambda_0`        0.6                0       1          Unmodulated pathways, amplitude
``param[1]``   :math:`\lambda_1`           0.2                0       1          1st modulated pathway, amplitude
``param[2]``   :math:`\lambda_2`           0.1                0       1          2nd modulated pathway, amplitude
``param[3]``   :math:`T_0^{(2)}`           1.5                0       20         2nd modulated pathway, refocusing time (μs)
``param[4]``   :math:`\lambda_3`           0.1                0       1          3rd modulated pathway, amplitude
``param[5]``   :math:`T_0^{(3)}`           3.5                0       20         3rd modulated pathway, refocusing time (μs)
============== ======================== ============= ============ ============ ================================================


Example
=========================================

Example of a simulated signal using the model evaluated at the start values of its parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.ex_4pdeer3
   t = np.linspace(-5,5,400)
   r = np.linspace(2,5,200)
   info = dl.dd_gauss()
   par0 = info['Start']
   P = dl.dd_gauss(r,par0)
   info = model()
   par0 = info['Start']
   paths = model(par0)
   K = dl.dipolarkernel(t,r,paths,lambda t,lam: dl.bg_hom3d(t,80,lam))
   V = K@P
   plt.figure(figsize=[6,3])
   plt.plot(t,V)
   plt.xlabel('t (µs)',fontsize=13)
   plt.ylabel('V',fontsize=13)
   plt.grid(alpha=0.4)
   plt.tick_params(labelsize=12)
   plt.tick_params(labelsize=12)
   plt.tight_layout()
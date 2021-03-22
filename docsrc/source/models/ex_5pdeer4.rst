.. highlight:: python
.. _ex_5pdeer4:


***********************
:mod:`ex_5pdeer4`
***********************

.. autofunction:: deerlab.ex_models.ex_5pdeer4

Model
=========================================

This experiment model has four modulated pathways and an unmodulated contribution. The kernel is 

.. math::

    K(t,r) =
    [\Lambda_0 + \sum^4_{p=1} \lambda_p K_0(t-T_0^{(p)},r)]
    \prod^4_{p=1} B(t - T_0^{(p)},\lambda_p)

where :math:`T_0^{(1)}=0\;\mu s`, :math:`T_0^{(2)}`, :math:`T_0^{(3)}`, and :math:`T_0^{(4)}` are the refocusing times of the four modulated dipolar pathways at positive evolution times.


============== ======================== ============= ============ ============ ================================================
 Variable        Symbol                  Start Values     Lower        Upper                Description
============== ======================== ============= ============ ============ ================================================
``param[0]``   :math:`\varLambda_0`        0.4                0       1          Unmodulated pathways, amplitude
``param[1]``   :math:`\lambda_1`           0.3                0       1          1st modulated pathway, amplitude
``param[2]``   :math:`\lambda_2`           0.2                0       1          2nd modulated pathway, amplitude
``param[3]``   :math:`T_0^{(2)}`           5.0                0       20         2nd modulated pathway, refocusing time (μs)
``param[4]``   :math:`\lambda_3`           0.05               0       1          3rd modulated pathway, amplitude
``param[5]``   :math:`T_0^{(3)}`           2.0                0       20         3rd modulated pathway, refocusing time (μs)
``param[6]``   :math:`\lambda_4`           0.05               0       1          3rd modulated pathway, amplitude
``param[7]``   :math:`T_0^{(4)}`           3.0                0       20         3rd modulated pathway, refocusing time (μs)
============== ======================== ============= ============ ============ ================================================



Example
=========================================

Example of a simulated signal using the model evaluated at the start values of its parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.ex_5pdeer4
   t = np.linspace(-0.5,9,400)
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
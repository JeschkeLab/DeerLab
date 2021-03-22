.. highlight:: python
.. _ex_5pdeer1:


***********************
:mod:`ex_5pdeer1`
***********************

.. autofunction:: deerlab.ex_models.ex_5pdeer1


Model
=========================================

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



Example
=========================================

Example of a simulated signal using the model evaluated at the start values of its parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.ex_5pdeer1
   t = np.linspace(-0.5,8,400)
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
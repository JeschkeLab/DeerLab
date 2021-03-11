.. highlight:: python
.. _ex_ridme1:


***********************
:mod:`ex_ridme1`
***********************

.. autofunction:: deerlab.ex_models.ex_ridme1


Model
=========================================

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


Example
=========================================

Example of a simulated signal using the model evaluated at the start values of its parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.ex_ridme1
   t = np.linspace(-0.5,5,400)
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
   plt.xlabel('t [μs]',fontsize=13)
   plt.ylabel('V(t)',fontsize=13)
   plt.grid(alpha=0.4)
   plt.tick_params(labelsize=12)
   plt.tick_params(labelsize=12)
   plt.tight_layout()
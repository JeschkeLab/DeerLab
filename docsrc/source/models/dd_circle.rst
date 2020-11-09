.. highlight:: python
.. _dd_circle:


***********************
:mod:`dd_circle`
***********************

.. autofunction:: deerlab.dd_models.dd_circle


Model
=========================================

This provides a `semi-circle distribution <https://en.wikipedia.org/wiki/Wigner_semicircle_distribution>`_, defined by  :math:`P(r) = 2\pi\sqrt{(r-r_0)^2/R^2+1}` for :math:`r_0-R\le r\le r_0+R` and zero otherwise.


============== ================= ============= ============= ============= =================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`r_0`            3.0          0.1              20          Center (nm)
``param[1]``   :math:`R`              0.5          0.1               5          Radius (nm)
============== ================= ============= ============= ============= =================================



Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_circle
   r = np.linspace(2,5,400)
   info = model() 
   par0 = info['Start']
   P = model(r,par0)
   plt.figure(figsize=[6,3])
   plt.plot(r,P)
   plt.xlabel('r [nm]',fontsize=13)
   plt.ylabel('P(r) [nm$^{-1}$]',fontsize=13)
   plt.grid(alpha=0.4)
   plt.tick_params(labelsize=12)
   plt.tick_params(labelsize=12)
   plt.tight_layout()
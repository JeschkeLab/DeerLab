.. highlight:: python
.. _dd_triangle:


***********************
:mod:`dd_triangle`
***********************

.. autofunction:: deerlab.dd_models.dd_triangle

Model
=========================================


This provides a simple triangular distribution.

============== ======================== ============= ============= ============= =========================
 Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =========================
``param[0]``   :math:`r_0`                 3.5            1.0              20         Mode (nm)
``param[1]``   :math:`w_\mathrm{L}`        0.3            0.1              5          Left width (nm)
``param[2]``   :math:`w_\mathrm{R}`        0.3            0.1              5          Right width (nm)
============== ======================== ============= ============= ============= =========================

Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_triangle
   r = np.linspace(2.5,4.5,400)
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
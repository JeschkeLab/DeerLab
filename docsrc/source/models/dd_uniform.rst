.. highlight:: python
.. _dd_uniform:


***********************
:mod:`dd_uniform`
***********************

.. autofunction:: deerlab.dd_models.dd_uniform


Model
=========================================


This provides a simple uniform distribution.

============== ======================== ============= ============= ============= =========================
 Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =========================
``param[0]``   :math:`r_\mathrm{L}`         2.5             0.1            6           Left edge (nm)
``param[1]``   :math:`r_\mathrm{R}`         3.0             0.2            20          Right edge (nm)
============== ======================== ============= ============= ============= =========================

Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_uniform
   r = np.linspace(0.5,4,400)
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


.. highlight:: python
.. _dd_randcoil:

***********************
:mod:`dd_randcoil`
***********************

.. autofunction:: deerlab.dd_models.dd_randcoil


Model
=========================================

.. image:: ../images/model_scheme_dd_randcoil.png
   :width: 35%

:math:`P(r) = \frac{3}{(2\pi\nu_0)^{3/2}}4\pi r^2\exp(-\frac{3 r^2}{\nu_0})`

where :math:`\nu_0 = 3/(12\pi r_0 N \nu)^{3/2}`

============== ============= ============= ============= ============= =======================================
 Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============= ============= ============= ============= =======================================
``param[0]``   :math:`N`         50                2         1000           Number of residues
``param[1]``   :math:`R_0`       0.20           0.10         0.40           Segment length (nm)
``param[2]``   :math:`\nu`       0.60           0.33         1.00           Scaling exponent
============== ============= ============= ============= ============= =======================================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_randcoil
   r = np.linspace(0.5,8,400)
   info = model() 
   par0 = info['Start']
   P = model(r,par0)
   plt.figure(figsize=[6,3])
   plt.plot(r,P)
   plt.xlabel('r (nm)',fontsize=13)
   plt.ylabel('P (nm⁻¹)',fontsize=13)
   plt.grid(alpha=0.4)
   plt.tick_params(labelsize=12)
   plt.tick_params(labelsize=12)
   plt.tight_layout()


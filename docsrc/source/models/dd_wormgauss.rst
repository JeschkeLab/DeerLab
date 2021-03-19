.. highlight:: python
.. _dd_wormgauss:

***********************
:mod:`dd_wormgauss`
***********************


.. autofunction:: deerlab.dd_models.dd_wormgauss


Model
=========================================

============== ============== ============= ============= ============= ==================================
 Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= ==================================
``param[0]``   :math:`L`         3.7            1.5           10         Contour length (nm)
``param[1]``   :math:`L_p`       10              2            100        Persistence length (nm)
``param[2]``   :math:`\sigma`    0.2            0.01           2         Gaussian standard deviation (nm)
============== ============== ============= ============= ============= ==================================

Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_wormgauss
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


References
=========================================

[1] J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
DOI:  `10.1103/PhysRevLett.77.2581 <https://doi.org/10.1103/PhysRevLett.77.2581>`_
.. highlight:: python
.. _dd_gauss:


***********************
:mod:`dd_gauss`
***********************

.. autofunction:: deerlab.dd_models.dd_gauss

Model
=========================================

:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sigma^2}\right)`

============== ======================== ============= ============= ============= ===========================
 Variable       Symbol                   Start value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= ===========================
``param[0]``   :math:`\left<r\right>`     3.5            1.0              20        Mean (nm)
``param[1]``   :math:`\sigma`             0.2            0.05             2.5       Standard deviation (nm)
============== ======================== ============= ============= ============= ===========================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_gauss
   r = np.linspace(2,5,400)
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
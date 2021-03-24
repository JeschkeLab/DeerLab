.. highlight:: python
.. _dd_rice:


***********************
:mod:`dd_rice`
***********************


.. autofunction:: deerlab.dd_models.dd_rice


Model
=========================================

:math:`P(r) = \frac{\nu^{n/2-1}}{\sigma^2}r^{n/2}\exp\left(-\frac{(r^2+\nu^2)}{2\sigma^2}\right)I_{n/2-1}\left(\frac{r\nu}{\sigma^2} \right)`

where :math:`n=3` and :math:`I_{n/2-1}(x)` is the modified Bessel function of the first kind with order :math:`n/2-1`. This is a three-dimensional non-central chi distribution, the 3D generalization of the 2D Rice distribution.

============== ======================== ============= ============= ============= =======================================
 Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =======================================
``param[0]``     :math:`\nu`                3.5           1.0              10          Location (nm)
``param[1]``     :math:`\sigma`             0.7           0.1              5           Spread (nm)
============== ======================== ============= ============= ============= =======================================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_rice
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

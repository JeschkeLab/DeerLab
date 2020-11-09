.. highlight:: python
.. _dd_gauss2:


***********************
:mod:`dd_gauss2`
***********************

.. autofunction:: deerlab.dd_models.dd_gauss2

Model
=========================================

:math:`P(r) = a_1\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_1}\exp\left(-\frac{(r-\left<r_1\right>)^2}{\sigma_1^2}\right) + a_2\sqrt{\frac{2}{\pi}}\frac{1}{\sigma_2}\exp\left(-\frac{(r-\left<r_2\right>)^2}{\sigma_2^2}\right)`


============== ========================= ============= ============= ============= ======================================
 Variable         Symbol                  Start Value   Lower bound   Upper bound      Description
============== ========================= ============= ============= ============= ======================================
``param[0]``   :math:`\left<r_1\right>`     2.5            1.0            20        1st Gaussian mean distance (nm)
``param[1]``   :math:`\sigma_1`             0.2            0.05          2.5        1st Gaussian standard deviation (nm)
``param[2]``   :math:`a_1`                  0.5            0              1         1st Gaussian amplitude
``param[3]``   :math:`\left<r_2\right>`     3.5            1.0            20        2nd Gaussian mean distance (nm)
``param[4]``   :math:`\sigma_2`             0.2            0.05          2.5        2nd Gaussian standard deviation (nm)
``param[5]``   :math:`a_2`                  0.5            0              1         2nd Gaussian amplitude
============== ========================= ============= ============= ============= ======================================

Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_gauss2
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
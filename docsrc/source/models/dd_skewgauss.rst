.. highlight:: python
.. _dd_skewgauss:


***********************
:mod:`dd_skewgauss`
***********************

.. autofunction:: deerlab.dd_models.dd_skewgauss

Model
=========================================

.. image:: ../images/model_scheme_dd_skewgauss.png
   :width: 650px


:math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\sigma}\exp\left(-\frac{(r-\left<r\right>)^2}{\sqrt(2)\sigma^2}\right)\frac{1}{2}\left(1 + erf\left(\frac{(r-\left<r\right>)}{\sqrt{2}\sigma}\right) \right)`

============== ============== ============= ============= ============= =========================
 Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``   :math:`r_0`          3.5         1.0            20         Location (nm)
``param[1]``   :math:`\sigma`       0.5         0.2             5         Spread (nm)
``param[2]``   :math:`\alpha`       5.0         -15            15         Skewness
============== ============== ============= ============= ============= =========================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_skewgauss
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
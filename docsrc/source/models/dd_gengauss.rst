.. highlight:: python
.. _dd_gengauss:


***********************
:mod:`dd_gengauss`
***********************


.. autofunction:: deerlab.dd_models.dd_gengauss

Model
=========================================

.. image:: ../images/model_scheme_dd_gengauss.png
   :width: 650px

:math:`P(r) = \frac{\beta}{2\sigma\Gamma(1/\beta)}\exp\left(-\left(\frac{(r-\left<r\right>)}{\sigma}\right)^\beta \right)`

============== ========================== ============= ============= ============= =======================================
 Variable         Symbol                   Start Value   Lower bound   Upper bound      Description
============== ========================== ============= ============= ============= =======================================
``param[0]``   :math:`r_0`                     3.5          1.0              20         Mean (nm)
``param[1]``   :math:`\sigma`                  0.2          0.05             2.5        Spread (nm)
``param[2]``   :math:`\beta`                   5.0          0.25             15         kurtosis
============== ========================== ============= ============= ============= =======================================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_gengauss
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
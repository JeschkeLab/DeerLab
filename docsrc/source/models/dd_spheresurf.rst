.. highlight:: python
.. _dd_spheresurf:


************************
:mod:`dd_spheresurf`
************************

.. autofunction:: deerlab.dd_models.dd_spheresurf


Model
=========================================

.. image:: ../images/model_scheme_dd_spheresurf.png
   :width: 25%

:math:`P(r) = \begin{cases} \frac{r}{2R^2} \quad \text{for} \quad 0 \leq r < 2R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`


============== ============== ============= ============= ============= =========================
 Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]`     :math:`R`       2.5               0.1          20        Sphere radius (nm)
============== ============== ============= ============= ============= =========================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_spheresurf
   r = np.linspace(0.5,6,400)
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

[1] D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63.
`DOI:  10.1016/j.jmr.2013.01.007 <http://doi.org/10.1016/j.jmr.2013.01.007>`_
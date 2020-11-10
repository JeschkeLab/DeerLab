.. highlight:: python
.. _dd_spherepoint:


************************
:mod:`dd_spherepoint`
************************

.. autofunction:: deerlab.dd_models.dd_spherepoint

Model
=========================================

.. image:: ../images/model_scheme_dd_spherepoint.png
   :width: 25%

:math:`P(r) = \begin{cases} \frac{3r(R^2-(d-r)^2)}{4dR^3} \quad \text{for} \quad d-R \leq r < d+R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`


============== ============== ============= ============= ============= =========================
 Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =========================
``param[0]``     :math:`R`       1.5            0.1            20        Sphere radius (nm)
``param[1]``     :math:`d`       3.5            0.1            20        Distance to point (nm)
============== ============== ============= ============= ============= =========================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_spherepoint
   r = np.linspace(0.5,6,400)
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




References
=========================================

[1] D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63.
`DOI:  10.1016/j.jmr.2013.01.007 <http://doi.org/10.1016/j.jmr.2013.01.007>`_
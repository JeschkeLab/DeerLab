.. highlight:: python
.. _dd_shellshell:


************************
:mod:`dd_shellshell`
************************

.. autofunction:: deerlab.dd_models.dd_shellshell


Model
=========================================

.. image:: ../images/model_scheme_dd_shellshell.png
   :width: 25%

:math:`P(r) = (R_1^3(R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) - R_1^3(R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - R_2^3(R_3^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_3))/((R_3^3 - R_2^3)(R_2^3 - R_1^3))`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w_1`

:math:`R_3 = R + w_1 + w_2`



============== ============== ============= ============= ============= =======================================
 Variable         Symbol       Start Value   Lower bound   Upper bound      Description
============== ============== ============= ============= ============= =======================================
``param[0]``     :math:`R`       1.5             0.1           20         Inner shell radius (nm)
``param[1]``     :math:`w_1`     0.5             0.1           20         1st Shell thickness (nm)
``param[2]``     :math:`w_2`     0.5             0.1           20         2nd Shell thickness (nm)
============== ============== ============= ============= ============= =======================================


Example
=========================================

Example of the model evaluated at the start values of the parameters:

.. plot::

   import deerlab as dl
   import matplotlib.pyplot as plt 
   import numpy as np 
   model = dl.dd_shellshell
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
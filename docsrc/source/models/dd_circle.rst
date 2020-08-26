.. highlight:: python
.. _dd_circle:


***********************
:mod:`dd_circle`
***********************

.. autofunction:: deerlab.dd_models.dd_circle


Model
=========================================

This provides a `semi-circle distribution <https://en.wikipedia.org/wiki/Wigner_semicircle_distribution>`_, defined by  :math:`P(r) = 2\pi\sqrt{(r-r_0)^2/R^2+1}` for :math:`r_0-R\le r\le r_0+R` and zero otherwise.


============== ================= ============= ============= ============= =================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`r_0`            3.0          0.1              20          Center (nm)
``param[1]``   :math:`R`              0.5          0.1               5          Radius (nm)
============== ================= ============= ============= ============= =================================


Example using the start values:

.. image:: ../images/model_dd_circle.png
   :width: 650px

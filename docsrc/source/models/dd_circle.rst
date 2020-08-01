.. highlight:: python
.. _dd_circle:


***********************
:mod:`dd_circle`
***********************

.. autofunction:: deerlab.dd_models.dd_circle


Model
=========================================

This provides a `semi-circle distribution <https://en.wikipedia.org/wiki/Wigner_semicircle_distribution>`_, defined by  :math:`P(r) = 2\pi\sqrt{(r-r_0)^2/R^2+1}` for :math:`r_0-R\le r\le r_0+R` and zero otherwise.


============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_0`                 3.0       0.1              20          center
``param(2)``   :math:`R`                   0.5       0.1               5          radius
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_circle.png
   :width: 650px

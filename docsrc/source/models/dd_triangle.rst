.. highlight:: python
.. _dd_triangle:


***********************
:mod:`dd_triangle`
***********************

.. autofunction:: deerlab.dd_models.dd_triangle

Model
=========================================


This provides a simple triangular distribution.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_0`                3.5       1.0              20         mode
``param(2)``   :math:`w_\mathrm{L}`       0.3       0.1              5          left width
``param(3)``   :math:`w_\mathrm{R}`       0.3       0.1              5          right width
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_triangle.png
   :width: 650px

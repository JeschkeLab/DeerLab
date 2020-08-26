.. highlight:: python
.. _dd_triangle:


***********************
:mod:`dd_triangle`
***********************

.. autofunction:: deerlab.dd_models.dd_triangle

Model
=========================================


This provides a simple triangular distribution.

============== ======================== ============= ============= ============= =========================
 Variable         Symbol                 Start Value   Lower bound   Upper bound      Description
============== ======================== ============= ============= ============= =========================
``param[0]``   :math:`r_0`                 3.5            1.0              20         Mode (nm)
``param[1]``   :math:`w_\mathrm{L}`        0.3            0.1              5          Left width (nm)
``param[2]``   :math:`w_\mathrm{R}`        0.3            0.1              5          Right width (nm)
============== ======================== ============= ============= ============= =========================


Example using start values:

.. image:: ../images/model_dd_triangle.png
   :width: 650px

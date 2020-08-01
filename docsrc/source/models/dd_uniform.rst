.. highlight:: python
.. _dd_uniform:


***********************
:mod:`dd_uniform`
***********************

.. autofunction:: deerlab.dd_models.dd_uniform


Model
=========================================


This provides a simple uniform distribution.

============== ======================== ========= ============= ============= ========================
 Variable       Symbol                    Default   Lower bound   Upper bound      Description
============== ======================== ========= ============= ============= ========================
``param(1)``   :math:`r_\mathrm{L}`         2.5       0.1              6           left edge
``param(2)``   :math:`r_\mathrm{R}`         3.0       0.2              20          right edge
============== ======================== ========= ============= ============= ========================


Example using default parameters:

.. image:: ../images/model_dd_uniform.png
   :width: 650px


.. highlight:: python
.. _bg_homfractal:

***********************
:mod:`bg_homfractal`
***********************

.. autofunction:: deerlab.bg_models.bg_homfractal


Model
=========================================

This implements the background due to a homogeneous distribution of spins in a d-dimensional space, with d-dimensional spin concentration ``c_d``.


============= ============= ========= ============= ============= ==========================================================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ==========================================================
``param(1)``   :math:`c_d`     50          0.01          5000          Pumped spin fractal concentration (μmol/dm\ :sup:`d`)
``param(2)``   :math:`d`       3           0                6          Fractal dimension
============= ============= ========= ============= ============= ==========================================================

.. highlight:: python
.. _bg_poly2:

***********************
:mod:`bg_poly2`
***********************

.. autofunction:: deerlab.bg_models.bg_poly2


Model
=========================================

:math:`B(t) = p_0 + p_1(\lambda t) + p_2(\lambda t)^2`

============= ============= ========= ============= ============= ==============================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ==============================
``param(1)``   :math:`p_0`     1          0            200          Intercept
``param(2)``   :math:`p_1`     -1         -200         200          1st order weight
``param(3)``   :math:`p_2`     -1         -200         200          2nd order weight
============= ============= ========= ============= ============= ==============================

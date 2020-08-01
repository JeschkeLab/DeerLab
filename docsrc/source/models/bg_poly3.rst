.. highlight:: python
.. _bg_poly3:

***********************
:mod:`bg_poly3`
***********************

.. autofunction:: deerlab.bg_models.bg_poly3


Model
=========================================

:math:`B(t) = p_0 + p_1(\lambda t) + p_2(\lambda t)^2 + p_3(\lambda t)^3`

============= ============= ========= ============= ============= ==============================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ==============================
``param(1)``   :math:`p_0`     1          0            200          Intercept
``param(2)``   :math:`p_1`     -1         -200         200          1st order weight
``param(3)``   :math:`p_2`     -1         -200         200          2nd order weight
``param(3)``   :math:`p_3`     -1         -200         200          3rd order weight
============= ============= ========= ============= ============= ==============================

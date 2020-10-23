.. highlight:: python
.. _bg_sumstrexp:

***********************
:mod:`bg_sumstrexp`
***********************

.. autofunction:: deerlab.bg_models.bg_sumstrexp

Model
=========================================

:math:`B(t) = A_1\exp \left(-\kappa_1 \vert t \vert^{d_1}\right) + (1-A_1)\exp\left(-\kappa_2 \vert t \vert^{d_2}\right)`

============== ================= ============= ============= ============= ========================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= ========================================
``param[0]``   :math:`\kappa_1`     3.5             0            200         1st strexp decay rate (μs\ :sup:`-d`)
``param[1]``   :math:`d_1`          1               0            6           1st strexp stretch factor
``param[2]``   :math:`\kappa_2`     3.5             0            200         2nd strexp decay rate (μs\ :sup:`-d`)
``param[3]``   :math:`d_2`          1               0            6           2nd strexp stretch factor
``param[4]``   :math:`A_1`          0.5             0            1           Relative amplitude
============== ================= ============= ============= ============= ========================================
.. highlight:: python
.. _bg_prodstrexp:


***********************
:mod:`bg_prodstrexp`
***********************

.. autofunction:: deerlab.bg_models.bg_prodstrexp

Model
=========================================

:math:`B(t) = \exp\left(-\lambda\kappa_1 \vert t \vert^{d_1}\right) \exp\left(-\lambda\kappa_2 \vert t\vert^{d_2}\right)`

============= ================== ========= ============= ============= ==============================
 Variable         Symbol          Default   Lower bound   Upper bound      Description
============= ================== ========= ============= ============= ==============================
``param(1)``   :math:`\kappa_1`    3.5         0            200         1st strexp decay rate
``param(2)``   :math:`d_1`         1           0            6           1st strexp stretch factor
``param(3)``   :math:`\kappa_2`    3.5         0            200         2nd strexp decay rate
``param(4)``   :math:`d_2`         1           0            6           2nd strexp stretch factor
============= ================== ========= ============= ============= ==============================

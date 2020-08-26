.. highlight:: python
.. _bg_strexp:


***********************
:mod:`bg_strexp`
***********************

.. autofunction:: deerlab.bg_models.bg_strexp

Model
=========================================

.. math::

    B(t) = \exp\left(-\lambda \kappa \vert t\vert^{d}\right)

============== ================= ============= ============= ============= =================================
 Variable         Symbol          Start Value   Lower bound   Upper bound      Description
============== ================= ============= ============= ============= =================================
``param[0]``   :math:`\kappa`      3.5               0              200           decay rate (1/us)
``param[1]``   :math:`d`           1                 0              6             stretch factor
============== ================= ============= ============= ============= =================================

Although the ``bg_strexp`` model has the same functional form as ``bg_homfractal``, it is distinct since its first parameter is a decay rate constant and not a spin concentration like for ``bg_homfractal``.

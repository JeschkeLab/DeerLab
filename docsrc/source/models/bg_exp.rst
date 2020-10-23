.. highlight:: python
.. _bg_exp:

***********************
:mod:`bg_exp`
***********************

.. autofunction:: deerlab.bg_models.bg_exp

Model
=====

.. math::

   B(t) = \exp\left(-\kappa \vert t \vert\right)

============== =============== ============= ============= ============= ================================
 Variable         Symbol        Start Value   Lower bound   Upper bound      Description
============== =============== ============= ============= ============= ================================
``param[0]``   :math:`\kappa`      0.35         0            200          Decay rate (μs\ :sup:`-1`)
============== =============== ============= ============= ============= ================================


Although the ``bg_exp`` model has the same functional form as ``bg_hom3d``, it is distinct since its parameter is a decay rate constant and not a spin concentration like for ``bg_hom3d``. 


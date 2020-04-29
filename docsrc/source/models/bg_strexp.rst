.. highlight:: matlab
.. _bg_strexp:


***********************
:mod:`bg_strexp`
***********************

Stretched-exponential background model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_strexp()
        P = bg_strexp(r,param)
        P = bg_strexp(r,param,lambda)

Inputs
    *   ``t`` -- Time axis (N-array)
    *   ``param`` -- Model parameters
    *   ``lambda`` -- Modulation amplitude (between 0 and 1)

Outputs
    *   ``B`` -- Model background (N-array)
    *   ``info`` -- Model information (struct)


-----------------------------

Model
=========================================

:math:`B(t) = \exp\left(-\lambda \kappa \vert t\vert^{d}\right)`

============= ================= ========= ============= ============= ========================
 Variable       Symbol            Default   Lower bound   Upper bound      Description
============= ================= ========= ============= ============= ========================
``param(1)``   :math:`\kappa`      3.5      0              200           Decay rate
``param(2)``   :math:`d`           1        0              6             Stretch factor
============= ================= ========= ============= ============= ========================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_strexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_strexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_strexp(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda``.

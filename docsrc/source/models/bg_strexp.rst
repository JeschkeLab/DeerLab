.. highlight:: matlab
.. _bg_strexp:


***********************
:mod:`bg_strexp`
***********************

Stretched exponential background parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_strexp()
        P = bg_strexp(r,param)
        P = bg_strexp(r,param,lambda)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
    *   ``lambda`` -Dipolar pathway amplitude

Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`B(t) = \exp\left(-\lambda\kappa_d\vert t\vert^{d/3}\right)`

============= ================= ========= ============= ============= ========================
 Variable       Symbol            Default   Lower bound   Upper bound      Description
============= ================= ========= ============= ============= ========================
``param(1)``   :math:`\kappa_d`    3.5      0              200           Decay rate
``param(2)``   :math:`d`           3        0              6             Fractal dimension
============= ================= ========= ============= ============= ========================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_strexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_strexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a dipolar pathway amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_strexp(t,param,lambda)

Computes the background model ``B`` for a given dipolar pathway amplitude ``lambda``.

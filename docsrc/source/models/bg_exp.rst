.. highlight:: matlab
.. _bg_exp:

***********************
:mod:`bg_exp`
***********************

Product of two stretched exponentials background parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_exp()
        P = bg_exp(r,param)
        P = bg_exp(r,param,lambda)

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

:math:`B(t) = \exp\left(-\lambda\kappa \vert t \vert\right)`

============== =============== ========= ============= ============= ==============================
 Variable         Symbol        Default   Lower bound   Upper bound      Description
============== =============== ========= ============= ============= ==============================
``param(1)``   :math:`\kappa`   0.35         0            200          Decay rate
============== =============== ========= ============= ============= ==============================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_exp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_exp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a dipolar pathway amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_exp(t,param,lambda)

Computes the background model ``B`` for a given dipolar pathway amplitude ``lambda``.


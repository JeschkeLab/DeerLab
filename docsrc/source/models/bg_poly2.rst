.. highlight:: matlab
.. _bg_poly2:

***********************
:mod:`bg_poly2`
***********************

Second order polynomial parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_poly2()
        P = bg_poly2(r,param)
        P = bg_poly2(r,param,lambda)

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

:math:`B(t) = p_0 + p_1(\lambda t) + p_2(\lambda t)^2`

============= ============= ========= ============= ============= ==============================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ==============================
``param(1)``   :math:`p_0`     1          0            200          Intercept
``param(2)``   :math:`p_1`     -1         -200         200          1st order weight
``param(3)``   :math:`p_2`     -1         -200         200          2nd order weight
============= ============= ========= ============= ============= ==============================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_poly2()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_poly2(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a dipolar pathway amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_poly2(t,param,lambda)

Computes the background model ``B`` for a given dipolar pathway amplitude ``lambda``.
.. highlight:: matlab
.. _bg_poly3:

***********************
:mod:`bg_poly3`
***********************

Third-order polynomial (cubic) background model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_poly3()
        P = bg_poly3(r,param)
        P = bg_poly3(r,param,lambda)

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

:math:`B(t) = p_0 + p_1(\lambda t) + p_2(\lambda t)^2 + p_3(\lambda t)^3`

============= ============= ========= ============= ============= ==============================
 Variable       Symbol        Default   Lower bound   Upper bound      Description
============= ============= ========= ============= ============= ==============================
``param(1)``   :math:`p_0`     1          0            200          Intercept
``param(2)``   :math:`p_1`     -1         -200         200          1st order weight
``param(3)``   :math:`p_2`     -1         -200         200          2nd order weight
``param(3)``   :math:`p_3`     -1         -200         200          3rd order weight
============= ============= ========= ============= ============= ==============================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_poly3()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -- Full name of the parametric model.
* ``info.nparam`` -- Total number of adjustable parameters.
* ``info.parameters`` -- Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_poly3(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param`` for a modulation amplitude ``lambda=1``. The required parameters can also be found in the ``info`` structure.

-----------------------------

.. code-block:: matlab

    B = bg_poly3(t,param,lambda)

Computes the background model ``B`` for a given modulation amplitude ``lambda``.
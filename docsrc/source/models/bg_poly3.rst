.. highlight:: matlab
.. _bg_poly3:

***********************
:mod:`bg_poly3`
***********************

Third order polynomial parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_poly3()
        P = bg_poly3(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

:math:`B(t) = p_0 + p_1t + p_2t^2 + p_3t^3`

========== ============= ========= ============= ============= ==============================
 Variable   Symbol        Default   Lower bound   Upper bound      Description
========== ============= ========= ============= ============= ==============================
param(1)    :math:`p_0`     1          0            200          Intercept
param(2)    :math:`p_1`     -1         -200         200          1st order weight
param(3)    :math:`p_2`     -1         -200         200          2nd order weight
param(3)    :math:`p_3`     -1         -200         200          3rd order weight
========== ============= ========= ============= ============= ==============================

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

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


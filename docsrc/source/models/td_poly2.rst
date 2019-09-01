.. highlight:: matlab
.. _td_poly1:

***********************
:mod:`td_poly2`
***********************

Second order polynomial parametric model

Syntax
=========================================

.. code-block:: matlab

        info = td_poly2()
        P = td_poly2(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

Model equation: :math:`B(t) = p_0 + p_1t + p_2t^2`

========== ============= ========= ============= ============= ==============================
 Variable   Symbol        Default   Lower bound   Upper bound      Description
========== ============= ========= ============= ============= ==============================
param(1)    :math:`p_0`     1          0            200          Intercept
param(2)    :math:`p_1`     -1         -200         200          1st order weight
param(3)    :math:`p_2`     -1         -200         200          2nd order weight
========== ============= ========= ============= ============= ==============================

Description
=========================================

.. code-block:: matlab

        info = td_exp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    B = td_exp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


.. highlight:: matlab
.. _td_exp:

***********************
:mod:`td_exp`
***********************

Product of two stretched exponentials background parametric model

Syntax
=========================================

.. code-block:: matlab

        info = td_exp()
        P = td_exp(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

Model equation: :math:`B(t) = e^{-kt}`

========== ============= ========= ============= ============= ==============================
 Variable   Symbol        Default   Lower bound   Upper bound      Description
========== ============= ========= ============= ============= ==============================
param(1)    :math:`k`       3.5         0            200          Decay rate
========== ============= ========= ============= ============= ==============================

Description
=========================================

.. code-block:: matlab

        info = td_exp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    B = td_exp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


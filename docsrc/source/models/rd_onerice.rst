.. highlight:: matlab
.. _rd_onerice:


***********************
:mod:`rd_onerice`
***********************

Rician distribution parametric model

Syntax
=========================================

.. code-block:: matlab

        info = rd_onerice()
        P = rd_onerice(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

Model equation: :math:`P(r) = \frac{r}{\sigma^2}\exp\left(-\frac{(r^2+\left<r\right>^2)}{2\sigma^2}\right)I_0\left(\frac{r\left<r\right>}{\sigma^2} \right)`

where :math:`I_0(x)` is the modified Bessel function of the first kind with order zero.

========== ======================== ========= ============= ============= ========================
 Variable   Symbol                    Default   Lower bound   Upper bound      Description
========== ======================== ========= ============= ============= ========================
param(1)   :math:`\left<r\right>`     3.5     1.0              10         Mean distance
param(2)   :math:`\sigma`             0.7     0.1              5          Standard deviation
========== ======================== ========= ============= ============= ========================

Description
=========================================

.. code-block:: matlab

        info = rd_onerice()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = rd_onerice(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


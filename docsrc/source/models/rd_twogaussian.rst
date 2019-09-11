.. highlight:: matlab
.. _rd_twogaussian:


***********************
:mod:`rd_twogaussian`
***********************

Sum of two Gaussian distributions parametric model

Syntax
=========================================

.. code-block:: matlab

        info = rd_twogaussian()
        P = rd_twogaussian(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

Model equation: :math:`P(r) = A_1\exp\left(-\frac{(r-\left<r_1\right>)^2}{w_1^2}\right) + (1 - A_1)\exp\left(-\frac{(r-\left<r_2\right>)^2}{w_2^2}\right)`

========== ======================== ========= ======== ========= ===================================
 Variable   Symbol                    Default   Lower    Upper       Description
========== ======================== ========= ======== ========= ===================================
param(1)   :math:`\left<r_1\right>`     2.5     1.0        20         1st Gaussian mean distance
param(2)   :math:`w_1`                  0.5     0.2        5          1st Gaussian FWHM
param(3)   :math:`\left<r_2\right>`     3.5     1.0        20         2nd Gaussian mean distance
param(4)   :math:`w_2`                  0.5     0.2        5          2nd Gaussian FWHM
param(5)   :math:`A_1`                  0.5     0          1          1st Gaussian relative amplitude
========== ======================== ========= ======== ========= ===================================

Description
=========================================

.. code-block:: matlab

        info = rd_twogaussian()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = rd_twogaussian(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


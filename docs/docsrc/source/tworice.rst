.. highlight:: matlab
.. _tworice:

***********************
:mod:`tworice`
***********************

Sum of two rician distributions parametric model

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = tworice(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **r** - Distance axis (N-array)
    *   **param** - Model parameters
Returns
    *   **P** - Model distance distribution (N-array)

Model equation: :math:`P(r) = \frac{r}{\sigma_1^2}\exp\left(-\frac{(r^2+\left<r_1\right>^2)}{2\sigma_1^2}\right)I_0\left(\frac{r\left<r_1\right>}{\sigma_1^2} \right) + (1 - A_1) \frac{r}{\sigma_2^2}\exp\left(-\frac{(r^2+\left<r_2\right>^2)}{2\sigma_2^2}\right)I_0\left(\frac{r\left<r_2\right>}{\sigma_2^2} \right)`

where :math:`I_0(x)` is the modified Bessel function of the first kind with order zero.

========== ======================== ========= ======== ======== ===============================
 Variable   Symbol                    Default   Lower   Upper       Description
========== ======================== ========= ======== ======== ===============================
param(1)   :math:`\left<r_1\right>`     2.5     1.0      10      1st Rician mean distance
param(2)   :math:`\sigma_1`             0.4     0.1      5       1st Rician standard deviation
param(3)   :math:`\left<r_2\right>`     4.0     1.0      10      2nd Rician mean distance
param(4)   :math:`\sigma_2`             0.4     0.1      5       2nd Rician standard deviation
param(5)   :math:`A_1`                  0.5     0        1       1st Rician relative amplitude
========== ======================== ========= ======== ======== ===============================

Usage
=========================================

.. code-block:: matlab

        info = tworice()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = tworice(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


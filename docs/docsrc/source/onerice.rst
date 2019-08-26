.. highlight:: matlab
.. _onerice:


***********************
:mod:`onerice`
***********************

Rician distribution parametric model

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = onerice(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **r** - Distance axis (N-array)
    *   **param** - Model parameters
Returns
    *   **P** - Model distance distribution (N-array)

Model equation: :math:`P(r) = \frac{r}{\sigma^2}\exp\left(-\frac{(r^2+\left<r\right>^2)}{2\sigma^2}\right)I_0\left(\frac{r\left<r\right>}{\sigma^2} \right)`

where :math:`I_0(x)` is the modified Bessel function of the first kind with order zero.

========== ======================== ========= ============= ============= ========================
 Variable   Symbol                    Default   Lower bound   Upper bound      Description
========== ======================== ========= ============= ============= ========================
param(1)   :math:`\left<r\right>`     3.5     1.0              10         Mean distance
param(2)   :math:`\sigma`             0.7     0.1              5          Standard deviation
========== ======================== ========= ============= ============= ========================

Usage
=========================================

.. code-block:: matlab

        info = onerice()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = onerice(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


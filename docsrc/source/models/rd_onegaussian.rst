.. highlight:: matlab
.. _rd_onegaussian:


***********************
:mod:`rd_onegaussian`
***********************

Gaussian distribution parametric model

Syntax
=========================================

.. code-block:: matlab

        info = rd_onegaussian()
        P = rd_onegaussian(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

Model equation: :math:`P(r) = \exp\left(-\frac{(r-\left<r\right>)^2}{(\sqrt{2}\sigma)^2}\right)`

========== ======================== ========= ============= ============= ========================
 Variable   Symbol                    Default   Lower bound   Upper bound      Description
========== ======================== ========= ============= ============= ========================
param(1)   :math:`\left<r\right>`     3.5     1.0              20         Mean distance
param(2)   :math:`\sigma`             0.5     0.2              5          Standard deviation
========== ======================== ========= ============= ============= ========================

Description
=========================================

.. code-block:: matlab

        info = rd_onegaussian()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = rd_onegaussian(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


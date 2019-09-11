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

Model equation: :math:`P(r) = \sqrt{\frac{2}{\pi}}\frac{1}{\Gamma}\exp\left(-\frac{(r-\left<r\right>)^2}{\Gamma^2}\right)`

with :math:`\Gamma = w/\sqrt{2ln(2)}`

========== ======================== ========= ============= ============= ========================
 Variable   Symbol                    Default   Lower bound   Upper bound      Description
========== ======================== ========= ============= ============= ========================
param(1)   :math:`\left<r\right>`     3.5     1.0              20         Mean distance
param(2)   :math:`w`                  0.5     0.2              5          FWHM
========== ======================== ========= ============= ============= ========================

Description
=========================================

.. code-block:: matlab

        info = rd_onegaussian()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = rd_onegaussian(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


.. highlight:: matlab
.. _onegaussian:


***********************
:mod:`onegaussian`
***********************

Gaussian distribution parametric model

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = onegaussian(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **r** - Distance axis (N-array)
    *   **param** - Model parameters
Returns
    *   **P** - Model distance distribution (N-array)

Model equation: :math:`P(r) = \exp\left(-\frac{(r-\left<r\right>)^2}{(\sqrt{2}\sigma)^2}\right)`

========== ======================== ========= ============= ============= ========================
 Variable   Symbol                    Default   Lower bound   Upper bound      Description
========== ======================== ========= ============= ============= ========================
param(1)   :math:`\left<r\right>`     3.5     1.0              20         Mean distance
param(2)   :math:`\sigma`             0.5     0.02             5          Standard deviation
========== ======================== ========= ============= ============= ========================

Usage
=========================================

.. code-block:: matlab

        info = onegaussian()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = onegaussian(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


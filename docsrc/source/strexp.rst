.. highlight:: matlab
.. _strexp:


***********************
:mod:`strexp`
***********************

Stretched exponential background parametric model

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`B = strexp(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **t** - Time axis (N-array)
    *   **param** - Model parameters
Returns
    *   **B** - Model background (N-array)

Model equation: :math:`B(t) = e^{-(kt)^{d/3}}`

========== ========== ========= ============= ============= ========================
 Variable   Symbol     Default   Lower bound   Upper bound      Description
========== ========== ========= ============= ============= ========================
param(1)   :math:`k`      3.5      0              200           Decay rate
param(2)   :math:`d`      3        0              6             Fractal dimension
========== ========== ========= ============= ============= ========================

Usage
=========================================

.. code-block:: matlab

        info = strexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    B = strexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


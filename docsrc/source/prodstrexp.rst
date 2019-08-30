.. highlight:: matlab
.. _prodstrexp:


***********************
:mod:`prodstrexp`
***********************

Product of two stretched exponentials background parametric model

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`B = prodstrexp(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **t** - Time axis (N-array)
    *   **param** - Model parameters
Returns
    *   **B** - Model background (N-array)

Model equation: :math:`B(t) = e^{-(k_1t)^{d_1/3}}e^{-(k_2t)^{d_2/3}}`

========== ============= ========= ============= ============= ==============================
 Variable   Symbol        Default   Lower bound   Upper bound      Description
========== ============= ========= ============= ============= ==============================
param(1)   :math:`k_1`      3.5         0            200         1st strexp decay rate
param(2)   :math:`d_1`      3           0            6           1st strexp fractal dimension
param(3)   :math:`k_2`      3.5         0            200         2nd strexp decay rate
param(4)   :math:`d_2`      3           0            6           2nd strexp fractal dimension
========== ============= ========= ============= ============= ==============================

Usage
=========================================

.. code-block:: matlab

        info = prodstrexp()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    B = prodstrexp(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


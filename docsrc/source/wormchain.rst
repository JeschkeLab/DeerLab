.. highlight:: matlab
.. _wormchain:

***********************
:mod:`wormchain`
***********************

Worm-like chain model near the rigid limit

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = wormchain(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **r** - Distance axis (N-array)
    *   **param** - Model parameters
Returns
    *   **P** - Model distance distribution (N-array)

========== =========== ======== ======== ======== ===============================
 Variable   Symbol     Default   Lower   Upper       Description
========== =========== ======== ======== ======== ===============================
param(1)   :math:`L`      3.7     1.5      10       Length of the worm-like chain
param(2)   :math:`L_p`    10      2        100      Persistence length
========== =========== ======== ======== ======== ===============================

Usage
=========================================

.. code-block:: matlab

        info = wormchain()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = wormchain(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


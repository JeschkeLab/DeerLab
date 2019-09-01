.. highlight:: matlab
.. _rd_wormchain:

***********************
:mod:`rd_wormchain`
***********************

Worm-like chain model near the rigid limit

Syntax
=========================================

.. code-block:: matlab

        info = rd_wormchain()
        P = rd_wormchain(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

========== =========== ======== ======== ======== ===============================
 Variable   Symbol     Default   Lower   Upper       Description
========== =========== ======== ======== ======== ===============================
param(1)   :math:`L`      3.7     1.5      10       Length of the worm-like chain
param(2)   :math:`L_p`    10      2        100      Persistence length
========== =========== ======== ======== ======== ===============================

Description
=========================================

.. code-block:: matlab

        info = rd_wormchain()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = rd_wormchain(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

